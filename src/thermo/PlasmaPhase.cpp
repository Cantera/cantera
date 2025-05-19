//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include <boost/math/special_functions/gamma.hpp>
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"
#include "cantera/numerics/funcs.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/Reaction.h"
#include <boost/polymorphic_pointer_cast.hpp>
#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"

namespace Cantera {

namespace {
    const double gamma = sqrt(2 * ElectronCharge / ElectronMass);
}

PlasmaPhase::PlasmaPhase(const string& inputFile, const string& id_)
{
    initialize();

    initThermoFile(inputFile, id_);

    // initial grid
    m_electronEnergyLevels = Eigen::ArrayXd::LinSpaced(m_nPoints, 0.0, 1.0);

    // initial electron temperature
    m_electronTemp = temperature();
}

void PlasmaPhase::initialize()
{
    m_ncs = 0;
    m_f0_ok = false;
    m_EN = 0.0;
    m_E = 0.0;
    m_F = 0.0;
    m_ionDegree = 0.0;
}

PlasmaPhase::~PlasmaPhase()
{
    if (shared_ptr<Solution> soln = m_soln.lock()) {
        soln->removeChangedCallback(this);
        soln->kinetics()->removeReactionAddedCallback(this);
    }
    for (size_t k = 0; k < nCollisions(); k++) {
        // remove callback
        m_collisions[k]->removeSetRateCallback(this);
    }
}

void PlasmaPhase::setTemperature(const double temp)
{
    Phase::setTemperature(temp);
    m_kT = Boltzmann * temp / ElectronCharge;
}

void PlasmaPhase::updateElectronEnergyDistribution()
{
    if (m_distributionType == "discretized") {
        throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
            "Invalid for discretized electron energy distribution.");
    } else if (m_distributionType == "isotropic") {
        setIsotropicElectronEnergyDistribution();
    } else if (m_distributionType == "TwoTermApproximation") {
        auto ierr = ptrEEDFSolver->calculateDistributionFunction();
        if (ierr == 0) {
            auto x = ptrEEDFSolver->getGridEdge();
            auto y = ptrEEDFSolver->getEEDFEdge();
            m_nPoints = x.size();
            m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(x.data(), m_nPoints);
            m_electronEnergyDist = Eigen::Map<const Eigen::ArrayXd>(y.data(), m_nPoints);
        } else {
            throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
                "Call to calculateDistributionFunction failed.");
        }
        bool validEEDF = (
            m_electronEnergyDist.size() == m_nPoints &&
            m_electronEnergyDist.allFinite() &&
            m_electronEnergyDist.maxCoeff() > 0.0 &&
            m_electronEnergyDist.sum() > 0.0
        );

        if (validEEDF) {
            updateElectronTemperatureFromEnergyDist();
        } else {
            writelog("Skipping Te update: EEDF is empty, non-finite, or unnormalized.\n");
        }
    }
    updateElectronEnergyDistDifference();
    electronEnergyDistributionChanged();

}

void PlasmaPhase::normalizeElectronEnergyDistribution() {
    Eigen::ArrayXd eps32 = m_electronEnergyLevels.pow(3./2.);
    double norm = 2./3. * numericalQuadrature(m_quadratureMethod,
                                              m_electronEnergyDist, eps32);
    if (norm < 0.0) {
        throw CanteraError("PlasmaPhase::normalizeElectronEnergyDistribution",
                           "The norm is negative. This might be caused by bad "
                           "electron energy distribution");
    }
    m_electronEnergyDist /= norm;
}

void PlasmaPhase::setElectronEnergyDistributionType(const string& type)
{
    if (type == "discretized" ||
        type == "isotropic" ||
        type == "TwoTermApproximation") {
        m_distributionType = type;
    } else {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistributionType",
            "Unknown type for electron energy distribution.");
    }
}

void PlasmaPhase::setIsotropicElectronEnergyDistribution()
{
    m_electronEnergyDist.resize(m_nPoints);
    double x = m_isotropicShapeFactor;
    double gamma1 = boost::math::tgamma(3.0 / 2.0 / x);
    double gamma2 = boost::math::tgamma(5.0 / 2.0 / x);
    double c1 = x * std::pow(gamma2, 1.5) / std::pow(gamma1, 2.5);
    double c2 = std::pow(gamma2 / gamma1, x);
    m_electronEnergyDist =
        c1  / std::pow(meanElectronEnergy(), 1.5) *
        (-c2 * (m_electronEnergyLevels /
        meanElectronEnergy()).pow(x)).exp();
    checkElectronEnergyDistribution();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    m_electronTemp = Te;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setMeanElectronEnergy(double energy) {
    m_electronTemp = 2.0 / 3.0 * energy * ElectronCharge / Boltzmann;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::setElectronEnergyLevels(const double* levels, size_t length)
{
    m_nPoints = length;
    m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(levels, length);
    checkElectronEnergyLevels();
    electronEnergyLevelChanged();
    updateElectronEnergyDistribution();
}

void PlasmaPhase::electronEnergyDistributionChanged()
{
    m_distNum++;
}

void PlasmaPhase::electronEnergyLevelChanged()
{
    // Cross sections are interpolated on the energy levels
    if (m_collisions.size() > 0) {
        updateInterpolatedCrossSections();
    }
    m_levelNum++;
}

void PlasmaPhase::checkElectronEnergyLevels() const
{
    Eigen::ArrayXd h = m_electronEnergyLevels.tail(m_nPoints - 1) -
                       m_electronEnergyLevels.head(m_nPoints - 1);
    if (m_electronEnergyLevels[0] < 0.0 || (h <= 0.0).any()) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
            "Values of electron energy levels need to be positive and "
            "monotonically increasing.");
    }
}

void PlasmaPhase::checkElectronEnergyDistribution() const
{
    Eigen::ArrayXd h = m_electronEnergyLevels.tail(m_nPoints - 1) -
                       m_electronEnergyLevels.head(m_nPoints - 1);
    if ((m_electronEnergyDist < 0.0).any()) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
            "Values of electron energy distribution cannot be negative.");
    }
    if (m_electronEnergyDist[m_nPoints - 1] > 0.01) {
        warn_user("PlasmaPhase::checkElectronEnergyDistribution",
        "The value of the last element of electron energy distribution exceed 0.01. "
        "This indicates that the value of electron energy level is not high enough "
        "to contain the isotropic distribution at mean electron energy of "
        "{} eV", meanElectronEnergy());
    }
}

void PlasmaPhase::setDiscretizedElectronEnergyDist(const double* levels,
                                                const double* dist,
                                                size_t length)
{
    m_distributionType = "discretized";
    m_nPoints = length;
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::ArrayXd>(levels, length);
    m_electronEnergyDist =
        Eigen::Map<const Eigen::ArrayXd>(dist, length);
    checkElectronEnergyLevels();
    if (m_do_normalizeElectronEnergyDist) {
        normalizeElectronEnergyDistribution();
    }
    checkElectronEnergyDistribution();
    updateElectronEnergyDistDifference();
    updateElectronTemperatureFromEnergyDist();
    electronEnergyLevelChanged();
    electronEnergyDistributionChanged();
}

void PlasmaPhase::updateElectronTemperatureFromEnergyDist()
{
    // calculate mean electron energy and electron temperature
    Eigen::ArrayXd eps52 = m_electronEnergyLevels.pow(5./2.);
    double epsilon_m = 2.0 / 5.0 * numericalQuadrature(m_quadratureMethod,
                                                       m_electronEnergyDist, eps52);
    if (epsilon_m < 0.0 && m_quadratureMethod == "simpson") {
        // try trapezoidal method
        epsilon_m = 2.0 / 5.0 * numericalQuadrature(
            "trapezoidal", m_electronEnergyDist, eps52);
    }

    if (epsilon_m < 0.0) {
        throw CanteraError("PlasmaPhase::updateElectronTemperatureFromEnergyDist",
            "The electron energy distribution produces negative electron temperature.");
    }

    m_electronTemp = 2.0 / 3.0 * epsilon_m * ElectronCharge / Boltzmann;
}

void PlasmaPhase::setIsotropicShapeFactor(double x) {
    m_isotropicShapeFactor = x;
    updateElectronEnergyDistribution();
    //setIsotropicElectronEnergyDistribution();
}

void PlasmaPhase::getParameters(AnyMap& phaseNode) const
{
    IdealGasPhase::getParameters(phaseNode);
    AnyMap eedf;
    eedf["type"] = m_distributionType;
    vector<double> levels(m_nPoints);
    Eigen::Map<Eigen::ArrayXd>(levels.data(), m_nPoints) = m_electronEnergyLevels;
    eedf["energy-levels"] = levels;
    if (m_distributionType == "isotropic") {
        eedf["shape-factor"] = m_isotropicShapeFactor;
        eedf["mean-electron-energy"].setQuantity(meanElectronEnergy(), "eV");
    } else if (m_distributionType == "discretized") {
        vector<double> dist(m_nPoints);
        Eigen::Map<Eigen::ArrayXd>(dist.data(), m_nPoints) = m_electronEnergyDist;
        eedf["distribution"] = dist;
        eedf["normalize"] = m_do_normalizeElectronEnergyDist;
    }
    phaseNode["electron-energy-distribution"] = std::move(eedf);
}


void PlasmaPhase::setParameters(const AnyMap& phaseNode, const AnyMap& rootNode)
{
    IdealGasPhase::setParameters(phaseNode, rootNode);
    m_root = rootNode;

    if (phaseNode.hasKey("electron-energy-distribution")) {
        const AnyMap eedf = phaseNode["electron-energy-distribution"].as<AnyMap>();
        m_distributionType = eedf["type"].asString();

        if (m_distributionType == "isotropic") {
            if (eedf.hasKey("shape-factor")) {
                setIsotropicShapeFactor(eedf["shape-factor"].asDouble());
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "isotropic type requires shape-factor key.");
            }
            if (eedf.hasKey("mean-electron-energy")) {
                double energy = eedf.convert("mean-electron-energy", "eV");
                setMeanElectronEnergy(energy);
            } else {
                throw CanteraError("PlasmaPhase::setParameters",
                    "isotropic type requires electron-temperature key.");
            }
            if (eedf.hasKey("energy-levels")) {
                setElectronEnergyLevels(eedf["energy-levels"].asVector<double>().data(),
                                        eedf["energy-levels"].asVector<double>().size());
            }
            setIsotropicElectronEnergyDistribution();
        } 
        
        else if (m_distributionType == "discretized") {
            if (!eedf.hasKey("energy-levels")) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cannot find key energy-levels.");
            }
            if (!eedf.hasKey("distribution")) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "Cannot find key distribution.");
            }
            if (eedf.hasKey("normalize")) {
                enableNormalizeElectronEnergyDist(eedf["normalize"].asBool());
            }
            setDiscretizedElectronEnergyDist(eedf["energy-levels"].asVector<double>().data(),
                                             eedf["distribution"].asVector<double>().data(),
                                             eedf["energy-levels"].asVector<double>().size());
        } 
        
        else if (m_distributionType == "TwoTermApproximation") {
            bool foundCrossSections = false;

            // Check for 'cross-sections' block (Old Format)
            if (rootNode.hasKey("cross-sections")) {
                for (const auto& item : rootNode["cross-sections"].asVector<AnyMap>()) {
                    addElectronCrossSection(newElectronCrossSection(item));
                }
                foundCrossSections = true;
            }

            // Check for 'electron-collision-plasma' reactions (New Format)
            if (rootNode.hasKey("reactions")) {
                for (const auto& reactionItem : rootNode["reactions"].asVector<AnyMap>()) {
                    if (reactionItem.hasKey("type") &&
                        reactionItem["type"].asString() == "electron-collision-plasma" &&
                        reactionItem.hasKey("energy-levels") &&
                        reactionItem.hasKey("cross-sections")) {

                        // Convert reaction into cross-section format
                        AnyMap newCrossSection;
                        newCrossSection["target"] = reactionItem.getString("target", "unknown");
                        newCrossSection["product"] = reactionItem.getString("product", "unknown");
                        newCrossSection["kind"] = reactionItem.getString("kind", "unknown");
                        newCrossSection["energy-levels"] = reactionItem["energy-levels"].asVector<double>();
                        newCrossSection["cross-sections"] = reactionItem["cross-sections"].asVector<double>();

                        addElectronCrossSection(newElectronCrossSection(newCrossSection));
                        foundCrossSections = true;
                    }
                }
            }

           if (rootNode.hasKey("collisions")) {
                for (const auto& collisionItem : rootNode["collisions"].asVector<AnyMap>()) {
                    if (collisionItem.hasKey("type") &&
                        collisionItem["type"].asString() == "electron-collision-plasma" &&
                        collisionItem.hasKey("energy-levels") &&
                        collisionItem.hasKey("cross-sections")) {

                        AnyMap newCrossSection;

                        // Extract target from equation if not explicitly defined
                        std::string equation = collisionItem.getString("equation", "");
                        std::string targetSpecies = "unknown";
                        std::vector<std::string> productSpeciesList;

                        if (!equation.empty()) {
                            size_t arrowPos = equation.find("=>");
                            if (arrowPos != std::string::npos) {
                                // Extract reactants
                                std::string reactantSide = equation.substr(0, arrowPos);
                                std::vector<std::string> reactants;
                                std::stringstream ssReactants(reactantSide);
                                std::string reactant;
                                while (ssReactants >> reactant) {
                                    if (reactant != "+" && reactant != "e") {  // Ignore '+' and 'e'
                                        reactants.push_back(reactant);
                                    }
                                }

                                // First reactant (non-electron) is the target species
                                if (!reactants.empty()) {
                                    targetSpecies = reactants[0];
                                }

                                // Extract products
                                std::string productSide = equation.substr(arrowPos + 2);
                                std::stringstream ssProducts(productSide);
                                std::string product;
                                while (ssProducts >> product) {
                                    if (product != "+" && product != "e") {
                                        productSpeciesList.push_back(product);
                                    }
                                }

                            }
                        }

                        std::string productListStr = "{ ";
                        for (const auto& p : productSpeciesList) {
                            productListStr += p + " ";
                        }
                        productListStr += "}";

                        std::string kind = "excitation"; // Default type
                        if (productSpeciesList.size() == 1 && productSpeciesList[0] == targetSpecies) {
                            kind = "effective";  // Elastic collision (momentum transfer)
                        } else {
                            for (const auto& p : productSpeciesList) {
                                if (p.back() == '+') {
                                    kind = "ionization";
                                    break;
                                }
                                if (p.back() == '-') {
                                    kind = "attachment";
                                    break;
                                }
                            }
                        }

                        // store the correctly identified data
                        newCrossSection["target"] = targetSpecies;
                        newCrossSection["products"] = productSpeciesList;  
                        newCrossSection["product"] = productSpeciesList.empty() ? "unknown" : productSpeciesList[0]; 
                        newCrossSection["kind"] = kind;
                        newCrossSection["energy-levels"] = collisionItem["energy-levels"].asVector<double>();
                        newCrossSection["cross-sections"] = collisionItem["cross-sections"].asVector<double>();

                        // Convert 'energy-levels' and 'cross-sections' into 'data' format
                        std::vector<std::vector<double>> dataPairs;
                        std::vector<double> energyLevels = collisionItem["energy-levels"].asVector<double>();
                        std::vector<double> crossSections = collisionItem["cross-sections"].asVector<double>();

                        if (energyLevels.size() != crossSections.size()) {
                            throw CanteraError("PlasmaPhase::setParameters",
                                                "Mismatch: `energy-levels` and `cross-sections` must have the same length.");
                        }

                        // Properly format 'data' field
                        for (size_t i = 0; i < energyLevels.size(); i++) {
                            dataPairs.push_back({energyLevels[i], crossSections[i]});
                        }

                        newCrossSection["data"] = dataPairs;

                        addElectronCrossSection(newElectronCrossSection(newCrossSection));
                        foundCrossSections = true;
                    }
                }
            }

            // Throw error if no valid cross-section data is found
            if (!foundCrossSections) {
                throw CanteraError("PlasmaPhase::setParameters",
                    "No valid electron collision cross-section data found.");
            }

            // Initialize the Boltzmann Solver
            ptrEEDFSolver = make_unique<EEDFTwoTermApproximation>(*this);

            // Set Energy Grid (Hardcoded Defaults for Now)
            double kTe_max = 60;
            size_t nGridCells = 301;
            m_nPoints = nGridCells + 1;
            ptrEEDFSolver->setLinearGrid(kTe_max, nGridCells);
        }
    }
}


bool PlasmaPhase::addElectronCrossSection(shared_ptr<ElectronCrossSection> ecs)
{
    // ecs->validate();
    m_ecss.push_back(ecs);

    m_energyLevels.push_back(ecs->energyLevel);
    m_crossSections.push_back(ecs->crossSection);

    // shift factor
    if (ecs->kind == "ionization") {
        m_shiftFactor.push_back(2);
    } else {
        m_shiftFactor.push_back(1);
    }

    // scattering-in factor
    if (ecs->kind == "ionization") {
        m_inFactor.push_back(2);
    } else if (ecs->kind == "attachment") {
        m_inFactor.push_back(0);
    } else {
        m_inFactor.push_back(1);
    }

    if (ecs->kind == "effective" || ecs->kind == "elastic") {
        for (size_t k = 0; k < m_ncs; k++) {
            if (target(k) == ecs->target)
                if (kind(k) == "elastic" || kind(k) == "effective") {
                    throw CanteraError("PlasmaPhase::addElectronCrossSection"
                                       "Already contains a data of effective/ELASTIC cross section for '{}'.",
                                       ecs->target);
            }
        }
        m_kElastic.push_back(m_ncs);
    } else {
        m_kInelastic.push_back(m_ncs);
    }

    // add one to number of cross sections
    m_ncs++;

    m_f0_ok = false;

    return true;
}

bool PlasmaPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = IdealGasPhase::addSpecies(spec);
    size_t k = m_kk - 1;

    if ((spec->name == "e" || spec->name == "Electron") ||
        (spec->composition.find("E") != spec->composition.end() &&
         spec->composition.size() == 1 &&
         spec->composition["E"] == 1)) {
        if (m_electronSpeciesIndex == npos) {
            m_electronSpeciesIndex = k;
        } else {
            throw CanteraError("PlasmaPhase::addSpecies",
                               "Cannot add species, {}. "
                               "Only one electron species is allowed.", spec->name);
        }
    }
    return added;
}

void PlasmaPhase::initThermo()
{
    IdealGasPhase::initThermo();

    // Check electron species
    if (m_electronSpeciesIndex == npos) {
        throw CanteraError("PlasmaPhase::initThermo",
                           "No electron species found.");
    }

    // Initialize kinetics
    m_kinetics = newKinetics("bulk");
    m_kinetics->addThermo(shared_from_this());

    vector<shared_ptr<Reaction>> reactions;
    for (AnyMap R : reactionsAnyMapList(*m_kinetics, m_input, m_root)) {
        shared_ptr<Reaction> reaction = newReaction(R, *m_kinetics);

        // Check if this is an 'electron-collision-plasma' reaction
        if (reaction->type() == "electron-collision-plasma") {
            auto rate = std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(reaction->rate());

            // If the reaction already has 'energy-levels' and 'cross-sections', use them
            if (R.hasKey("energy-levels") && R.hasKey("cross-sections")) {
                rate->set_energyLevels(R["energy-levels"].asVector<double>());
                rate->set_crossSections(R["cross-sections"].asVector<double>());
                rate->set_threshold(R.getDouble("threshold", 0.0));
                rate->set_cs_ok();  // Mark as valid cross-section data

            } else {
                // Try to match with preloaded cross-sections in 'm_ecss[]' (old format)
                for (size_t k = 0; k < m_ncs; k++) {
                    if (rate->target() == m_ecss[k]->target &&
                        rate->product() == m_ecss[k]->product &&
                        rate->kind() == m_ecss[k]->kind) {

                        rate->set_crossSections(m_ecss[k]->crossSection);
                        rate->set_energyLevels(m_ecss[k]->energyLevel);
                        rate->set_threshold(m_ecss[k]->threshold);
                        rate->set_cs_ok();
                    }
                }
            }

            // Check if cross-section data exists
            if (!(rate->get_cs_ok())) {   
                throw CanteraError("PlasmaPhase::initThermo",
                                   "Energy levels and cross-sections are undefined for an electron-collision-plasma reaction.");
            }
        }

        reactions.push_back(reaction);
    }

    // Add reactions to the kinetics object
    addReactions(*m_kinetics, reactions);

    // Initialize 'm_collisions' and identify elastic collisions
    m_collisions.clear();
    size_t i = 0;
    for (shared_ptr<Reaction> reaction : reactions) {
        if (reaction->type() == "electron-collision-plasma") {
            m_collisions.push_back(reaction);

            // Check if the reaction is elastic (reactants = products)
            if (reaction->reactants == reaction->products) {
                m_elasticCollisionIndices.push_back(i);
            }
            i++;  // Count only 'electron-collision-plasma' reactions
        }
    }

    // Interpolate cross-sections
    updateInterpolatedCrossSections();
}

void PlasmaPhase::setSolution(std::weak_ptr<Solution> soln) {
    ThermoPhase::setSolution(soln);
    // register callback function to be executed
    // when the thermo or kinetics object changed
    if (shared_ptr<Solution> soln = m_soln.lock()) {
        soln->registerChangedCallback(this, [&]() {
            setCollisions();
        });
    }
}

void PlasmaPhase::setCollisions()
{
    m_collisions.clear();
    m_collisionRates.clear();
    m_targetSpeciesIndices.clear();

    if (shared_ptr<Solution> soln = m_soln.lock()) {
        shared_ptr<Kinetics> kin = soln->kinetics();
        if (!kin) {
            return;
        }

        // add collision from the initial list of reactions
        for (size_t i = 0; i < kin->nReactions(); i++) {
            std::shared_ptr<Reaction> R = kin->reaction(i);
            if (R->rate()->type() != "electron-collision-plasma") {
                continue;
            }
            addCollision(R);
        }

        // register callback when reaction is added later
        // Modifying collision reactions is not supported
        kin->registerReactionAddedCallback(this, [this, kin]() {
            size_t i = kin->nReactions() - 1;
            if (kin->reaction(i)->type() == "electron-collision-plasma") {
                addCollision(kin->reaction(i));
            }
        });
    }
}

void PlasmaPhase::addCollision(std::shared_ptr<Reaction> collision)
{
    size_t i = nCollisions();

    // setup callback to signal updating the cross-section-related
    // parameters
    collision->registerSetRateCallback(this, [this, i, collision]() {
        m_interp_cs_ready[i] = false;
        m_collisionRates[i] =
            std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(collision->rate());
    });

    // Identify target species for electron-collision reactions
    for (const auto& [name, _] : collision->reactants) {
        // Reactants are expected to be electrons and the target species
        if (name != electronSpeciesName()) {
            m_targetSpeciesIndices.emplace_back(speciesIndex(name));
            break;
        }
    }

    m_collisions.emplace_back(collision);
    m_collisionRates.emplace_back(
        std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(collision->rate()));
    m_interp_cs_ready.emplace_back(false);

    // resize parameters
    m_elasticElectronEnergyLossCoefficients.resize(nCollisions());
}

bool PlasmaPhase::updateInterpolatedCrossSection(size_t i)
{
    if (m_interp_cs_ready[i]) {
        return false;
    }
    vector<double> levels(m_nPoints);
    Eigen::Map<Eigen::ArrayXd>(levels.data(), m_nPoints) = m_electronEnergyLevels;
    m_collisionRates[i]->updateInterpolatedCrossSection(levels);
    m_interp_cs_ready[i] = true;
    return true;
}

void PlasmaPhase::updateInterpolatedCrossSections()
{
    for (shared_ptr<Reaction> collision : m_collisions) {
        auto rate = boost::polymorphic_pointer_downcast
            <ElectronCollisionPlasmaRate>(collision->rate());

        vector<double> cs_interp;
        for (double level : m_electronEnergyLevels) {
            cs_interp.push_back(linearInterp(level,
                rate->energyLevels(), rate->crossSections()));
        }

        // Set interpolated cross-section
        rate->setCrossSectionInterpolated(cs_interp);
    }
}

size_t PlasmaPhase::targetSpeciesIndex(shared_ptr<Reaction> R)
{
    if (R->type() != "electron-collision-plasma") {
        throw CanteraError("PlasmaPhase::targetSpeciesIndex",
            "Invalid reaction type. Type electron-collision-plasma is needed.");
    }
    for (const auto& [name, stoich] : R->reactants) {
        if (name != electronSpeciesName()) {
            return speciesIndex(name);
        }
    }
    throw CanteraError("PlasmaPhase::targetSpeciesIndex",
        "No target found. Target cannot be electron.");
}

vector<double> PlasmaPhase::crossSection(shared_ptr<Reaction> reaction)
{
    if (reaction->type() != "electron-collision-plasma") {
        throw CanteraError("PlasmaPhase::crossSection",
            "Invalid reaction type. Type electron-collision-plasma is needed.");
    } else {
        auto rate = boost::polymorphic_pointer_downcast
            <ElectronCollisionPlasmaRate>(reaction->rate());
        std::vector<double> cs_interp;
        for (double level : m_electronEnergyLevels) {
            cs_interp.push_back(linearInterp(level,
                                rate->energyLevels(),
                                rate->crossSections()));
        }
        return cs_interp;
    }
}

void PlasmaPhase::updateElectronEnergyDistDifference()
{
    m_electronEnergyDistDiff.resize(nElectronEnergyLevels());
    // Forward difference for the first point
    m_electronEnergyDistDiff[0] =
        (m_electronEnergyDist[1] - m_electronEnergyDist[0]) /
        (m_electronEnergyLevels[1] - m_electronEnergyLevels[0]);

    // Central difference for the middle points
    for (size_t i = 1; i < m_nPoints - 1; i++) {
        double h1 = m_electronEnergyLevels[i+1] - m_electronEnergyLevels[i];
        double h0 = m_electronEnergyLevels[i] - m_electronEnergyLevels[i-1];
        m_electronEnergyDistDiff[i] = (h0 * h0 * m_electronEnergyDist[i+1] +
                (h1 * h1 - h0 * h0) * m_electronEnergyDist[i] -
                h1 * h1 * m_electronEnergyDist[i-1]) /
                (h1 * h0) / (h1 + h0);
    }

    // Backward difference for the last point
    m_electronEnergyDistDiff[m_nPoints-1] =
        (m_electronEnergyDist[m_nPoints-1] -
        m_electronEnergyDist[m_nPoints-2]) /
        (m_electronEnergyLevels[m_nPoints-1] -
        m_electronEnergyLevels[m_nPoints-2]);
}

void PlasmaPhase::updateElasticElectronEnergyLossCoefficients()
{
    // cache of cross section plus distribution plus energy-level number
    static const int cacheId = m_cache.getId();
    CachedScalar last_stateNum = m_cache.getScalar(cacheId);

    // combine the distribution and energy level number
    int stateNum = m_distNum + m_levelNum;

    vector<bool> interpChanged(m_collisions.size());
    for (size_t i = 0; i < m_collisions.size(); i++) {
        interpChanged[i] = updateInterpolatedCrossSection(i);
    }

    if (last_stateNum.validate(temperature(), stateNum)) {
        // check each cross section, and only update coefficients that
        // the interpolated cross sections change
        for (size_t i = 0; i < m_collisions.size(); i++) {
            if (interpChanged[i]) {
                updateElasticElectronEnergyLossCoefficient(i);
            }
        }
    } else {
        // update every coefficient if distribution, temperature,
        // or energy levels change.
        for (size_t i = 0; i < m_collisions.size(); i++) {
            updateElasticElectronEnergyLossCoefficient(i);
        }
    }
}

void PlasmaPhase::updateElasticElectronEnergyLossCoefficient(size_t i)
{
    // @todo exclude attachment collisions
    size_t k = m_targetSpeciesIndices[i];

    // Map cross sections to Eigen::ArrayXd
    auto cs_array = Eigen::Map<const Eigen::ArrayXd>(
        m_collisionRates[i]->crossSectionInterpolated().data(),
        m_collisionRates[i]->crossSectionInterpolated().size()
    );

    // Mass ratio calculation
    double mass_ratio = ElectronMass / molecularWeight(k) * Avogadro;

    // Calculate the rate using Simpson's rule or trapezoidal rule
    Eigen::ArrayXd f0_plus = m_electronEnergyDist + Boltzmann * temperature() /
                                ElectronCharge * m_electronEnergyDistDiff;
    m_elasticElectronEnergyLossCoefficients[i] = 2.0 * mass_ratio * gamma *
        numericalQuadrature(
            m_quadratureMethod, 1.0 / 3.0 * f0_plus.cwiseProduct(cs_array),
            m_electronEnergyLevels.pow(3.0));
}

double PlasmaPhase::elasticPowerLoss()
{
    updateElasticElectronEnergyLossCoefficients();
    // The elastic power loss includes the contributions from inelastic
    // collisions (inelastic recoil effects).
    double rate = 0.0;
    for (size_t i = 0; i < nCollisions(); i++) {
        rate += concentration(m_targetSpeciesIndices[i]) *
                m_elasticElectronEnergyLossCoefficients[i];
    }

    return Avogadro * Avogadro * ElectronCharge *
        concentration(m_electronSpeciesIndex) * rate;
}

void PlasmaPhase::updateThermo() const
{
    IdealGasPhase::updateThermo();
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tempNow = temperature();
    double electronTempNow = electronTemperature();
    size_t k = m_electronSpeciesIndex;
    // If the electron temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tempNow || cached.state2 != electronTempNow) {
        m_spthermo.update_single(k, electronTemperature(),
                &m_cp0_R[k], &m_h0_RT[k], &m_s0_R[k]);
        cached.state1 = tempNow;
        cached.state2 = electronTempNow;
    }
    // update the species Gibbs functions
    m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
    // update the nDensity array
}

double PlasmaPhase::enthalpy_mole() const {
    double value = IdealGasPhase::enthalpy_mole();
    value += GasConstant * (electronTemperature() - temperature()) *
             moleFraction(m_electronSpeciesIndex) *
             m_h0_RT[m_electronSpeciesIndex];
    return value;
}

void PlasmaPhase::getGibbs_ref(double* g) const
{
    IdealGasPhase::getGibbs_ref(g);
    g[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getStandardVolumes_ref(double* vol) const
{
    IdealGasPhase::getStandardVolumes_ref(vol);
    vol[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEnthalpies(double* hbar) const
{
    IdealGasPhase::getPartialMolarEnthalpies(hbar);
    hbar[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEntropies(double* sbar) const
{
    IdealGasPhase::getPartialMolarEntropies(sbar);
    double logp = log(pressure());
    double logpe = log(electronPressure());
    sbar[m_electronSpeciesIndex] += GasConstant * (logp - logpe);
}

void PlasmaPhase::getPartialMolarIntEnergies(double* ubar) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
    size_t k = m_electronSpeciesIndex;
    ubar[k] = RTe() * (_h[k] - 1.0);
}

void PlasmaPhase::getChemPotentials(double* mu) const
{
    IdealGasPhase::getChemPotentials(mu);
    size_t k = m_electronSpeciesIndex;
    double xx = std::max(SmallNumber, moleFraction(k));
    mu[k] += (RTe() - RT()) * log(xx);
}

void PlasmaPhase::getStandardChemPotentials(double* muStar) const
{
    IdealGasPhase::getStandardChemPotentials(muStar);
    size_t k = m_electronSpeciesIndex;
    muStar[k] -= log(pressure() / refPressure()) * RT();
    muStar[k] += log(electronPressure() / refPressure()) * RTe();
}

void PlasmaPhase::getEntropy_R(double* sr) const
{
    const vector<double>& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            sr[k] -= tmp;
        } else {
            sr[k] -= log(electronPressure() / refPressure());
        }
    }
}

void PlasmaPhase::getGibbs_RT(double* grt) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            grt[k] += tmp;
        } else {
            grt[k] += log(electronPressure() / refPressure());
        }
    }
}

}
