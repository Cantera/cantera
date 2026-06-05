//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include <boost/math/special_functions/gamma.hpp>
#include "cantera/thermo/Species.h"
#include "cantera/base/global.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/numerics/funcs.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/kinetics/Reaction.h"
#include <boost/polymorphic_pointer_cast.hpp>
#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"

namespace Cantera {

namespace {
    const double gamma = sqrt(2 * ElectronCharge / ElectronMass);
}

PlasmaPhase::PlasmaPhase(const string& inputFile, const string& id_)
{
    initThermoFile(inputFile, id_);

    // initial electron temperature
    m_electronTemp = temperature();

    // Initialize the Boltzmann Solver
    m_eedfSolver = make_unique<EEDFTwoTermApproximation>(this);

    // Set Energy Grid (Hardcoded Defaults for Now)
    double kTe_max = 60;
    size_t nGridCells = 301;
    m_nPoints = nGridCells + 1;
    m_eedfSolver->setLinearGrid(kTe_max, nGridCells);
    m_electronEnergyLevels = asVectorXd(m_eedfSolver->getGridEdge());
    m_electronEnergyDist.setZero(m_nPoints);
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

void PlasmaPhase::initThermo()
{
    IdealGasPhase::initThermo();

    // Check if there is an electron species in the phase.
    if (m_electronSpeciesIndex == npos) {
        throw CanteraError("PlasmaPhase::initThermo",
                           "No electron species found.");
    }
}

void PlasmaPhase::updateThermo() const
{
    // Update the heavy species thermodynamic properties
    // before updating the electron species properties.
    IdealGasPhase::updateThermo();
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tempNow = temperature();
    double electronTempNow = electronTemperature();
    size_t k = m_electronSpeciesIndex;
    // If the electron temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tempNow || cached.state2 != electronTempNow) {
        // Evaluate the electron species thermodynamic properties
        // at the electron temperature.
        m_spthermo.update_single(k, electronTemperature(),
                m_cp0_R[k], m_h0_RT[k], m_s0_R[k]);
        cached.state1 = tempNow;
        cached.state2 = electronTempNow;

        // Update the electron Gibbs functions, with the electron temperature.
        m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
    }
}

// ================================================================= //
//           Overridden from IdealGasPhase or ThermoPhase            //
// ================================================================= //

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

void PlasmaPhase::setSolution(std::weak_ptr<Solution> soln) {
    ThermoPhase::setSolution(soln);
    // Register callback function to be executed
    // when the thermo or kinetics object changed.
    if (shared_ptr<Solution> soln = m_soln.lock()) {
        soln->registerChangedCallback(this, [&]() {
            setCollisions();
        });
    }
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
                auto levels = eedf["energy-levels"].asVector<double>();
                setElectronEnergyLevels(levels);
            }
            setIsotropicElectronEnergyDistribution();
        } else if (m_distributionType == "discretized") {
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
            auto levels = eedf["energy-levels"].asVector<double>();
            auto distribution = eedf["distribution"].asVector<double>(levels.size());
            setDiscretizedElectronEnergyDist(levels, distribution);
        }
    }

    if (rootNode.hasKey("electron-collisions")) {
        for (const auto& item : rootNode["electron-collisions"].asVector<AnyMap>()) {
            auto rate = make_shared<ElectronCollisionPlasmaRate>(item);
            Composition reactants, products;
            reactants[item["target"].asString()] = 1;
            reactants[electronSpeciesName()] = 1;
            if (item.hasKey("product")) {
                products[item["product"].asString()] = 1;
            } else {
                products[item["target"].asString()] = 1;
            }
            products[electronSpeciesName()] = 1;
            if (rate->kind() == "ionization") {
                products[electronSpeciesName()] += 1;
            } else if (rate->kind() == "attachment") {
                products[electronSpeciesName()] -= 1;
            }
            auto R = make_shared<Reaction>(reactants, products, rate);
            addCollision(R);
        }
    }
}

// ================================================================= //
//               Electron Energy Distribution Functions              //
// ================================================================= //

void PlasmaPhase::updateElectronEnergyDistribution()
{
    if (m_distributionType == "discretized") {
        throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
            "Invalid for discretized electron energy distribution.");
    } else if (m_distributionType == "isotropic") {
        setIsotropicElectronEnergyDistribution();
    } else if (m_distributionType == "Boltzmann-two-term") {
        auto ierr = m_eedfSolver->calculateDistributionFunction();
        if (ierr == 0) {
            auto y = m_eedfSolver->getEEDFEdge();
            m_electronEnergyDist = Eigen::Map<const Eigen::ArrayXd>(y.data(), m_nPoints);
        } else {
            throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
                "Call to calculateDistributionFunction failed.");
        }
        bool validEEDF = (
            static_cast<size_t>(m_electronEnergyDist.size()) == m_nPoints &&
            m_electronEnergyDist.allFinite() &&
            m_electronEnergyDist.maxCoeff() > 0.0 &&
            m_electronEnergyDist.sum() > 0.0
        );

        if (validEEDF) {
            updateElectronTemperatureFromEnergyDist();
        } else {
            writelog("Skipping Te update: EEDF is empty, non-finite, or unnormalized.\n");
        }
    } else {
        throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
            "Unknown method '{}' for determining EEDF", m_distributionType);
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
        type == "Boltzmann-two-term") {
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
    if (Te < 0.0) {
        throw CanteraError("PlasmaPhase::setElectronTemperature",
            "Electron temperature cannot be negative.");
    }
    m_electronTemp = Te;
    updateElectronEnergyDistribution();
}

void PlasmaPhase::beginEquilibrate()
{
    ThermoPhase::beginEquilibrate();

    if (!m_inEquilibrate) {
        m_inEquilibrate = true;
        // Remember current Te and lock Te -> T for the duration
        m_electronTempEquil = electronTemperature();
        setElectronTemperature(temperature());
    }
}

void PlasmaPhase::endEquilibrate()
{
    if (m_inEquilibrate) {
        // Restore Te to the pre-equilibrate value
        setElectronTemperature(m_electronTempEquil);
        m_inEquilibrate = false;
    }

    ThermoPhase::endEquilibrate();
}

void PlasmaPhase::setMeanElectronEnergy(double energy) {
    setElectronTemperature(2.0 / 3.0 * energy * ElectronCharge / Boltzmann);
}

void PlasmaPhase::setElectronEnergyLevels(span<const double> levels)
{
    m_nPoints = levels.size();
    m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(levels.data(), m_nPoints);
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
    m_levelNum++;
    // Cross sections are interpolated on the energy levels
    if (m_collisions.size() > 0) {
        for (shared_ptr<Reaction> collision : m_collisions) {
            const auto& rate = boost::polymorphic_pointer_downcast
                <ElectronCollisionPlasmaRate>(collision->rate());
            rate->updateInterpolatedCrossSection(asSpan(m_electronEnergyLevels));
        }
    }
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

void PlasmaPhase::setDiscretizedElectronEnergyDist(span<const double> levels,
                                                  span<const double> dist)
{
    m_distributionType = "discretized";
    m_nPoints = levels.size();
    m_electronEnergyLevels = asVectorXd(levels);
    m_electronEnergyDist = asVectorXd(dist);
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
}

void PlasmaPhase::setCollisions()
{
    if (shared_ptr<Solution> soln = m_soln.lock()) {
        shared_ptr<Kinetics> kin = soln->kinetics();
        if (!kin) {
            return;
        }

        // add collision from the initial list of reactions. Only add reactions we
        // haven't seen before
        set<Reaction*> existing;
        for (auto& R : m_collisions) {
            existing.insert(R.get());
        }
        for (size_t i = 0; i < kin->nReactions(); i++) {
            shared_ptr<Reaction> R = kin->reaction(i);
            if (R->rate()->type() != "electron-collision-plasma"
                || existing.count(R.get())) {
                continue;
            }
            addCollision(R);
        }

        // Register callback when reaction is added later.
        // Modifying collision reactions is not supported.
        kin->registerReactionAddedCallback(this, [this, kin]() {
            size_t i = kin->nReactions() - 1;
            if (kin->reaction(i)->type() == "electron-collision-plasma") {
                addCollision(kin->reaction(i));
            }
        });
    }
}

void PlasmaPhase::addCollision(shared_ptr<Reaction> collision)
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
    string target;
    for (const auto& [name, _] : collision->reactants) {
        // Reactants are expected to be electrons and the target species
        if (name != electronSpeciesName()) {
            m_targetSpeciesIndices.emplace_back(speciesIndex(name, true));
            target = name;
            break;
        }
    }
    if (target.empty()) {
        throw CanteraError("PlasmaPhase::addCollision", "Error identifying target for"
            " collision with equation '{}'", collision->equation());
    }

    m_collisions.emplace_back(collision);
    m_collisionRates.emplace_back(
        std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(collision->rate()));
    m_interp_cs_ready.emplace_back(false);

    // resize parameters
    m_elasticElectronEnergyLossCoefficients.resize(nCollisions());
    updateInterpolatedCrossSection(i);

    // Set up data used by Boltzmann solver
    auto& rate = *m_collisionRates.back();
    string kind = m_collisionRates.back()->kind();

    if ((kind == "effective" || kind == "elastic")) {
        for (size_t k = 0; k < m_collisions.size() - 1; k++) {
            if (m_collisions[k]->reactants == collision->reactants &&
                (m_collisionRates[k]->kind() == "elastic" ||
                 m_collisionRates[k]->kind() == "effective") && !collision->duplicate)
            {
                throw CanteraError("PlasmaPhase::addCollision", "Phase already contains"
                    " an effective/elastic cross section for '{}'.", target);
            }
        }
        m_kElastic.push_back(i);
    } else {
        m_kInelastic.push_back(i);
    }

    auto levels = rate.energyLevels();
    m_energyLevels.emplace_back(levels.begin(), levels.end());
    auto sections = rate.crossSections();
    m_crossSections.emplace_back(sections.begin(), sections.end());
    m_eedfSolver->setGridCache();
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
    if (m_electronEnergyDist.size() != m_nPoints
        || m_electronEnergyDistDiff.size() != m_nPoints) {
        throw CanteraError("PlasmaPhase::elasticPowerLoss:",
            "EEDF not initialized");
    }

    updateElasticElectronEnergyLossCoefficients();
    // The elastic power loss includes the contributions from inelastic
    // collisions (inelastic recoil effects).
    double rate = 0.0;
    for (size_t i = 0; i < nCollisions(); i++) {
        rate += concentration(m_targetSpeciesIndices[i]) *
                m_elasticElectronEnergyLossCoefficients[i];
    }
    const double q_elastic = Avogadro * Avogadro * ElectronCharge *
                concentration(m_electronSpeciesIndex) * rate;

    if (!std::isfinite(q_elastic)) {
        throw CanteraError("PlasmaPhase::elasticPowerLoss:",
            "Non-finite elastic power loss");
    }

    return q_elastic;
}

double PlasmaPhase::electronMobility() const
{
    // Only implemented when using the Boltzmann two-term EEDF
    if (m_distributionType == "Boltzmann-two-term") {
        return m_eedfSolver->getElectronMobility();
    } else {
        throw NotImplementedError("PlasmaPhase::electronMobility",
            "Electron mobility is only available for 'Boltzmann-two-term' "
            "electron energy distributions.");
    }
}

// ================================================================= //
//           Molar Thermodynamic Properties of the Solution          //
// ================================================================= //

double PlasmaPhase::enthalpy_mole() const
{
    m_work.resize(m_kk);
    getPartialMolarEnthalpies(m_work);
    double h = 0.0;
    for (size_t k = 0; k < m_kk; ++k) {
        h += moleFraction(k) * m_work[k];
    }
    return h;
}

double PlasmaPhase::intEnergy_mole() const
{
    m_work.resize(m_kk);
    getPartialMolarIntEnergies(m_work);
    double u = 0.0;
    for (size_t k = 0; k < m_kk; ++k) {
        u += moleFraction(k) * m_work[k];
    }
    return u;
}

double PlasmaPhase::entropy_mole() const
{
    m_work.resize(m_kk);
    getPartialMolarEntropies(m_work);
    double s = 0.0;
    for (size_t k = 0; k < m_kk; ++k) {
        s += moleFraction(k) * m_work[k];
    }
    return s;
}

double PlasmaPhase::gibbs_mole() const
{
    m_work.resize(m_kk);
    getChemPotentials(m_work);
    double g = 0.0;
    for (size_t k = 0; k < m_kk; ++k) {
        g += moleFraction(k) * m_work[k];
    }
    return g;
}

// ================================================================= //
//                     Mechanical Equation of State                  //
// ================================================================= //

double PlasmaPhase::meanTemperature() const
{
    double T_g = temperature();
    double T_e = electronTemperature();
    double X_e = moleFraction(m_electronSpeciesIndex);
    return T_g + X_e * (T_e - T_g);
}

double PlasmaPhase::pressure() const {
    return GasConstant * meanTemperature() * density() / meanMolecularWeight();
}


// ================================================================= //
//                Chemical Potentials and Activities                 //
// ================================================================= //

double PlasmaPhase::standardConcentration(size_t k) const
{
    return pressure() / (GasConstant * meanTemperature());
}


// ================================================================= //
//              Partial Molar Properties of the Solution             //
// ================================================================= //

void PlasmaPhase::getChemPotentials(span<double> mu) const
{
    IdealGasPhase::getChemPotentials(mu);
    size_t k = m_electronSpeciesIndex;
    double xx = std::max(SmallNumber, moleFraction(k));
    mu[k] += (RTe() - RT()) * log(xx);
}

void PlasmaPhase::getPartialMolarEnthalpies(span<double> hbar) const
{
    // Since the `updateThermo` is overriden in `PlasmaPhase`,
    // `enthalpy_RT_ref` returns \tilde{h}_k(T_k) / (R * T_k).
    // When calling `IdealGasPhase::getPartialMolarEnthalpies(hbar)`,
    // the `hbar` array is equal to \tilde{h}_k(T_k) * (R * T) / (R * T_k).
    // For all heavy species, T_k == T, so we get \tilde{h}_k(T).
    // For electrons, we need to multiply by T_e/T to get \tilde{h}_k(T_e).
    IdealGasPhase::getPartialMolarEnthalpies(hbar);
    hbar[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEntropies(span<double> sbar) const
{
    // Since the `updateThermo` is overriden in `PlasmaPhase`,
    // `entropy_R_ref` returns s^\text{ref}_k(T_k)/R.
    // When calling `IdealGasPhase::getPartialMolarEntropies(hbar)`,
    // the `sbar` array is equal to s^\text{ref}_k(T_k)*R/R - R ln(X_k P/P^ref).
    // Therefore, there is no need to correct for temperature.
    IdealGasPhase::getPartialMolarEntropies(sbar);
}

void PlasmaPhase::getPartialMolarIntEnergies(span<double> ubar) const
{
    checkArraySize("PlasmaPhase::getPartialMolarIntEnergies", ubar.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
    // Redefine it for the electron species.
    size_t k = m_electronSpeciesIndex;
    ubar[k] = RTe() * (_h[k] - 1.0);
}

void PlasmaPhase::getPartialMolarVolumes(span<double> vbar) const
{
    double vol = RT() / pressure();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] = vol;
    }
    vbar[m_electronSpeciesIndex] = RTe() / pressure();
}

// ================================================================= //
//  Properties of the Standard State of the Species in the Solution  //
// ================================================================= //

void PlasmaPhase::getStandardChemPotentials(span<double> muStar) const
{
    // After calling PlasmaPhase::getGibbs_ref, muStar = mu^\text{ref}_k(T_k)(T_k).
    // mu^\text{ref} is evaluated at T for heavy species and at Te for electrons.
    getGibbs_ref(muStar);

    // Then, we need to add R*T_k*ln(P/Pref) to mu^\text{ref}.
    // .. For heavy species, mu_star = mu^\text{ref}(T) + R*T*ln(P/Pref)
    double tmp = log(pressure() / refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp;
    }
    // .. For electrons, mu_star = mu^\text{ref}(Te) + R*T_e*ln(P/Pref)
    size_t k = m_electronSpeciesIndex;
    muStar[k] -= log(pressure() / refPressure()) * RT();
    muStar[k] += log(pressure() / refPressure()) * RTe();
}

void PlasmaPhase::getStandardVolumes(span<double> vol) const
{
    double tmp = RT() / pressure();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
    vol[m_electronSpeciesIndex] = RTe() / pressure();
}

// ================================================================= //
//       Thermodynamic Values for the Species Reference States       //
// ================================================================= //

void PlasmaPhase::getGibbs_ref(span<double> g) const
{
    // Since the `updateThermo` is overriden in `PlasmaPhase`,
    // `gibbs_RT_ref` returns \mu^\text{ref}_k(T_k) / (R * T_k).
    // When calling `IdealGasPhase::getGibbs_ref(g)`,
    // the `g` array is equal to \mu^\text{ref}_k(T_k) * (R * T) / (R * T_k).
    // For all heavy species, T_k == T, so we get \mu^\text{ref}_k(T).
    // For electrons, we need to multiply by T_e/T to get \mu^\text{ref}_k(T_e).
    IdealGasPhase::getGibbs_ref(g);
    g[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getStandardVolumes_ref(span<double> vol) const
{
    IdealGasPhase::getStandardVolumes_ref(vol);
    vol[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

// ================================================================= //
//                        Setting the State                          //
// ================================================================= //

void PlasmaPhase::setState(const AnyMap& input_state)
{
    AnyMap state = input_state;

    // Set electron temperature first.
    if (state.hasKey("electron-temperature")) {
        state["Te"] = state["electron-temperature"];
    }

    if (state.hasKey("Te")) {
        setElectronTemperature(state.convert("Te", "K"));
    }

    // Remap allowable synonyms for gas temperature after setting electron temperature,
    if (state.hasKey("gas-temperature")) {
        state["T"] = state["gas-temperature"];
    }
    if (state.hasKey("Tg")) {
        state["T"] = state["Tg"];
    }
}

double PlasmaPhase::jouleHeatingPower() const
{
    // sigma = e * n_e * mu_e   [S/m];   q_J = sigma * E^2   [W/m^3]
    const double mu_e = electronMobility();    // m^2 / (V·s)
    if (mu_e <= 0.0) {
        return 0.0;
    }
    const double ne = concentration(m_electronSpeciesIndex) * Avogadro; // m^-3
    if (ne <= 0.0) {
        return 0.0;
    }
    const double E  = electricField(); // V/m
    if (E <= 0.0) {
        return 0.0;
    }
    const double sigma = ElectronCharge * ne * mu_e; // S/m
    return sigma * E * E; // W/m^3
}

double PlasmaPhase::intrinsicHeating()
{
    // Joule heating: sigma * E^2 [W/m^3]
    const double qJ = jouleHeatingPower();
    checkFinite(qJ);

    // Elastic + inelastic recoil power loss [W/m^3]
    double qElastic = elasticPowerLoss();
    checkFinite(qElastic);

    return qJ + qElastic;
}


}
