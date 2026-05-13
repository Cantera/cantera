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

// Initial constructor for PlasmaPhase. Might have generated some crashes so was changed by the implementation below:
// PlasmaPhase::PlasmaPhase(const string& inputFile, const string& id_)
// {
//     initThermoFile(inputFile, id_);

//     // initial electron temperature
//     m_electronTemp = temperature();

//     // Initialize the Boltzmann Solver
//     m_eedfSolver = make_unique<EEDFTwoTermApproximation>(this);

//     // Set Energy Grid (Hardcoded Defaults for Now)
//     double kTe_max = 60;
//     size_t nGridCells = 301;
//     m_nPoints = nGridCells + 1;
//     m_eedfSolver->setLinearGrid(kTe_max, nGridCells);
//     m_electronEnergyLevels = MappedVector(m_eedfSolver->getGridEdge().data(), m_nPoints);
//     m_electronEnergyDist.setZero(m_nPoints);
// }

// New PlasmaPhase constructor
PlasmaPhase::PlasmaPhase(const string& inputFile, const string& id_)
{
    // Initialize electron temperature before setParameters() can trigger EEDF updates.
    m_electronTemp = temperature();

    // The EEDF solver must exist before initThermoFile(), because setParameters()
    // may add electron collisions and addCollision() updates the solver cache.
    m_eedfSolver = make_unique<EEDFTwoTermApproximation>(this);

    // Safe default grid used before / unless the input file overrides it.
    double kTe_max = 100.0;
    size_t nGridCells = 1001;
    m_eedfSolver->setInitialGridParameters(kTe_max, nGridCells);
    m_eedfSolver->setLinearGrid(kTe_max, nGridCells);

    auto levels = m_eedfSolver->getGridEdge();
    m_nPoints = levels.size();
    m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(
        levels.data(), m_nPoints);
    m_electronEnergyDist.resize(m_nPoints);
    m_electronEnergyDist.setZero();

    if (!inputFile.empty()) {
        initThermoFile(inputFile, id_);
    }

    // If no EEDF was supplied by input, initialize a valid default isotropic EEDF.
    if (m_distributionType == "isotropic" &&
        static_cast<size_t>(m_electronEnergyDist.size()) == m_nPoints &&
        m_electronEnergyDist.sum() == 0.0) {
        updateElectronEnergyDistribution();
    }
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
            auto levels = m_eedfSolver->getGridEdge();
            auto y = m_eedfSolver->getEEDFEdge();

            if (levels.size() != y.size()) {
                throw CanteraError("PlasmaPhase::updateElectronEnergyDistribution",
                    "Inconsistent EEDF solver output: grid edge size and EEDF edge size differ.");
            }

            m_nPoints = levels.size();

            m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(
                levels.data(), m_nPoints);

            m_electronEnergyDist = Eigen::Map<const Eigen::ArrayXd>(
                y.data(), m_nPoints);

            electronEnergyLevelChanged();
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

void PlasmaPhase::normalizeElectronEnergyDistribution()
{
    checkElectronEnergyLevels();
    checkElectronEnergyDistribution();

    Eigen::ArrayXd eps32 = m_electronEnergyLevels.pow(3.0 / 2.0);
    double norm = 2.0 / 3.0 * numericalQuadrature(
        m_quadratureMethod, m_electronEnergyDist, eps32);

    if (!std::isfinite(norm) || norm <= 0.0) {
        throw CanteraError("PlasmaPhase::normalizeElectronEnergyDistribution",
            "Electron energy distribution has invalid norm: {}.", norm);
    }

    m_electronEnergyDist /= norm;
    checkElectronEnergyDistribution();
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

void PlasmaPhase::setElectronTemperature(const double Te)
{
    if (!std::isfinite(Te) || Te <= 0.0) {
        throw CanteraError("PlasmaPhase::setElectronTemperature",
            "Electron temperature must be finite and positive.");
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

void PlasmaPhase::setMeanElectronEnergy(double energy)
{
    if (!std::isfinite(energy) || energy <= 0.0) {
        throw CanteraError("PlasmaPhase::setMeanElectronEnergy",
            "Mean electron energy must be finite and positive.");
    }

    m_electronTemp = 2.0 / 3.0 * energy * ElectronCharge / Boltzmann;
    updateElectronEnergyDistribution();
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
        vector<double> energyLevels(m_nPoints);
        MappedVector(energyLevels.data(), m_nPoints) = m_electronEnergyLevels;
        for (shared_ptr<Reaction> collision : m_collisions) {
            const auto& rate = boost::polymorphic_pointer_downcast
                <ElectronCollisionPlasmaRate>(collision->rate());
            rate->updateInterpolatedCrossSection(energyLevels);
        }
    }
}

void PlasmaPhase::checkElectronEnergyLevels() const
{
    if (m_nPoints < 2 ||
        static_cast<size_t>(m_electronEnergyLevels.size()) != m_nPoints) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
            "Electron energy levels must contain at least two points and match m_nPoints.");
    }

    for (size_t i = 0; i < m_nPoints; i++) {
        const double eps = m_electronEnergyLevels[i];
        if (!std::isfinite(eps)) {
            throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
                "Electron energy level {} is non-finite.", i);
        }
        if (eps < 0.0) {
            throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
                "Electron energy levels must be non-negative.");
        }
        if (i > 0 && eps <= m_electronEnergyLevels[i - 1]) {
            throw CanteraError("PlasmaPhase::checkElectronEnergyLevels",
                "Electron energy levels must be strictly increasing.");
        }
    }
}


void PlasmaPhase::checkElectronEnergyDistribution() const
{
    if (m_nPoints < 2 ||
        static_cast<size_t>(m_electronEnergyDist.size()) != m_nPoints ||
        static_cast<size_t>(m_electronEnergyLevels.size()) != m_nPoints) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
            "Electron energy distribution and energy levels must have matching "
            "sizes of at least two points.");
    }

    double sum = 0.0;
    double maxVal = -BigNumber;

    for (size_t i = 0; i < m_nPoints; i++) {
        const double f = m_electronEnergyDist[i];
        if (!std::isfinite(f)) {
            throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
                "Electron energy distribution contains a non-finite value at index {}.",
                i);
        }
        if (f < 0.0) {
            throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
                "Electron energy distribution cannot contain negative values.");
        }
        sum += f;
        maxVal = std::max(maxVal, f);
    }

    if (!(sum > 0.0) || !(maxVal > 0.0)) {
        throw CanteraError("PlasmaPhase::checkElectronEnergyDistribution",
            "Electron energy distribution must contain at least one positive value.");
    }

    if (m_electronEnergyDist[m_nPoints - 1] > 0.01) {
        warn_user("PlasmaPhase::checkElectronEnergyDistribution",
            "The last value of the electron energy distribution exceeds 0.01. "
            "The electron energy grid may be too short for a mean electron energy "
            "of {} eV.", meanElectronEnergy());
    }
}

void PlasmaPhase::setDiscretizedElectronEnergyDist(span<const double> levels,
                                                  span<const double> dist)
{
    if (levels.size() != dist.size()) {
        throw CanteraError("PlasmaPhase::setDiscretizedElectronEnergyDist",
            "Energy levels and electron energy distribution must have the same size. "
            "Got {} levels and {} distribution values.", levels.size(), dist.size());
    }

    if (levels.size() < 2) {
        throw CanteraError("PlasmaPhase::setDiscretizedElectronEnergyDist",
            "A discretized electron energy distribution requires at least two points.");
    }
    m_distributionType = "discretized";
    m_nPoints = levels.size();
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::ArrayXd>(levels.data(), m_nPoints);
    m_electronEnergyDist =
        Eigen::Map<const Eigen::ArrayXd>(dist.data(), m_nPoints);
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

    if (!std::isfinite(epsilon_m) || epsilon_m <= 0.0) {
        throw CanteraError("PlasmaPhase::updateElectronTemperatureFromEnergyDist",
            "The electron energy distribution produces an invalid mean electron energy: {}.",
            epsilon_m);
    }

    m_electronTemp = 2.0 / 3.0 * epsilon_m * ElectronCharge / Boltzmann;
}

void PlasmaPhase::setIsotropicShapeFactor(double x)
{
    if (!std::isfinite(x) || x <= 0.0) {
        throw CanteraError("PlasmaPhase::setIsotropicShapeFactor",
            "The isotropic shape factor must be finite and positive.");
    }

    m_isotropicShapeFactor = x;
    updateElectronEnergyDistribution();
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
        } else if (m_distributionType == "Boltzmann-two-term") {
            m_eedfSolver = make_unique<EEDFTwoTermApproximation>(this);

            if (eedf.hasKey("energy-levels")) {
                // Mode A: user-provided grid edges.
                // In this mode, the grid is considered explicit and fixed by default.
                auto levels = eedf["energy-levels"].asVector<double>();

                m_eedfSolver->setCustomGrid(levels);
                m_eedfSolver->enableGridAdaptation(false);

                m_nPoints = levels.size();

            } else {
                // Mode B: generated initial grid.
                if (!eedf.hasKey("initial_max_energy_level")) {
                    throw CanteraError("PlasmaPhase::setParameters",
                        "Boltzmann-two-term requires either 'energy-levels' or "
                        "'initial_max_energy_level'.");
                }

                if (!eedf.hasKey("initial_number_of_energy_grid_cells")) {
                    throw CanteraError("PlasmaPhase::setParameters",
                        "Boltzmann-two-term requires either 'energy-levels' or "
                        "'initial_number_of_energy_grid_cells'.");
                }

                double initialMaxEnergy = eedf["initial_max_energy_level"].asDouble();
                size_t nGridCells = static_cast<size_t>(
                    eedf["initial_number_of_energy_grid_cells"].asInt());

                if (!std::isfinite(initialMaxEnergy) || initialMaxEnergy <= 0.0) {
                    throw CanteraError("PlasmaPhase::setParameters",
                        "initial_max_energy_level must be finite and greater than zero.");
                }

                if (nGridCells == 0) {
                    throw CanteraError("PlasmaPhase::setParameters",
                        "initial_number_of_energy_grid_cells must be greater than zero.");
                }

                string energyLevelsDistribution = "Linear";
                if (eedf.hasKey("energy_levels_distribution")) {
                    energyLevelsDistribution =
                        eedf["energy_levels_distribution"].asString();
                } else {
                    writelog("No energy_levels_distribution key found in the input file. "
                            "Defaulting to linear grid.\n");
                }

                m_eedfSolver->setGridType(energyLevelsDistribution);
                m_eedfSolver->setInitialGridParameters(initialMaxEnergy, nGridCells);

                if (energyLevelsDistribution == "Linear") {
                    m_eedfSolver->setLinearGrid(initialMaxEnergy, nGridCells);
                } else if (energyLevelsDistribution == "Quadratic") {
                    m_eedfSolver->setQuadraticGrid(initialMaxEnergy, nGridCells);
                } else if (energyLevelsDistribution == "Geometric") {
                    m_eedfSolver->setGeometricGrid(initialMaxEnergy, nGridCells);
                } else {
                    throw CanteraError("PlasmaPhase::setParameters",
                        "energy_levels_distribution should be Linear, Quadratic or Geometric.");
                }

                if (eedf.hasKey("energy_grid_adaptation")) {
                    const AnyMap adapt = eedf["energy_grid_adaptation"].as<AnyMap>();

                    bool enabled = true;
                    if (adapt.hasKey("enabled")) {
                        enabled = adapt["enabled"].asBool();
                    }

                    double minDecayDecades = 8.0;
                    if (adapt.hasKey("min_decay_decades")) {
                        minDecayDecades = adapt["min_decay_decades"].asDouble();
                    }

                    double maxDecayDecades = 14.0;
                    if (adapt.hasKey("max_decay_decades")) {
                        maxDecayDecades = adapt["max_decay_decades"].asDouble();
                    }

                    double updateFactor = 0.25;
                    if (adapt.hasKey("update_factor")) {
                        updateFactor = adapt["update_factor"].asDouble();
                    }

                    size_t maxIterations = 5;
                    if (adapt.hasKey("max_iterations")) {
                        maxIterations = static_cast<size_t>(
                            adapt["max_iterations"].asInt());
                    }

                    m_eedfSolver->setGridAdaptationParameters(
                        enabled, minDecayDecades, maxDecayDecades,
                        updateFactor, maxIterations);
                } else {
                    m_eedfSolver->enableGridAdaptation(false);
                }

                // In PlasmaPhase, m_nPoints is the number of grid edges.
                m_nPoints = nGridCells + 1;
            }

            m_electronEnergyLevels = Eigen::Map<const Eigen::ArrayXd>(
                m_eedfSolver->getGridEdge().data(), m_eedfSolver->getGridEdge().size());

            m_nPoints = m_eedfSolver->getGridEdge().size();

            m_electronEnergyDist.resize(m_nPoints);
            m_electronEnergyDist.setZero();

            checkElectronEnergyLevels();
            electronEnergyLevelChanged();
        }

 
    }

    // writelog("[plasma-collision-debug] PlasmaPhase::setParameters: scanning electron-collisions\n");
    if (rootNode.hasKey("electron-collisions")) {
        for (const auto& item : rootNode["electron-collisions"].asVector<AnyMap>()) {
            if (item.hasKey("name")) {
                string name = item["name"].asString();
                if (name.empty()) {
                    throw InputFileError("setParameters", rootNode, "Empty electron-collision name.");
                }
                if (m_electronCollisionDefinitions.count(name)) {
                    throw InputFileError("setParameters", rootNode, "Duplicate electron-collision name '{}'.", name);
                }
                m_electronCollisionDefinitions[name] = item;
            }
        //     if (item.hasKey("name")) {
        //         writelog("[plasma-collision-debug]   found named electron-collision: '{}'\n",
        //             item["name"].asString());
        //     } else {
        //         writelog("[plasma-collision-debug]   found anonymous electron-collision entry\n");
        //     }

        //     if (item.hasKey("target")) {
        //         writelog("[plasma-collision-debug]     target: '{}'\n",
        //             item["target"].asString());
        //     }

        //     if (item.hasKey("kind")) {
        //         writelog("[plasma-collision-debug]     kind: '{}'\n",
        //             item["kind"].asString());
        //     }
        }
        // writelog("[plasma-collision-debug] PlasmaPhase::setParameters: registered {} named electron-collision definitions\n",
        // m_electronCollisionDefinitions.size());
    }

    // in the case of a wrong combination of the two entry formats, ensure that each collision is only loaded once
    // writelog("[plasma-collision-debug] PlasmaPhase::setParameters: scanning reactions for collision references\n");
    if (rootNode.hasKey("reactions")) {
        for (const auto& rxnNode : rootNode["reactions"].asVector<AnyMap>()) {
            if (rxnNode.hasKey("type")
                && rxnNode["type"].asString() == "electron-collision-plasma"
                && rxnNode.hasKey("collision")) {
                string name = rxnNode["collision"].asString();
                // writelog("[plasma-collision-debug]   reaction references collision '{}'\n", name);

                if (!m_electronCollisionDefinitions.count(name)) {
                    // writelog("[plasma-collision-debug]     reference '{}' found in electron-collisions\n",
                    //     name);
                    throw InputFileError("setParameters", rootNode,
                        "Reaction references unknown electron collision '{}'.", name);
                }
                m_referencedElectronCollisions.insert(name);
            }
        }
    }

    // writelog("[plasma-collision-debug] PlasmaPhase::setParameters: {} electron-collisions are referenced by reactions\n",
    // m_referencedElectronCollisions.size());

    if (rootNode.hasKey("electron-collisions")) {
        for (const auto& item : rootNode["electron-collisions"].asVector<AnyMap>()) {
            if (item.hasKey("name")) {
                string name = item["name"].asString();

                if (m_referencedElectronCollisions.count(name)) {
                    // writelog("[plasma-collision-debug]   skipping electron-collision '{}' because it is referenced by a reaction\n",
                    //     name);
                     continue;
                 }

                // writelog("[plasma-collision-debug]   adding named standalone electron-collision '{}'\n",
                //     name);
            } else {
                // writelog("[plasma-collision-debug]   adding anonymous standalone electron-collision\n");
            }
            bool hasName = item.hasKey("name");

            if (hasName) {
                string name = item["name"].asString();

                // Nouveau format : si la collision est référencée par une réaction,
                // ne pas créer de réaction synthétique ici.
                if (m_referencedElectronCollisions.count(name)) {
                    continue;
                }
            }

            // Ancien format anonyme, ou nouvelle collision nommée non référencée :
            // elle reste une collision EEDF autonome.
            addStandaloneElectronCollision(item);
            
        }
    }
}


void PlasmaPhase::addStandaloneElectronCollision(const AnyMap& item)
{
    auto rate = make_shared<ElectronCollisionPlasmaRate>(item);

    if (!rate->hasCrossSectionData()) {
        throw InputFileError("addStandaloneElectronCollision", item,
            "Cross-sections are required in the reaction if you are using the deprecated input format!");
    }

    Composition reactants;
    Composition products;

    string target = item["target"].asString();

    reactants[target] = 1;
    reactants[electronSpeciesName()] = 1;

    if (item.hasKey("product")) {
        products[item["product"].asString()] = 1;
    } else {
        products[target] = 1;
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

void PlasmaPhase::addCollision(shared_ptr<Reaction> collision)
{
    // writelog("[plasma-collision-debug] Adding electron collision '{}'\n", collision->equation());
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
    // writelog("[plasma-collision-debug]   identified target species: '{}'\n", target);

    // management of the new data file format

    auto ratePtr = std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(collision->rate());
    // writelog("[plasma-collision-debug]   rate before resolution: collisionName='{}', kind='{}', "
    //      "target='{}', product='{}', threshold={}, hasCrossSectionData={}, nLevels={}, nSections={}\n",
    // ratePtr->collisionName(), ratePtr->kind(), ratePtr->target(), ratePtr->product(),
    // ratePtr->threshold(), ratePtr->hasCrossSectionData(),
    // ratePtr->energyLevels().size(), ratePtr->crossSections().size());

    if (!ratePtr) {
        throw CanteraError("addCollision", "The rate is not initialised");
    }

    if (!ratePtr->collisionName().empty() && !ratePtr->hasCrossSectionData()) {
        auto it = m_electronCollisionDefinitions.find(ratePtr->collisionName());
        if (it == m_electronCollisionDefinitions.end()) {
            throw CanteraError("addCollision", "Unknown electron collision '{}'.", ratePtr->collisionName());
        }

    //     writelog("[plasma-collision-debug]   resolving collision reference '{}'\n",
    // ratePtr->collisionName());
        ratePtr->applyCollisionData(it->second);
    //     writelog("[plasma-collision-debug]   rate after resolution: collisionName='{}', kind='{}', "
    //      "target='{}', product='{}', threshold={}, hasCrossSectionData={}, nLevels={}, nSections={}\n",
    // ratePtr->collisionName(), ratePtr->kind(), ratePtr->target(), ratePtr->product(),
    // ratePtr->threshold(), ratePtr->hasCrossSectionData(),
    // ratePtr->energyLevels().size(), ratePtr->crossSections().size());
    }

    if (!ratePtr->hasCrossSectionData()) {
        throw CanteraError("addCollision",
            "ElectronCollisionPlasmaRate requires either inline cross-section data "
            "or a valid 'collision' reference.");
    }

    if (!ratePtr->target().empty() && ratePtr->target() != target) {
        throw CanteraError("PlasmaPhase::addCollision",
            "Electron collision '{}' targets '{}', but reaction '{}' uses target '{}'.",
            ratePtr->collisionName(), ratePtr->target(),
            collision->equation(), target);
    }
    // writelog("[plasma-collision-debug]   target validation OK: reaction target='{}', collision target='{}'\n",
    // target, ratePtr->target());

    string kindFromReaction = inferElectronCollisionKind(collision);
    string kindFromCollision = ratePtr->kind();

    bool compatibleKind = kindFromReaction == kindFromCollision;

    if ((kindFromReaction == "elastic" || kindFromReaction == "effective") &&
        (kindFromCollision == "elastic" || kindFromCollision == "effective")) {
        compatibleKind = true;
    }

    if (!compatibleKind) {
        throw CanteraError("PlasmaPhase::addCollision",
            "Electron collision '{}' has kind '{}', but reaction '{}' is inferred as '{}'.",
            ratePtr->collisionName(), kindFromCollision,
            collision->equation(), kindFromReaction);
    }
    // writelog("[plasma-collision-debug]   kind validation OK: reaction inferred kind='{}', collision kind='{}'\n",
    // kindFromReaction, ratePtr->kind());

    // writelog("[plasma-collision-debug]   adding electron collision at index {}\n", i);
    m_collisions.emplace_back(collision);
    m_collisionRates.emplace_back(
        std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(collision->rate()));
    m_interp_cs_ready.emplace_back(false);

    // resize parameters
    m_elasticElectronEnergyLossCoefficients.resize(nCollisions());
    updateInterpolatedCrossSection(i);

    // Set up data used by Boltzmann solver
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
//     if ((kind == "effective" || kind == "elastic")) {
//     writelog("[plasma-collision-debug]   classified as elastic/effective collision, index={}\n", i);
// } else {
//     writelog("[plasma-collision-debug]   classified as inelastic collision, index={}\n", i);
// }

    auto levels = ratePtr->energyLevels();
    m_energyLevels.emplace_back(levels.begin(), levels.end());
    auto sections = ratePtr->crossSections();
    m_crossSections.emplace_back(sections.begin(), sections.end());
    m_eedfSolver->setGridCache();

    // writelog("[plasma-collision-debug]   addCollision done: nCollisions={}, nElastic={}, nInelastic={}\n",
    // nCollisions(), m_kElastic.size(), m_kInelastic.size());
}

string PlasmaPhase::inferElectronCollisionKind(
    const shared_ptr<Reaction>& collision) const
{
    const string eName = electronSpeciesName();

    double nReactantElectrons = 0.0;
    double nProductElectrons = 0.0;

    auto reactantElectron = collision->reactants.find(eName);
    if (reactantElectron != collision->reactants.end()) {
        nReactantElectrons = reactantElectron->second;
    }

    auto productElectron = collision->products.find(eName);
    if (productElectron != collision->products.end()) {
        nProductElectrons = productElectron->second;
    }

    if (nProductElectrons > nReactantElectrons) {
        return "ionization";
    }

    if (nProductElectrons < nReactantElectrons) {
        return "attachment";
    }

    if (collision->reactants == collision->products) {
        return "elastic";
    }

    return "excitation";
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
    if (m_elasticElectronEnergyLossCoefficients.size() != nCollisions()) {
        m_elasticElectronEnergyLossCoefficients.resize(nCollisions(), 0.0);
    }

    const string kind = m_collisionRates[i]->kind();

    if (kind == "attachment") {
        m_elasticElectronEnergyLossCoefficients[i] = 0.0;
        return;
    }

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
        const string kind = m_collisionRates[i]->kind();

        if (kind == "attachment") {
            continue;
        }

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
                m_cp0_R[k], m_h0_RT[k], m_s0_R[k]);
        cached.state1 = tempNow;
        cached.state2 = electronTempNow;
    }
    // update the species Gibbs functions
    m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
}

double PlasmaPhase::enthalpy_mole() const {
    double value = IdealGasPhase::enthalpy_mole();
    value += GasConstant * (electronTemperature() - temperature()) *
             moleFraction(m_electronSpeciesIndex) *
             m_h0_RT[m_electronSpeciesIndex];
    return value;
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

void PlasmaPhase::getGibbs_ref(span<double> g) const
{
    IdealGasPhase::getGibbs_ref(g);
    g[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getStandardVolumes_ref(span<double> vol) const
{
    IdealGasPhase::getStandardVolumes_ref(vol);
    vol[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEnthalpies(span<double> hbar) const
{
    IdealGasPhase::getPartialMolarEnthalpies(hbar);
    hbar[m_electronSpeciesIndex] *= electronTemperature() / temperature();
}

void PlasmaPhase::getPartialMolarEntropies(span<double> sbar) const
{
    IdealGasPhase::getPartialMolarEntropies(sbar);
    double logp = log(pressure());
    double logpe = log(electronPressure());
    sbar[m_electronSpeciesIndex] += GasConstant * (logp - logpe);
}

void PlasmaPhase::getPartialMolarIntEnergies(span<double> ubar) const
{
    checkArraySize("PlasmaPhase::getPartialMolarIntEnergies", ubar.size(), m_kk);
    auto _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
    size_t k = m_electronSpeciesIndex;
    ubar[k] = RTe() * (_h[k] - 1.0);
}

void PlasmaPhase::getChemPotentials(span<double> mu) const
{
    IdealGasPhase::getChemPotentials(mu);
    size_t k = m_electronSpeciesIndex;
    double xx = std::max(SmallNumber, moleFraction(k));
    mu[k] += (RTe() - RT()) * log(xx);
}

void PlasmaPhase::getStandardChemPotentials(span<double> muStar) const
{
    IdealGasPhase::getStandardChemPotentials(muStar);
    size_t k = m_electronSpeciesIndex;
    muStar[k] -= log(pressure() / refPressure()) * RT();
    muStar[k] += log(electronPressure() / refPressure()) * RTe();
}

void PlasmaPhase::getEntropy_R(span<double> sr) const
{
    checkArraySize("PlasmaPhase::getEntropy_R", sr.size(), m_kk);
    auto _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr.begin());
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            sr[k] -= tmp;
        } else {
            sr[k] -= log(electronPressure() / refPressure());
        }
    }
}

void PlasmaPhase::getGibbs_RT(span<double> grt) const
{
    checkArraySize("PlasmaPhase::getGibbs_RT", grt.size(), m_kk);
    auto gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt.begin());
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_electronSpeciesIndex) {
            grt[k] += tmp;
        } else {
            grt[k] += log(electronPressure() / refPressure());
        }
    }
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

double PlasmaPhase::inelasticPower()
{
    // Joule heating: sigma * E^2 [W/m^3]
    const double qJ = jouleHeatingPower();
    checkFinite(qJ);

    // Elastic + inelastic recoil power loss [W/m^3]
    double qElastic = elasticPowerLoss();
    checkFinite(qElastic);

    return qJ - qElastic;
}

double PlasmaPhase::intrinsicHeating()
{
    // Joule heating: sigma * E^2 [W/m^3]
    const double qJ = jouleHeatingPower();
    checkFinite(qJ);

    return qJ;
}


}
