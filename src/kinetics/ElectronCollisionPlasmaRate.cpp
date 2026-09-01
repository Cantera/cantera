//! @file ElectronCollisionPlasmaRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/numerics/funcs.h"

namespace Cantera
{

ElectronCollisionPlasmaData::ElectronCollisionPlasmaData()
{
    energyLevels.assign(1, 0.0);
    distribution.assign(1, 0.0);
}

bool ElectronCollisionPlasmaData::update(const ThermoPhase& phase, const Kinetics& kin)
{
    const PlasmaPhase& pp = dynamic_cast<const PlasmaPhase&>(phase);

    // The distribution number dictates whether the rate should be updated.
    // Three scenarios involving changes of the distribution number:
    // 1. Change of the electron energy levels
    // 2. Change of the electron energy distribution
    // 3. Combined changes of one and two
    if (pp.distributionNumber() == m_dist_number) {
        return false;
    }

    // Update distribution
    m_dist_number = pp.distributionNumber();
    distribution.resize(pp.nElectronEnergyLevels());
    pp.getElectronEnergyDistribution(distribution);

    // Update energy levels
    if (pp.levelNumber() != levelNumber || energyLevels.empty()) {
        levelNumber = pp.levelNumber();
        energyLevels.resize(pp.nElectronEnergyLevels());
        pp.getElectronEnergyLevels(energyLevels);
    }
    return true;
}

void ElectronCollisionPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);
    
    if (!node.hasKey("collision")) {
        throw InputFileError("ElectronCollisionPlasmaRate::setParameters", node,
            "Electron-collision reactions require a named 'collision' reference. "
            "Tabulated cross-section data must be declared in the root "
            "'electron-collisions' section.");
    }

    m_collisionName = node["collision"].asString();

    if (m_collisionName.empty()) {
        throw InputFileError("ElectronCollisionPlasmaRate::setParameters", node,
            "The 'collision' reference cannot be empty.");
    }

    for (const string& key : {
        "name", "kind", "target", "product", "threshold",
        "energy-levels", "cross-sections"
    }) {
        if (node.hasKey(key)) {
            throw InputFileError("ElectronCollisionPlasmaRate::setParameters", node,
                "Electron-collision reaction entries cannot contain '{}'. "
                "Collision metadata and tabulated data must be declared in the "
                "referenced root 'electron-collisions' entry.",
                key);
        }
    }
}

void ElectronCollisionPlasmaRate::getParameters(AnyMap& node) const {
    node["type"] = type();
    node["collision"] = m_collisionName;
}

void ElectronCollisionPlasmaRate::updateInterpolatedCrossSection(
    span<const double> sharedLevels)
{
    m_crossSectionsInterpolated.clear();
    m_crossSectionsInterpolated.reserve(sharedLevels.size());
    for (double level : sharedLevels) {
        m_crossSectionsInterpolated.emplace_back(linearInterp(level,
                                            m_energyLevels, m_crossSections));
    }
}

double ElectronCollisionPlasmaRate::evalFromStruct(
    const ElectronCollisionPlasmaData& shared_data)
{
    // Interpolate cross-sections data to the energy levels of
    // the electron energy distribution function when the EEDF from the phase changes
    if (m_levelNumber != shared_data.levelNumber) {
        m_crossSectionsInterpolated.clear();
        for (double level : shared_data.energyLevels) {
            m_crossSectionsInterpolated.push_back(linearInterp(level,
                                                  m_energyLevels, m_crossSections));
        }
        m_levelNumber = shared_data.levelNumber;
    }

    AssertThrowMsg(m_crossSectionsInterpolated.size() == shared_data.distribution.size(),
        "ECPR:evalFromStruct", "Size mismatch: len(interp) = {}, len(distrib) = {}",
        m_crossSectionsInterpolated.size(), shared_data.distribution.size());

    // Map cross sections to Eigen::ArrayXd
    auto cs_array = Eigen::Map<const Eigen::ArrayXd>(
        m_crossSectionsInterpolated.data(), m_crossSectionsInterpolated.size()
    );

    // Map energyLevels in Eigen::ArrayXd
    auto eps = Eigen::Map<const Eigen::ArrayXd>(
        shared_data.energyLevels.data(), shared_data.energyLevels.size()
    );

    // Map energyLevels in Eigen::ArrayXd
    auto distribution = Eigen::Map<const Eigen::ArrayXd>(
        shared_data.distribution.data(), shared_data.distribution.size()
    );

    // unit in kmol/m3/s
    return pow(2.0 * ElectronCharge / ElectronMass, 0.5) * Avogadro *
           simpson(eps.cwiseProduct(distribution.cwiseProduct(cs_array)), eps);
}

void ElectronCollisionPlasmaRate::modifyRateConstants(
    const ElectronCollisionPlasmaData& shared_data, double& kf, double& kr)
{
    if (kr == 0.0) {
        // The reverse rate constant is only for reversible reactions
        // kr = 0.0 indicates that the reaction is irreversible
        return;
    }

    // Interpolate cross-sections data to the energy levels of
    // the electron energy distribution function
    if (m_levelNumberSuperelastic != shared_data.levelNumber) {
        // super elastic collision energy levels and cross-sections
        vector<double> superElasticEnergyLevels{0.0};
        m_crossSectionsOffset.resize(shared_data.energyLevels.size());
        for (size_t i = 1; i < m_energyLevels.size(); i++) {
            // The energy levels are offset by the first energy level (threshold)
            superElasticEnergyLevels.push_back(m_energyLevels[i] - m_energyLevels[0]);
        }
        for (size_t i = 0; i < shared_data.energyLevels.size(); i++) {
            // The interpolated super-elastic cross section is evaluated
            // at the shared energy grid
            m_crossSectionsOffset[i] = linearInterp(shared_data.energyLevels[i],
                                                        superElasticEnergyLevels,
                                                        m_crossSections);
        }
        m_levelNumberSuperelastic = shared_data.levelNumber;
    }

    // Map energyLevels in Eigen::ArrayXd
    auto eps = Eigen::Map<const Eigen::ArrayXd>(
        shared_data.energyLevels.data(), shared_data.energyLevels.size()
    );

    // Map energyLevels in Eigen::ArrayXd
    auto distribution = Eigen::Map<const Eigen::ArrayXd>(
        shared_data.distribution.data(), shared_data.distribution.size()
    );

    // unit in kmol/m3/s
    kr = pow(2.0 * ElectronCharge / ElectronMass, 0.5) * Avogadro *
         simpson((eps + m_energyLevels[0]).cwiseProduct(
         distribution.cwiseProduct(m_crossSectionsOffset)), eps);
}

void ElectronCollisionPlasmaRate::setContext(const Reaction& rxn, const Kinetics& kin)
{
    const ThermoPhase& thermo = kin.thermo();
    // get electron species name
    string electronName;
    if (thermo.type() == "plasma") {
        electronName = dynamic_cast<const PlasmaPhase&>(thermo).electronSpeciesName();
    } else {
        throw CanteraError("ElectronCollisionPlasmaRate::setContext",
                           "ElectronCollisionPlasmaRate requires plasma phase");
    }

    // Number of reactants needs to be two
    if (rxn.reactants.size() != 2) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires exactly two reactants");
    }

    // Must have only one electron
    // @todo add electron-electron collision rate
    if (rxn.reactants.at(electronName) != 1) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires one and only one electron");
    }

    // Determine the "kind" of collision if not specified explicitly
    if (m_kind.empty()) {
        m_kind = "excitation"; // default
        if (rxn.reactants == rxn.products) {
            m_kind = "effective";
        } else {
            for (const auto& [p, stoich] : rxn.products) {
                if (p == electronName) {
                    continue;
                }
                double q = thermo.charge(thermo.speciesIndex(p, true));
                if (q > 0) {
                    m_kind = "ionization";
                } else if (q < 0) {
                    m_kind = "attachment";
                }
            }
        }
    }

    if (m_threshold == 0.0 &&
        (m_kind == "excitation" || m_kind == "ionization" || m_kind == "attachment"))
    {
        for (size_t i = 0; i < m_energyLevels.size(); i++) {
            if (m_energyLevels[i] > 0.0) {  // Look for first non-zero cross-section
                m_threshold = m_energyLevels[i];
                break;
            }
        }
    }

    if (!rxn.reversible) {
        return; // end checking of forward reaction
    }

    // For super-elastic collisions
    if (rxn.products.size() != 2) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires exactly two products"
            " if the reaction is reversible (super-elastic collisions)");
    }

    // Must have only one electron
    if (rxn.products.at(electronName) != 1) {
        throw InputFileError("ElectronCollisionPlasmaRate::setContext", rxn.input,
            "ElectronCollisionPlasmaRate requires one and only one electron in products"
            " if the reaction is reversible (super-elastic collisions)");
    }
}

void ElectronCollisionPlasmaRate::applyCollisionData(const AnyMap& node)
{
    const string routineName =
        "ElectronCollisionPlasmaRate::applyCollisionData";

    if (!node.hasKey("name")) {
        throw InputFileError(routineName, node,
            "Electron-collision definitions require a unique 'name'.");
    }

    const string name = node["name"].asString();

    if (name.empty()) {
        throw InputFileError(routineName, node,
            "Electron-collision definition names cannot be empty.");
    }

    if (!m_collisionName.empty() && m_collisionName != name) {
        throw InputFileError(routineName, node,
            "Reaction references electron collision '{}', but data for '{}' "
            "were supplied.",
            m_collisionName, name);
    }

    m_collisionName = name;

    if (!node.hasKey("kind")) {
        throw InputFileError(routineName, node,
            "Electron-collision definition '{}' requires 'kind'.", name);
    }

    if (!node.hasKey("target")) {
        throw InputFileError(routineName, node,
            "Electron-collision definition '{}' requires 'target'.", name);
    }

    if (!node.hasKey("energy-levels")) {
        throw InputFileError(routineName, node,
            "Electron-collision definition '{}' requires 'energy-levels'.",
            name);
    }

    if (!node.hasKey("cross-sections")) {
        throw InputFileError(routineName, node,
            "Electron-collision definition '{}' requires 'cross-sections'.",
            name);
    }

    m_kind = node["kind"].asString();
    m_target = node["target"].asString();
    m_product = node.getString("product", "");

    m_energyLevels = node["energy-levels"].asVector<double>();
    m_crossSections = node["cross-sections"].asVector<double>(m_energyLevels.size());

    m_threshold = node.getDouble("threshold", 0.0);
    setDefaultThreshold();
    validateCollisionData(node);

    // Invalidate all interpolated data after assigning a new table.
    m_levelNumber = -3;
    m_levelNumberSuperelastic = -2;
    m_crossSectionsInterpolated.clear();
    m_crossSectionsOffset.resize(0);

    m_hasCrossSectionData = true;
}

void ElectronCollisionPlasmaRate::validateCollisionData(
    const AnyMap& node) const
{
    const string routineName = "ElectronCollisionPlasmaRate::validateCollisionData";

    static const set<string> validKinds = {
        "effective", "elastic", "excitation", "ionization", "attachment"
    };

    if (!validKinds.count(m_kind)) {
        throw InputFileError(routineName, node,
            "Unknown electron-collision kind '{}'. Expected one of "
            "'effective', 'elastic', 'excitation', 'ionization', or "
            "'attachment'.",
            m_kind);
    }

    if (m_target.empty()) {
        throw InputFileError(routineName, node,
            "The electron-collision target cannot be empty.");
    }

    if (m_energyLevels.size() < 2) {
        throw InputFileError(routineName, node,
            "Electron-collision data require at least two energy levels.");
    }

    if (m_energyLevels.size() != m_crossSections.size()) {
        throw InputFileError(routineName, node,
            "The 'energy-levels' and 'cross-sections' arrays must have "
            "identical lengths.");
    }

    for (size_t i = 0; i < m_energyLevels.size(); i++) {
        if (!std::isfinite(m_energyLevels[i]) || m_energyLevels[i] < 0.0) {
            throw InputFileError(routineName, node,
                "Energy levels must be finite and non-negative.");
        }

        if (!std::isfinite(m_crossSections[i]) ||
            m_crossSections[i] < 0.0) {
            throw InputFileError(routineName, node,
                "Cross sections must be finite and non-negative.");
        }

        if (i > 0 && m_energyLevels[i] <= m_energyLevels[i - 1]) {
            throw InputFileError(routineName, node,
                "Energy levels must be strictly increasing.");
        }
    }

    if (!std::isfinite(m_threshold) || m_threshold < 0.0) {
        throw InputFileError(routineName, node,
            "The collision threshold must be finite and non-negative.");
    }
}

void ElectronCollisionPlasmaRate::setDefaultThreshold()
{
    if (m_threshold != 0.0 ||
        (m_kind != "excitation" &&
         m_kind != "ionization" &&
         m_kind != "attachment")) {
        return;
    }

    for (size_t i = 0; i < m_crossSections.size(); i++) {
        if (m_crossSections[i] > 0.0) {
            m_threshold = m_energyLevels[i];
            return;
        }
    }
}
}
