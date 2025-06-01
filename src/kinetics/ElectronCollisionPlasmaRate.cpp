//! @file ElectronCollisionPlasmaRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Kinetics.h"
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
    pp.getElectronEnergyDistribution(distribution.data());

    // Update energy levels
    // levelChanged = pp.levelNumber() != m_level_number;
    // if (levelChanged) {
        m_level_number = pp.levelNumber();
        energyLevels.resize(pp.nElectronEnergyLevels());
        pp.getElectronEnergyLevels(energyLevels.data());
    // }

    return true;
}

void ElectronCollisionPlasmaRate::setParameters(const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);

    //  **Extract kind, target, and product from reaction node**
    if (node.hasKey("kind")) {
        m_kind = node["kind"].asString();
    } /*else {
        throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
                           "Missing `kind` field in electron-collision-plasma reaction.");
    }*/

    if (node.hasKey("target")) {
        m_target = node["target"].asString();
    } /*else {
        throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
                           "Missing `target` field in electron-collision-plasma reaction.");
    }*/

    if (node.hasKey("product")) {
        m_product = node["product"].asString();
    } /*else {
        throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
                           "Missing `product` field in electron-collision-plasma reaction.");
    }*/

    //  **First, check if cross-sections are embedded in the reaction itself**
    if (node.hasKey("energy-levels") && node.hasKey("cross-sections")) {
        //writelog("Using embedded cross-section data from reaction definition.\n");

        m_energyLevels = node["energy-levels"].asVector<double>();
        m_crossSections = node["cross-sections"].asVector<double>();

        if (m_energyLevels.size() != m_crossSections.size()) {
            throw CanteraError("ElectronCollisionPlasmaRate::setParameters",
                               "Mismatch: `energy-levels` and `cross-sections` must have the same length.");
        }

        cs_ok = true;  // Mark as valid cross-section data
    }

    //  **If no cross-section data was found, defer to PlasmaPhase (old format)**
    else {
        //writelog("No cross-section data found in reaction, relying on PlasmaPhase initialization.\n");
        cs_ok = false;  //  This will be handled later in `PlasmaPhase::initThermo()`
    }
}

void ElectronCollisionPlasmaRate::getParameters(AnyMap& node) const {
    node["type"] = type();
    node["energy-levels"] = m_energyLevels;
    node["cross-sections"] = m_crossSections;
}

void ElectronCollisionPlasmaRate::updateInterpolatedCrossSection(
    const vector<double>& sharedLevels) {
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
    // the electron energy distribution function when the interpolated
    // cross section is empty
    // Note that the PlasmaPhase should handle the interpolated cross sections
    // for all ElectronCollisionPlasmaRate reactions
    if (m_crossSectionsInterpolated.size() == 0) {
        for (double level : shared_data.energyLevels) {
            m_crossSectionsInterpolated.push_back(linearInterp(level,
                                                  m_energyLevels, m_crossSections));
        }
    }

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

    if (m_crossSectionsOffset.size() != shared_data.energyLevels.size()) {
        m_crossSectionsOffset.resize(shared_data.energyLevels.size());
        vector<double> superElasticEnergyLevels{0.0};
        for (size_t i = 1; i < m_energyLevels.size(); i++) {
            superElasticEnergyLevels.push_back(m_energyLevels[i] - m_energyLevels[0]);
        }
        for (size_t i = 0; i < shared_data.energyLevels.size(); i++) {
            m_crossSectionsOffset[i] = linearInterp(shared_data.energyLevels[i],
                                                    superElasticEnergyLevels,
                                                    m_crossSections);
        }
    }

    // Interpolate cross-sections data to the energy levels of
    // the electron energy distribution function
    if (shared_data.levelChanged) {
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
    // get electron species name
    string electronName;
    if (kin.thermo().type() == "plasma") {
        electronName = dynamic_cast<const PlasmaPhase&>(kin.thermo()).electronSpeciesName();
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

}
