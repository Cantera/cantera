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

void ElectronCollisionPlasmaRate::setParameters(
    const AnyMap& node, const UnitStack& rate_units)
{
    ReactionRate::setParameters(node, rate_units);

    // Cross-section data are stored directly in the rate object.
    const bool hasLevels = node.hasKey("energy-levels");
    const bool hasSections = node.hasKey("cross-sections");

    // Preserve support for default-constructed rate objects.
    if (!hasLevels && !hasSections) {
        return;
    }

    if (hasLevels != hasSections) {
        throw InputFileError(
            "ElectronCollisionPlasmaRate::setParameters", node,
            "'energy-levels' and 'cross-sections' must be provided together.");
    }

    // The target of a reactive collision comes from its equation.
    if (node.hasKey("equation") && node.hasKey("target")) {
        throw InputFileError(
            "ElectronCollisionPlasmaRate::setParameters", node,
            "A reactive electron collision must infer its target from the equation.");
    }

    // Kind and target can be absent until the reaction context has been analyzed.
    m_kind = node.getString("kind", "");
    m_target = node.getString("target", "");
    m_product = node.getString("product", "");

    m_energyLevels = node["energy-levels"].asVector<double>();
    m_crossSections = node["cross-sections"].asVector<double>(m_energyLevels.size());

    m_hasExplicitThreshold = node.hasKey("threshold");
    m_threshold = node.getDouble("threshold", 0.0);
    m_isReactionRate = false;

    validateCollisionData(node);
    setDefaultThreshold();

    // Invalidate interpolation caches after loading new data.
    m_levelNumber = -3;
    m_levelNumberSuperelastic = -2;
    m_crossSectionsInterpolated.clear();
    m_crossSectionsOffset.resize(0);

    m_hasCrossSectionData = true;
}

void ElectronCollisionPlasmaRate::getParameters(
    AnyMap& node) const
{
    node["type"] = type();

    if (!m_hasCrossSectionData) {
        return;
    }

    node["energy-levels"] = m_energyLevels;
    node["cross-sections"] = m_crossSections;

    if (!m_kind.empty()) {
        node["kind"] = m_kind;
    }

    // Reactive entries infer the target from their equation.
    if (!m_isReactionRate && !m_target.empty()) {
        node["target"] = m_target;
    }

    if (!m_product.empty()) {
        node["product"] = m_product;
    }

    if (m_hasExplicitThreshold || m_threshold != 0.0) {
        node["threshold"] = m_threshold;
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

void ElectronCollisionPlasmaRate::setContext(
    const Reaction& rxn, const Kinetics& kin)
{
    const ThermoPhase& thermo = kin.thermo();

    if (thermo.type() != "plasma") {
        throw CanteraError(
            "ElectronCollisionPlasmaRate::setContext",
            "ElectronCollisionPlasmaRate requires a plasma phase.");
    }

    const auto& plasma = dynamic_cast<const PlasmaPhase&>(thermo);

    // All metadata inference is performed independently of Kinetics.
    setReactionContext(rxn, plasma.electronSpeciesName());
}

void ElectronCollisionPlasmaRate::setReactionContext(
    const Reaction& rxn, const string& electronName)
{
    const string routine =
        "ElectronCollisionPlasmaRate::setReactionContext";

    if (!m_hasCrossSectionData) {
        throw InputFileError(routine, rxn.input,
            "Electron-collision reaction '{}' has no tabulated cross sections.", rxn.equation());
    }

    if (rxn.reactants.size() != 2) {
        throw InputFileError(routine, rxn.input,
            "Electron-collision reactions require exactly two reactant species.");
    }

    auto electron = rxn.reactants.find(electronName);

    if (electron == rxn.reactants.end() || electron->second != 1.0) {
        throw InputFileError(routine, rxn.input,
            "Electron-collision reactions require exactly one electron reactant.");
    }

    // Infer the heavy-particle target from reactants.
    string target;
    for (const auto& [name, coefficient] : rxn.reactants) {
        if (name != electronName) {
            target = name;
        }
    }

    if (!m_target.empty() && m_target != target) {
        throw InputFileError(routine, rxn.input,
            "Collision target '{}' disagrees with '{}'.", m_target, target);
    }

    m_target = target;

    // Determine the number of electrons in the products.
    auto productElectron = rxn.products.find(electronName);
    const double productElectrons = (productElectron == rxn.products.end()) ? 0.0 : productElectron->second;

    string inferredKind;

    if (productElectrons > 1.0) {
        inferredKind = "ionization";
    } else if (productElectrons < 1.0) {
        inferredKind = "attachment";
    } else if (rxn.reactants != rxn.products) {
        inferredKind = "excitation";
    }

    // Identical reactants and products do not uniquely
    // determine the collision kind.
    if (inferredKind.empty() && m_kind.empty()) {
        throw InputFileError(routine, rxn.input,
            "The collision 'kind' is required when reactants and products are identical.");
    }

    // An unchanged electron-target reaction cannot
    // represent ionization or attachment.
    if (inferredKind.empty() && !m_kind.empty()
        && m_kind != "effective"
        && m_kind != "elastic"
        && m_kind != "excitation") {
        throw InputFileError(routine, rxn.input, "Collision kind '{}' is incompatible with an unchanged electron-target reaction.", m_kind);
    }

    if (!inferredKind.empty()) {
        if (!m_kind.empty() && m_kind != inferredKind) {
            throw InputFileError(routine, rxn.input,
                "Collision kind '{}' disagrees with the inferred kind '{}'.",
                m_kind, inferredKind);
        }

        m_kind = inferredKind;
    }

    // The method can be called first by PlasmaPhase and
    // again by Kinetics without changing the result.
    m_isReactionRate = true;
    setDefaultThreshold();

    // Validate superelastic reverse-reaction requirements.
    if (rxn.reversible && (rxn.products.size() != 2 ||
        productElectron == rxn.products.end() || productElectron->second != 1.0)) {
        throw InputFileError(routine, rxn.input,
            "A reversible electron-collision reaction requires exactly one electron and one other product species.");
    }
}

void ElectronCollisionPlasmaRate::validateCollisionData(
    const AnyMap& node) const
{
    const string routineName = "ElectronCollisionPlasmaRate::validateCollisionData";

    static const set<string> validKinds = {
        "effective", "elastic", "excitation", "ionization", "attachment"
    };

    // An empty kind is allowed before the reaction context is known.
    if (!m_kind.empty() && !validKinds.count(m_kind)) {
        throw InputFileError(routineName, node,
            "Unknown electron-collision kind '{}'. Expected one of "
            "'effective', 'elastic', 'excitation', 'ionization', or "
            "'attachment'.",
            m_kind);
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
