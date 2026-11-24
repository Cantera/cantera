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


namespace {

// Build a rate-integration grid containing every breakpoint of both
// the EEDF and the tabulated cross section.
//
// The cross section is evaluated as sigma(energy + sigmaArgumentShift).
// For super-elastic collisions, sigmaArgumentShift is the threshold U.
void buildMergedRateIntegrationGrid(
    const vector<double>& eedfEnergyLevels,
    const vector<double>& crossSectionEnergyLevels,
    const vector<double>& crossSections,
    double sigmaArgumentShift,
    vector<double>& integrationEnergy,
    vector<double>& integrationCrossSection)
{
    if (eedfEnergyLevels.size() < 2 ||
        !(eedfEnergyLevels.back() > eedfEnergyLevels.front())) {
        throw CanteraError("buildMergedRateIntegrationGrid",
            "Invalid EEDF energy grid.");
    }

    if (crossSectionEnergyLevels.size() < 2 ||
        crossSectionEnergyLevels.size() != crossSections.size()) {
        throw CanteraError("buildMergedRateIntegrationGrid",
            "Invalid cross-section tabulation.");
    }

    integrationEnergy.clear();
    integrationCrossSection.clear();

    // Integrate over the complete EEDF domain.
    // Cross sections use the extended table when enabled, and linearInterp()
    // applies constant extrapolation beyond the current table endpoints.
    const double lower = eedfEnergyLevels.front();
    const double upper = eedfEnergyLevels.back();

    integrationEnergy.reserve(eedfEnergyLevels.size() + crossSectionEnergyLevels.size() + 2);

    integrationEnergy.push_back(lower);
    integrationEnergy.push_back(upper);

    // Add every EEDF breakpoint inside the integration domain.
    for (double eps : eedfEnergyLevels) {
        if (eps > lower && eps < upper) {
            integrationEnergy.push_back(eps);
        }
    }

    // Add every cross-section breakpoint, shifted into the
    // electron-energy coordinates used by the EEDF.
    for (double epsSigma : crossSectionEnergyLevels) {
        const double eps = epsSigma - sigmaArgumentShift;

        if (eps > lower && eps < upper) {
            integrationEnergy.push_back(eps);
        }
    }

    std::sort(integrationEnergy.begin(), integrationEnergy.end());

    auto last = std::unique(
        integrationEnergy.begin(), integrationEnergy.end());

    integrationEnergy.erase(last, integrationEnergy.end());

    integrationCrossSection.resize(integrationEnergy.size());

    // linearInterp applies constant extrapolation outside the table.
    for (size_t i = 0; i < integrationEnergy.size(); i++) {
        integrationCrossSection[i] = linearInterp(
            integrationEnergy[i] + sigmaArgumentShift,
            crossSectionEnergyLevels,
            crossSections);
    }
}

// Integrate (eps + energyFactorShift) * F(eps) * sigma(eps)
// over a grid containing all EEDF and cross-section breakpoints.
//
// Within each interval, the EEDF and cross section are linear.
// Their product with the energy factor is a cubic polynomial.
// Two-point Gauss-Legendre quadrature is therefore exact for
// these piecewise-linear representations.
double integrateRateOnMergedGrid(
    const vector<double>& integrationEnergy,
    const vector<double>& integrationCrossSection,
    const vector<double>& eedfEnergyLevels,
    const vector<double>& eedfDistribution,
    double energyFactorShift)
{
    if (integrationEnergy.size() < 2) {
        return 0.0;
    }

    if (integrationEnergy.size() != integrationCrossSection.size()) {
        throw CanteraError("integrateRateOnMergedGrid",
            "Integration energy and cross-section arrays have "
            "different sizes.");
    }

    if (eedfEnergyLevels.size() != eedfDistribution.size() ||
        eedfEnergyLevels.size() < 2) {
        throw CanteraError("integrateRateOnMergedGrid",
            "Invalid EEDF arrays.");
    }

    // Two Gauss-Legendre nodes mapped from [-1, 1] to [0, 1].
    constexpr double invSqrt3 = 0.577350269189625764509148780501957456;
    constexpr double t1 = 0.5 * (1.0 - invSqrt3);
    constexpr double t2 = 0.5 * (1.0 + invSqrt3);

    double integral = 0.0;

    double fLeft = linearInterp(integrationEnergy.front(),
        eedfEnergyLevels, eedfDistribution);

    for (size_t i = 0; i + 1 < integrationEnergy.size(); i++) {
        const double a = integrationEnergy[i];
        const double b = integrationEnergy[i + 1];
        const double h = b - a;

        if (!(h > 0.0)) {
            throw CanteraError("integrateRateOnMergedGrid",
                "Integration grid is not strictly increasing.");
        }

        const double sigmaLeft = integrationCrossSection[i];
        const double sigmaRight = integrationCrossSection[i + 1];

        const double fRight = linearInterp(b, eedfEnergyLevels, eedfDistribution);

        // Both functions are linear throughout [a, b].
        const double f1 = fLeft + t1 * (fRight - fLeft);
        const double f2 = fLeft + t2 * (fRight - fLeft);

        const double sigma1 = sigmaLeft + t1 * (sigmaRight - sigmaLeft);
        const double sigma2 = sigmaLeft + t2 * (sigmaRight - sigmaLeft);

        const double eps1 = a + t1 * h;
        const double eps2 = a + t2 * h;

        const double integrand1 = (eps1 + energyFactorShift) * f1 * sigma1;
        const double integrand2 = (eps2 + energyFactorShift) * f2 * sigma2;

        integral += 0.5 * h * (integrand1 + integrand2);

        fLeft = fRight;
    }

    if (!std::isfinite(integral) || integral < 0.0) {
        throw CanteraError("integrateRateOnMergedGrid",
            "Invalid rate integral: {}.", integral);
    }

    return integral;
}

} // namespace

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

void ElectronCollisionPlasmaRate::invalidateRateIntegrationCaches()
{
    // Forward-rate integration cache.
    m_rateIntegrationEnergy.clear();
    m_rateIntegrationCrossSection.clear();
    m_rateIntegrationLevelNumber = -4;

    // Super-elastic integration cache.
    m_superElasticRateIntegrationEnergy.clear();
    m_superElasticRateIntegrationCrossSection.clear();
    m_superElasticRateIntegrationLevelNumber = -5;
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

    // Store the local preference without applying extrapolation yet.
    m_extrapolateCrossSections = node.getBool("extrapolate-cross-sections", false);

    m_crossSectionsExtrapolated = false;

    validateCollisionData(node);
    setDefaultThreshold();
    enforceThresholdCrossSections();

    // Keep track of the table before optional high-energy extrapolation.
    m_tabulatedPointCount = m_energyLevels.size();

    // Invalidate interpolation caches after loading new data.
    m_levelNumber = -3;
    m_levelNumberSuperelastic = -2;
    m_crossSectionsInterpolated.clear();
    m_crossSectionsOffset.resize(0);
    invalidateRateIntegrationCaches();

    m_hasCrossSectionData = true;
}

void ElectronCollisionPlasmaRate::getParameters(
    AnyMap& node) const
{
    node["type"] = type();

    if (!m_hasCrossSectionData) {
        return;
    }

    // Serialize the original (possibly threshold-normalized) table,
    // without the high-energy points generated at runtime if any extrapolation is present.
    node["energy-levels"] = vector<double>(m_energyLevels.begin(), m_energyLevels.begin() + m_tabulatedPointCount);

    node["cross-sections"] = vector<double>(m_crossSections.begin(), m_crossSections.begin() + m_tabulatedPointCount);

    // Preserve the opt-in setting across YAML round-trips.
    if (m_extrapolateCrossSections) {
        node["extrapolate-cross-sections"] = true;
    }

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

    // Rebuild the integration grid only when the EEDF grid changes.
    // Changes in the EEDF values do not affect this geometric cache.
    if (m_rateIntegrationLevelNumber != shared_data.levelNumber) {
        buildMergedRateIntegrationGrid(shared_data.energyLevels,
            m_energyLevels, m_crossSections, 0.0,
            m_rateIntegrationEnergy, m_rateIntegrationCrossSection);

        m_rateIntegrationLevelNumber = shared_data.levelNumber;
    }

    // Integrate eps * F0(eps) * sigma(eps) using all breakpoints.
    const double integral = integrateRateOnMergedGrid(
        m_rateIntegrationEnergy, m_rateIntegrationCrossSection,
        shared_data.energyLevels, shared_data.distribution, 0.0);

    // Unit: m3/kmol/s
    return pow(2.0 * ElectronCharge / ElectronMass, 0.5)
        * Avogadro * integral;

}

void ElectronCollisionPlasmaRate::modifyRateConstants(
    const ElectronCollisionPlasmaData& shared_data, double& kf, double& kr)
{
    if (kr == 0.0) {
        // The reverse rate constant is only for reversible reactions
        // kr = 0.0 indicates that the reaction is irreversible
        return;
    }

    // Use the physical excitation threshold for the energy shift.
    const double U = m_threshold;

    // Rebuild the super-elastic integration grid only when the
    // EEDF energy grid changes.
    //
    // The reverse integral samples sigma_exc(eps + U).
    // Native cross-section energies are therefore shifted by -U.
    if (m_superElasticRateIntegrationLevelNumber != shared_data.levelNumber) {
        buildMergedRateIntegrationGrid(shared_data.energyLevels,
            m_energyLevels, m_crossSections,
            U, m_superElasticRateIntegrationEnergy,
            m_superElasticRateIntegrationCrossSection);

        m_superElasticRateIntegrationLevelNumber =
            shared_data.levelNumber;
    }

    // Detailed balance gives:
    //
    // eps * sigma_super(eps)
    //     = (eps + U) * sigma_exc(eps + U).
    //
    // Integrate (eps + U) * F0(eps) * sigma_exc(eps + U).
    const double integral = integrateRateOnMergedGrid(
        m_superElasticRateIntegrationEnergy,
        m_superElasticRateIntegrationCrossSection,
        shared_data.energyLevels,
        shared_data.distribution,
        U);

    // Unit: m3/kmol/s
    kr = pow(2.0 * ElectronCharge / ElectronMass, 0.5)
        * Avogadro * integral;
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
    enforceThresholdCrossSections();
    // Update the original table size before phase-level extrapolation.
    if (!m_crossSectionsExtrapolated) {
        m_tabulatedPointCount = m_energyLevels.size();
    }

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
            // Use the last zero cross-section point before the first positive one.
            m_threshold = (i > 0 && m_crossSections[i - 1] == 0.0)
                ? m_energyLevels[i - 1] : m_energyLevels[i];
            return;
        }
    }
}

void ElectronCollisionPlasmaRate::enforceThresholdCrossSections()
{
    // Only excitation and ionization have an onset threshold in the current collision model.
    if (m_kind != "excitation" && m_kind != "ionization") {
        return;
    }

    // Remove all data points below the threshold and insert a zero cross section at the threshold itself.
    vector<double> levels{m_threshold};
    vector<double> sections{0.0};

    for (size_t i = 0; i < m_energyLevels.size(); i++) {
        if (m_energyLevels[i] > m_threshold) {
            levels.push_back(m_energyLevels[i]);
            sections.push_back(m_crossSections[i]);
        }
    }

    // A valid normalized table must contain at least one energy point strictly above the threshold.
    if (levels.size() < 2) {
        throw CanteraError("ElectronCollisionPlasmaRate::enforceThresholdCrossSections",
            "Collision '{}' requires a tabulated energy above its threshold ({} eV).",
            m_kind, m_threshold);
    }

    // Avoid modifying an already normalized table. This makes the operation idempotent.
    if (m_energyLevels == levels && m_crossSections == sections) {
        return;
    }

    m_energyLevels = std::move(levels);
    m_crossSections = std::move(sections);
    invalidateRateIntegrationCaches();

    // Invalidate interpolation caches after modifying the table.
    m_levelNumber = -3;
    m_levelNumberSuperelastic = -2;
    m_crossSectionsInterpolated.clear();
    m_crossSectionsOffset.resize(0);
}

void ElectronCollisionPlasmaRate::extrapolateCrossSections()
{
    if (m_crossSectionsExtrapolated) {
        return;
    }

    // Remember the original table size for YAML serialization.
    m_tabulatedPointCount = m_energyLevels.size();
    m_crossSectionsExtrapolated = true;

    const double lastEnergy = m_energyLevels.back();
    const double lastCrossSection = m_crossSections.back();

    // Extend the table up to an energy above 10 keV. A zero endpoint already gives a zero high-energy tail.
    if (lastEnergy >= 1.0e4 || lastCrossSection == 0.0) {
        return;
    }

    // Add logarithmically spaced energy points above the original table, using a factor of 1.5 between points.
    double energy = lastEnergy;

    do {
        energy *= 1.5;

        // Continuous high-energy tail with asymptotic behavior proportional to log(E) / E.
        const double crossSection = lastCrossSection * (lastEnergy / energy)
            * (1.0 + std::log(energy / lastEnergy));

        m_energyLevels.push_back(energy);
        m_crossSections.push_back(crossSection);

    } while (energy <= 1.0e4);

    // Invalidate cached interpolations after extending the table.
    m_levelNumber = -3;
    m_levelNumberSuperelastic = -2;
    m_crossSectionsInterpolated.clear();
    m_crossSectionsOffset.resize(0);

    invalidateRateIntegrationCaches();
}

void ElectronCollisionPlasmaRate::setCrossSectionExtrapolation(bool enabled)
{
    if (!m_hasCrossSectionData) {
        return;
    }

    if (enabled) {
        extrapolateCrossSections();
        return;
    }

    if (!m_crossSectionsExtrapolated) {
        return;
    }

    // Restore the table before high-energy extrapolation.
    // Keep any threshold normalization already performed.
    m_energyLevels.resize(m_tabulatedPointCount);
    m_crossSections.resize(m_tabulatedPointCount);

    m_crossSectionsExtrapolated = false;

    // Invalidate interpolated data after changing the table.
    m_levelNumber = -3;
    m_levelNumberSuperelastic = -2;
    m_crossSectionsInterpolated.clear();
    m_crossSectionsOffset.resize(0);

    invalidateRateIntegrationCaches();
}
}
