//! @file Coverage.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/InterfaceRate.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
CoverageBase::CoverageBase()
    : m_siteDensity(NAN)
    , m_acov(0.)
    , m_ecov(0.)
    , m_mcov(0.)
    , m_chargeTransfer(false)
    , m_beta(0.5)
    , m_exchangeCurrentDensityFormulation(false)
    , m_deltaG0(NAN)
    , m_prodStandardConcentrations(NAN)
{
}

void CoverageBase::setParameters(const AnyMap& node)
{
    if (node.hasKey("coverage-dependencies")) {
        setCoverageDependencies(
            node["coverage-dependencies"].as<AnyMap>(), node.units());
    }
    if (node.hasKey("beta")) {
        m_beta = node["beta"].asDouble();
    }
    m_exchangeCurrentDensityFormulation = node.getBool(
        "exchange-current-density-formulation", false);
}

void CoverageBase::getParameters(AnyMap& node) const
{
    if (!m_cov.empty()) {
        AnyMap deps;
        getCoverageDependencies(deps);
        node["coverage-dependencies"] = std::move(deps);
    }
    if (m_chargeTransfer) {
        if (m_beta != 0.5) {
            node["beta"] = m_beta;
        }
        if (m_exchangeCurrentDensityFormulation) {
            node["exchange-current-density-formulation"] = true;
        }
    }
}

void CoverageBase::setCoverageDependencies(const AnyMap& dependencies,
                                           const UnitSystem& units)
{
    m_cov.clear();
    m_ac.clear();
    m_ec.clear();
    m_mc.clear();
    for (const auto& item : dependencies) {
        double a, E, m;
        if (item.second.is<AnyMap>()) {
            auto& cov_map = item.second.as<AnyMap>();
            a = cov_map["a"].asDouble();
            m = cov_map["m"].asDouble();
            E = units.convertActivationEnergy(cov_map["E"], "K");
        } else {
            auto& cov_vec = item.second.asVector<AnyValue>(3);
            a = cov_vec[0].asDouble();
            m = cov_vec[1].asDouble();
            E = units.convertActivationEnergy(cov_vec[2], "K");
        }
        addCoverageDependence(item.first, a, m, E);
    }
}

void CoverageBase::getCoverageDependencies(AnyMap& dependencies, bool asVector) const
{
    for (size_t k = 0; k < m_cov.size(); k++) {
        if (asVector) {
            // this preserves the previous 'coverage_deps' units
            vector_fp dep(3);
            dep[0] = m_ac[k];
            dep[1] = m_mc[k];
            dep[2] = m_ec[k];
            dependencies[m_cov[k]] = std::move(dep);
        } else {
            AnyMap dep;
            dep["a"] = m_ac[k];
            dep["m"] = m_mc[k];
            dep["E"].setQuantity(m_ec[k], "K", true);
            dependencies[m_cov[k]] = std::move(dep);
        }
    }
}

void CoverageBase::addCoverageDependence(const std::string& sp,
                                         double a, double m, double e)
{
    if (std::find(m_cov.begin(), m_cov.end(), sp) == m_cov.end()) {
        m_cov.push_back(sp);
        m_ac.push_back(a);
        m_ec.push_back(e);
        m_mc.push_back(m);
        m_indices.clear();
    } else {
        throw CanteraError("CoverageBase::addCoverageDependence",
            "Coverage for species '{}' is already specified.", sp);
    }
}

void CoverageBase::setSpecies(const std::vector<std::string>& species)
{
    m_indices.clear();
    for (size_t k = 0; k < m_cov.size(); k++) {
        auto it = find(species.begin(), species.end(), m_cov[k]);
        if (it != species.end()) {
            m_indices.emplace(k, it - species.begin());
        } else {
            throw CanteraError("CoverageBase:setSpeciesIndices",
                "Species list does not contain '{}'.", m_cov[k]);
        }
    }
}

void CoverageBase::setContext(const Reaction& rxn, const Kinetics& kin)
{
    m_chargeTransfer = rxn.checkElectrochemistry(kin);
    setSpecies(kin.thermo().speciesNames());

    m_stoichCoeffs.clear();
    for (const auto& sp : rxn.reactants) {
        m_stoichCoeffs.emplace_back(kin.kineticsSpeciesIndex(sp.first), -sp.second);
    }
    for (const auto& sp : rxn.products) {
        m_stoichCoeffs.emplace_back(kin.kineticsSpeciesIndex(sp.first), sp.second);
    }
}

StickingCoverage::StickingCoverage()
    : m_motzWise(false)
    , m_explicitMotzWise(false)
    , m_stickingSpecies("")
    , m_explicitSpecies(false)
    , m_surfaceOrder(NAN)
    , m_multiplier(NAN)
    , m_factor(NAN)
{
}

void StickingCoverage::setStickingParameters(const AnyMap& node)
{
    m_motzWise = node.getBool("Motz-Wise", false);
    m_explicitMotzWise = node.hasKey("Motz-Wise");
    m_stickingSpecies = node.getString("sticking-species", "");
    m_explicitSpecies = node.hasKey("sticking-species");
}

void StickingCoverage::getStickingParameters(AnyMap& node) const
{
    if (m_explicitMotzWise) {
        node["Motz-Wise"] = m_motzWise;
    }
    if (m_explicitSpecies) {
        node["sticking-species"] = m_stickingSpecies;
    }
}

void StickingCoverage::setContext(const Reaction& rxn, const Kinetics& kin)
{
    // Ensure that site density is initialized
    const ThermoPhase& phase = kin.thermo(kin.surfacePhaseIndex());
    const auto& surf = dynamic_cast<const SurfPhase&>(phase);
    m_siteDensity = surf.siteDensity();
    if (!m_explicitMotzWise) {
        m_motzWise = kin.thermo().input().getBool("Motz-Wise", false);
    }

    // Identify the interface phase
    size_t iInterface = kin.reactionPhaseIndex();

    std::string sticking_species = m_stickingSpecies;
    if (sticking_species == "") {
        // Identify the sticking species if not explicitly given
        std::vector<std::string> gasSpecies;
        std::vector<std::string> anySpecies;
        for (const auto& sp : rxn.reactants) {
            size_t iPhase = kin.speciesPhaseIndex(kin.kineticsSpeciesIndex(sp.first));
            if (iPhase != iInterface) {
                // Non-interface species. There should be exactly one of these
                // (either in gas phase or other phase)
                if (kin.thermo(iPhase).phaseOfMatter() == "gas") {
                    gasSpecies.push_back(sp.first);
                }
                anySpecies.push_back(sp.first);
            }
        }
        if (gasSpecies.size() == 1) {
            // single sticking species in gas phase
            sticking_species = gasSpecies[0];
        } else if (anySpecies.size() == 1) {
            // single sticking species in any phase
            sticking_species = anySpecies[0];
        } else if (anySpecies.size() == 0) {
            throw InputFileError("StickingCoverage::setContext",
                rxn.input, "No non-interface species found "
                "in sticking reaction: '{}'", rxn.equation());
        } else {
            throw InputFileError("StickingCoverage::setContext",
                rxn.input, "Multiple non-interface species ({})\nfound in sticking "
                "reaction: '{}'.\nSticking species must be explicitly specified.",
                fmt::format("'{}'", fmt::join(anySpecies, "', '")), rxn.equation());
        }
    }
    m_stickingSpecies = sticking_species;

    double surface_order = 0.0;
    double multiplier = 1.0;
    // Adjust the A-factor
    for (const auto& sp : rxn.reactants) {
        size_t iPhase = kin.speciesPhaseIndex(kin.kineticsSpeciesIndex(sp.first));
        const ThermoPhase& p = kin.thermo(iPhase);
        size_t k = p.speciesIndex(sp.first);
        if (sp.first == sticking_species) {
            multiplier *= sqrt(GasConstant / (2 * Pi * p.molecularWeight(k)));
        } else {
            // Non-sticking species. Convert from coverages used in the
            // sticking probability expression to the concentration units
            // used in the mass action rate expression. For surface phases,
            // the dependence on the site density is incorporated when the
            // rate constant is evaluated, since we don't assume that the
            // site density is known at this time.
            double order = getValue(rxn.orders, sp.first, sp.second);
            if (&p == &surf) {
                multiplier *= pow(surf.size(k), order);
                surface_order += order;
            } else {
                multiplier *= pow(p.standardConcentration(k), -order);
            }
        }
    }
    m_surfaceOrder = surface_order;
    m_multiplier = multiplier;
}

}
