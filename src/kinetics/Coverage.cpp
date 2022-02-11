//! @file Coverage.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/Coverage.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
Coverage::Coverage()
    : m_siteDensity(NAN)
    , m_acov(0.)
    , m_ecov(0.)
    , m_mcov(0.)
{
}

void Coverage::setCoverageDependencies(const AnyMap& dependencies,
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

void Coverage::getCoverageDependencies(AnyMap& dependencies, bool asVector) const
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

void Coverage::addCoverageDependence(std::string sp, double a, double m, double e)
{
    if (std::find(m_cov.begin(), m_cov.end(), sp) == m_cov.end()) {
        m_cov.push_back(sp);
        m_ac.push_back(a);
        m_ec.push_back(e);
        m_mc.push_back(m);
        m_indices.clear();
    } else {
        throw CanteraError("Coverage::addCoverageDependence",
            "Coverage for species '{}' is already specified.", sp);
    }
}

void Coverage::setSpecies(const std::vector<std::string>& species)
{
    m_indices.clear();
    for (size_t k = 0; k < m_cov.size(); k++) {
        auto it = find(species.begin(), species.end(), m_cov[k]);
        if (it != species.end()) {
            m_indices.emplace(k, it - species.begin());
        } else {
            throw CanteraError("Coverage:setSpeciesIndices",
                "Species list does not contain '{}'.", m_cov[k]);
        }
    }
}

void Coverage::setSpecies(const Kinetics& kin) {
    setSpecies(kin.thermo().speciesNames());
}

StickCoverage::StickCoverage()
    : m_motz_wise(false)
    , m_stickingSpecies("")
    , m_explicitSpecies(false)
    , m_surfaceOrder(NAN)
    , m_multiplier(NAN)
    , m_factor(NAN)
{
}

void StickCoverage::setStickParameters(const AnyMap& node)
{
    m_motz_wise = node.getBool("Motz-Wise", false);
    m_explicitSpecies = node.hasKey("sticking-species");
    m_stickingSpecies = node.getString("sticking-species", "");
}

void StickCoverage::getStickParameters(AnyMap& node) const
{
    if (m_motz_wise) {
        node["Motz-Wise"] = true;
    }
    if (m_explicitSpecies) {
        node["sticking-species"] = m_stickingSpecies;
    }
}

void StickCoverage::buildStickCoefficients(const Reaction& rxn, const Kinetics& kin)
{
    // Identify the interface phase
    size_t iInterface = npos;
    size_t min_dim = 4;
    for (size_t n = 0; n < kin.nPhases(); n++) {
        if (kin.thermo(n).nDim() < min_dim) {
            iInterface = n;
            min_dim = kin.thermo(n).nDim();
        }
    }

    std::string sticking_species = m_stickingSpecies;
    if (sticking_species == "") {
        // Identify the sticking species if not explicitly given
        bool foundStick = false;
        for (const auto& sp : rxn.reactants) {
            size_t iPhase = kin.speciesPhaseIndex(kin.kineticsSpeciesIndex(sp.first));
            if (iPhase != iInterface) {
                // Non-interface species. There should be exactly one of these
                if (foundStick) {
                    throw InputFileError("StickCoverage::buildStickCoefficients",
                        rxn.input, "Multiple non-interface species ('{}' and '{}')\n"
                        "found in sticking reaction: '{}'.\nSticking species "
                        "must be explicitly specified.",
                        sticking_species, sp.first, rxn.equation());
                }
                foundStick = true;
                sticking_species = sp.first;
            }
        }
        if (!foundStick) {
            throw InputFileError("StickCoverage::buildStickCoefficients",
                rxn.input, "No non-interface species found "
                "in sticking reaction: '{}'", rxn.equation());
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
            const auto& surf = dynamic_cast<const SurfPhase&>(kin.thermo());
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
