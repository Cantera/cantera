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

}
