//! @file BulkRate.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/BulkRate.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

ThreeBodyBase::ThreeBodyBase()
    : m_thirdBodyConc(NAN)
    , m_defaultEfficiency(1.)
    , m_specifiedCollisionPartner(false)
    , m_massAction(false)
{}

void ThreeBodyBase::setParameters(const AnyMap& node)
{
    m_defaultEfficiency = node.getDouble("default-efficiency", 1.0);
    if (node.hasKey("efficiencies")) {
        m_efficiencyMap = node["efficiencies"].asMap<double>();
    }
    m_specifiedCollisionPartner = node.getBool("specified-collider", false);
}

void ThreeBodyBase::getParameters(AnyMap& node) const
{
    if (!m_specifiedCollisionPartner) {
        if (m_efficiencyMap.size()) {
            node["efficiencies"] = m_efficiencyMap;
            node["efficiencies"].setFlowStyle();
        }
        if (m_defaultEfficiency != 1.0) {
            node["default-efficiency"] = m_defaultEfficiency;
        }
    }
}

void ThreeBodyBase::getEfficiencies(AnyMap& efficiencies) const {
    efficiencies.clear();
    for (const auto& eff : m_efficiencyMap) {
        efficiencies[eff.first] = eff.second;
    }
}

double ThreeBodyBase::efficiency(const std::string& k) const {
    return getValue(m_efficiencyMap, k, m_defaultEfficiency);
}

void ThreeBodyBase::setContext(const Reaction& rxn, const Kinetics& kin)
{
    for (const auto& eff : m_efficiencyMap) {
        size_t k = kin.kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            m_efficiencies.emplace_back(k, eff.second - m_defaultEfficiency);
        } else if (!kin.skipUndeclaredThirdBodies()) {
            throw CanteraError("ThreeBodyBase::setContext",
                "Found third-body efficiency for undeclared species '{}' "
                "while adding reaction '{}'", eff.first, rxn.equation());
        }
    }
}

}
