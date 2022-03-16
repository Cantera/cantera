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
        m_efficiencies = node["efficiencies"].asMap<double>();
    }
    m_specifiedCollisionPartner = node.getBool("specified-collider", false);
    m_massAction = node.getBool("mass-action", true);
}

void ThreeBodyBase::getParameters(AnyMap& node) const
{
    if (!m_specifiedCollisionPartner) {
        if (m_efficiencies.size()) {
            node["efficiencies"] = m_efficiencies;
            node["efficiencies"].setFlowStyle();
        }
        if (m_defaultEfficiency != 1.0) {
            node["default-efficiency"] = m_defaultEfficiency;
        }
    }
}

void ThreeBodyBase::setEfficiencies(const Composition& efficiencies)
{
    if (m_efficiencyMap.size()) {
        throw CanteraError("ThreeBodyBase::setEfficiencies",
            "Unable to set efficiencies once they have been processed.");
    }
    m_efficiencies = efficiencies;
}

void ThreeBodyBase::getEfficiencyMap(std::map<size_t, double>& eff) const
{
    eff.clear();
    for (const auto& item : m_efficiencyMap) {
        eff[item.first] = item.second + m_defaultEfficiency;
    }
}

double ThreeBodyBase::efficiency(const std::string& k) const {
    return getValue(m_efficiencies, k, m_defaultEfficiency);
}

void ThreeBodyBase::setContext(const Reaction& rxn, const Kinetics& kin)
{
    for (const auto& eff : m_efficiencies) {
        size_t k = kin.kineticsSpeciesIndex(eff.first);
        if (k != npos) {
            m_efficiencyMap.emplace_back(k, eff.second - m_defaultEfficiency);
        } else if (!kin.skipUndeclaredThirdBodies()) {
            throw CanteraError("ThreeBodyBase::setContext",
                "Found third-body efficiency for undeclared species '{}' "
                "while adding reaction '{}'", eff.first, rxn.equation());
        }
    }
}

}
