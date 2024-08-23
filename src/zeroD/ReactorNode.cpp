//! @file ReactorNode.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ReactorNode.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"

namespace Cantera
{

ReactorNode::ReactorNode(shared_ptr<Solution> sol, const string& name) : m_name(name)
{
    if (!sol || !(sol->thermo())) {
        warn_deprecated("ReactorNode::ReactorNode",
            "Creation of empty reactor objects is deprecated in Cantera 3.1 and will "
            "raise\nexceptions thereafter; reactor contents should be provided in the "
            "constructor.");
        return;
    }
    m_solution = sol;
    m_solution->thermo()->addSpeciesLock();
}

ReactorNode::~ReactorNode()
{
    if (m_solution) {
        m_solution->thermo()->removeSpeciesLock();
    }
}

void ReactorNode::setDefaultName(map<string, int>& counts)
{
    if (m_defaultNameSet) {
        return;
    }
    m_defaultNameSet = true;
    string typ(type());
    if (m_name == "(none)" || m_name == "") {
        m_name = fmt::format("{}_{}", type(), counts[type()]);
    }
    counts[type()]++;
}

void ReactorNode::restoreState()
{
    throw NotImplementedError("ReactorNode::restoreState",
        "Method needs to be overloaded.");
}

void ReactorNode::syncState()
{
    throw NotImplementedError("ReactorNode::syncState",
        "Method needs to be overloaded.");
}

}
