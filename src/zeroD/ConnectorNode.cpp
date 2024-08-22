//! @file ConnectorNode.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/ConnectorNode.h"
#include "cantera/zeroD/ReactorBase.h"

namespace Cantera
{

void ConnectorNode::setDefaultName(map<string, int>& counts)
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

}
