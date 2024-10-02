//! @file Connector.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/Connector.h"
#include "cantera/zeroD/ReactorNode.h"

namespace Cantera
{

void Connector::setDefaultName(map<string, int>& counts)
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
