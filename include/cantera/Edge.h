//! @file Edge.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CXX_EDGE
#define CXX_EDGE

#include "thermo/ThermoFactory.h"
#include "kinetics/importKinetics.h"
#include "kinetics/EdgeKinetics.h"
#include "thermo/EdgePhase.h"

namespace Cantera
{

//! Convenience class which inherits from both EdgePhase and EdgeKinetics
class Edge :
    public EdgePhase, public EdgeKinetics
{
public:
    Edge(const std::string& infile, std::string id, std::vector<ThermoPhase*> phases)
        : m_ok(false), m_r(0)
    {
        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }

        XML_Node* x = get_XML_Node("#"+id, m_r);
        if (!x) {
            throw CanteraError("Edge::Edge","error in get_XML_Node");
        }

        importPhase(*x, this);
        phases.push_back(this);
        importKinetics(*x, phases, this);
        m_ok = true;
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }

protected:
    bool m_ok;
    XML_Node* m_r;
};
}

#endif
