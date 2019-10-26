//! @file Edge.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CXX_EDGE
#define CXX_EDGE

#pragma message("warning: Edge.h is deprecated and will be removed after Cantera 2.5.0.")

#include "thermo/ThermoFactory.h"
#include "kinetics/importKinetics.h"
#include "kinetics/EdgeKinetics.h"
#include "thermo/EdgePhase.h"

namespace Cantera
{

//! Convenience class which inherits from both EdgePhase and EdgeKinetics
/*!
 * @deprecated To be removed after Cantera 2.5.0.
 *             Replaceable with Solution and/or EdgePhase/EdgeKinetics.
 */
class Edge :
    public EdgePhase, public EdgeKinetics
{
public:
    Edge(const std::string& infile, std::string id, std::vector<ThermoPhase*> phases)
        : m_ok(false), m_r(0)
    {
        warn_deprecated("class Metal",
            "To be removed after Cantera 2.5.0. "
            "Replaceable with Solution and/or EdgePhase/EdgeKinetics.");
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
