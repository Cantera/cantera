//! @file IncompressibleSolid.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CXX_INCOMPRESSIBLE
#define CXX_INCOMPRESSIBLE

#include "thermo/ConstDensityThermo.h"
#include "kinetics/importKinetics.h"

namespace Cantera
{

//! Wrapper for ConstDensityThermo with constructor from file
class IncompressibleSolid : public ConstDensityThermo
{
public:
    IncompressibleSolid(const std::string& infile,
                        std::string id="") : m_ok(false), m_r(0)
    {
        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("IncompressibleSolid",
                                          "buildSolutionFromXML returned false");
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
