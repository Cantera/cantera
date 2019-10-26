//! @file Metal.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CXX_METAL
#define CXX_METAL

#pragma message("warning: Metal.h is deprecated and will be removed after Cantera 2.5.0.")

#include "thermo/MetalPhase.h"
#include "kinetics/importKinetics.h"

namespace Cantera
{

//! Wrapper for MetalPhase with constructor from file
/*!
 * @deprecated To be removed after Cantera 2.5.0.
 *             Replaceable with Solution and/or MetalPhase.
 */
class Metal : public MetalPhase
{
public:
    Metal(const std::string& infile, std::string id="") : m_ok(false), m_r(0)
    {
        warn_deprecated("class Metal",
            "To be removed after Cantera 2.5.0. "
            "Replaceable with Solution and/or MetalPhase.");
        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("Metal::Metal",
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
