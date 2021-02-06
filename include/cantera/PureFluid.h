//! @file PureFluid.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CXX_PUREFLUID
#define CXX_PUREFLUID

#pragma message("warning: PureFluid.h is deprecated and will be removed after Cantera 2.5.0.")

#include "thermo/PureFluidPhase.h"
#include "kinetics.h"

namespace Cantera
{

//! Wrapper for PureFluidPhase with constructor from file
/*!
 * @deprecated To be removed after Cantera 2.5.0. Replaceable with PureFluidPhase.
 */
class PureFluid : public PureFluidPhase
{
public:
    PureFluid() : m_ok(false), m_r(0) {
        warn_deprecated("class PureFluid", "To be removed after Cantera 2.5.0. "
            "Replaceable with PureFluidPhase.");
    }

    PureFluid(const std::string& infile, std::string id="") : m_ok(false), m_r(0) {
        warn_deprecated("class PureFluid", "To be removed after Cantera 2.5.0. "
            "Replaceable with PureFluidPhase::initThermoFile.");
        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("PureFluid",
                                          "buildSolutionFromXML returned false");
    }

    PureFluid(XML_Node& root, const std::string& id) : m_ok(false), m_r(0) {
        m_ok = buildSolutionFromXML(root, id, "phase", this, 0);
    }

    bool operator!() {
        return !m_ok;
    }

    bool ready() const {
        return m_ok;
    }

    friend std::ostream& operator<<(std::ostream& s, PureFluid& mix) {
        std::string r = mix.report(true);
        s << r;
        return s;
    }

protected:
    bool m_ok;
    XML_Node* m_r;
};


//! Water definition from liquidvapor input file
/*!
 * @deprecated To be removed after Cantera 2.5.0. Replaceable with PureFluidPhase.
 */
class Water : public PureFluid
{
public:
    Water() : PureFluid(std::string("liquidvapor.cti"),std::string("water")) {
        warn_deprecated("class Water", "To be removed after Cantera 2.5.0. "
        "Replaceable with PureFluidPhase.");
    }
    virtual ~Water() {}
};

}

#endif
