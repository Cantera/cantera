//! @file GRI30.h
#ifndef CXX_GRI30H
#define CXX_GRI30H

#include "thermo/IdealGasPhase.h"
#include "kinetics/GRI_30_Kinetics.h"
#include "kinetics/importKinetics.h"
#include "base/stringUtils.h"

namespace Cantera
{

/**
 * This class is a convenience class for use in C++ programs that
 * hard-wires the GRI 3.0 reaction mechanism. It derivees from
 * both Cantera::IdealGasPhase, which handles all composition and
 * state information, as well as thermodynamic properties, and
 * class GRI_30_Kinetics, which is the kinetics manager with
 * hard-wired replacements for some of the generic kinetics
 * methods like "getNetReactionRates."
 * @deprecated
 */
class GRI30 :
    public IdealGasPhase,
    public GRI_30_Kinetics
{
public:
    GRI30() : m_ok(false), m_r(0) {
        m_r = get_XML_File("gri30.xml");
        m_ok = buildSolutionFromXML(*m_r, "gri30",
                                    "phase", this, this);
        if (!m_ok) throw CanteraError("GRI30",
                                          "buildSolutionFromXML returned false");
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }
    friend std::ostream& operator<<(std::ostream& s, GRI30& mix) {
        std::string r = mix.report(true);
        s << r;
        return s;
    }

protected:
    bool m_ok;
    Cantera::XML_Node* m_r;

private:
};
}


#endif
