/*
 * $Id: GRI30.h,v 1.8 2008/02/16 21:35:30 hkmoffa Exp $
 */
#ifndef CXX_GRI30H
#define CXX_GRI30H

#include <string>

#include "kernel/IdealGasPhase.h"
#include "kernel/GRI_30_Kinetics.h"
#include "kernel/importKinetics.h"
#include "kernel/stringUtils.h"

namespace Cantera_CXX {

    /**
     * This class is a convenience class for use in C++ programs that
     * hard-wires the GRI 3.0 reaction mechanism. It derivees from
     * both Cantera::IdealGasPhase, which handles all composition and
     * state information, as well as thermodynamic properties, and
     * class GRI_30_Kinetics, which is the kinetics manager with
     * hard-wired replacements for some of the generic kinetics
     * methods like "getNetReactionRates."
     */
    class GRI30 : 
        public Cantera::IdealGasPhase, 
        public Cantera::GRI_30_Kinetics
    {
    public:
        GRI30() : m_ok(false), m_r(0) {
            m_r = Cantera::get_XML_File("gri30.xml");
            m_ok = Cantera::buildSolutionFromXML(*m_r, "gri30", 
                "phase", this, this);
            if (!m_ok) throw Cantera::CanteraError("GRI30",
                "buildSolutionFromXML returned false");
        }

        virtual ~GRI30() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }
        friend std::ostream& operator<<(std::ostream& s, GRI30& mix) {
            std::string r = Cantera::report(mix, true);
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
