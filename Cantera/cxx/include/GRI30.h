/*
 * $Id$
 */
#ifndef CXX_GRI30H
#define CXX_GRI30H

#include <string>

#include "kernel/IdealGasPhase.h"
#include "kernel/GRI_30_Kinetics.h"
#include "kernel/importKinetics.h"

namespace Cantera {

    class GRI30 : 
        public IdealGasPhase, public GRI_30_Kinetics
    {
    public:
        GRI30() : m_ok(false), m_r(0) {
            m_r = get_XML_File("gri30.xml");
            m_ok = buildSolutionFromXML(*m_r, "gri30", "phase", this, this);
            if (!m_ok) throw CanteraError("GRI30",
                "buildSolutionFromXML returned false");
        }


        virtual ~GRI30() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }
        friend std::ostream& operator<<(std::ostream& s, GRI30& mix) {
            std::string r = report(mix, true);
            s << r;
            return s;
        }

    protected:
        bool m_ok;
        XML_Node* m_r;

    private:
    };
}


#endif
