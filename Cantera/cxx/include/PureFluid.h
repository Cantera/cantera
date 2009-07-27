/*
 * $Id: PureFluid.h,v 1.7 2007/12/21 02:58:41 hkmoffa Exp $
 */
#ifndef CXX_PUREFLUID
#define CXX_PUREFLUID

#include <string>

#include "kernel/PureFluidPhase.h"
#include "kinetics.h"
#include "kernel/stringUtils.h"

namespace Cantera {


    class PureFluid : public PureFluidPhase
    {
    public:

        PureFluid() : m_ok(false), m_r(0) {}

        PureFluid(std::string infile, std::string id="") : m_ok(false), m_r(0) {
            
        m_r = get_XML_File(infile); 
        if (id == "-") id = "";
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("PureFluid",
            "buildSolutionFromXML returned false");
        }


        PureFluid(XML_Node& root, std::string id) : m_ok(false), m_r(0) {
            m_ok = buildSolutionFromXML(root, id, "phase", this, 0);
        }
        
        virtual ~PureFluid() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }
        friend std::ostream& operator<<(std::ostream& s, PureFluid& mix) {
            std::string r = Cantera::report(mix, true);
            s << r;
            return s;
        }

    protected:
        bool m_ok;
        XML_Node* m_r;

    private:
    };

    class Water : public PureFluid {
    public:
        Water() : PureFluid(std::string("liquidvapor.cti"),std::string("water")) {}
        virtual ~Water() {}
    };

}


#endif
