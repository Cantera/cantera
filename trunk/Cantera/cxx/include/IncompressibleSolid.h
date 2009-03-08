#ifndef CXX_INCOMPRESSIBLE
#define CXX_INCOMPRESSIBLE

#include <string>

#include "kernel/ConstDensityThermo.h"
#include "kernel/importCTML.h"

namespace Cantera {

    class IncompressibleSolid : public ConstDensityThermo
    {
    public:
        IncompressibleSolid(string infile, string id="") : m_ok(false), m_r(0) {
            
        m_r = get_XML_File(infile); 
        if (id == "-") id = "";
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("IncompressibleSolid",
            "buildSolutionFromXML returned false");
        }


        virtual ~IncompressibleSolid() {}

        bool operator!() { return !m_ok;}
        bool ready() { return m_ok; }

        //friend ostream& operator<<(ostream& s, IdealGasMix& mix) {
        //    string r = report(mix, true);
        //    s << r;
        //    return s;

    protected:
        bool m_ok;
        XML_Node* m_r;

    private:
    };
}


#endif
