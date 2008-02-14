#ifndef CXX_INCOMPRESSIBLE
#define CXX_INCOMPRESSIBLE

#include <string>

#include "kernel/ConstDensityThermo.h"
#include "kernel/importKinetics.h"

namespace Cantera_CXX {

    class IncompressibleSolid : public Cantera::ConstDensityThermo
    {
    public:
        IncompressibleSolid(std::string infile, 
            std::string id="") : m_ok(false), m_r(0) {
            
            m_r = Cantera::get_XML_File(infile); 
            if (id == "-") id = "";
            m_ok = Cantera::buildSolutionFromXML(*m_r, id, "phase", this, 0);
            if (!m_ok) throw Cantera::CanteraError("IncompressibleSolid",
                "buildSolutionFromXML returned false");
        }


        virtual ~IncompressibleSolid() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }

        //friend std::ostream& operator<<(std::ostream& s, IdealGasMix& mix) {
        //    std::string r = report(mix, true);
        //    s << r;
        //    return s;

    protected:
        bool m_ok;
        Cantera::XML_Node* m_r;

    private:
    };
}


#endif
