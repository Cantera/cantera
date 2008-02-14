#ifndef CXX_METAL
#define CXX_METAL

#include <string>

#include "kernel/MetalPhase.h"
#include "kernel/importKinetics.h"

namespace Cantera_CXX {

    class Metal : public Cantera::MetalPhase
    {
    public:
        Metal(std::string infile, std::string id="") : m_ok(false), m_r(0) {
            
            m_r = Cantera::get_XML_File(infile); 
            if (id == "-") id = "";
            m_ok = Cantera::buildSolutionFromXML(*m_r, id, "phase", this, 0);
            if (!m_ok) throw Cantera::CanteraError("Metal",
                "buildSolutionFromXML returned false");
        }
        
        virtual ~Metal() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }

    protected:
        bool m_ok;
        Cantera::XML_Node* m_r;

    private:
    };
}


#endif
