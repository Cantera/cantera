#ifndef CXX_IDEALGASMIX
#define CXX_IDEALGASMIX

#include <string>

#include "kernel/IdealGasPhase.h"
#include "kernel/GasKinetics.h"
#include "kernel/importKinetics.h"
#include "kernel/stringUtils.h"

namespace Cantera_CXX {

    class IdealGasMix : 
        public Cantera::IdealGasPhase, 
        public Cantera::GasKinetics
    {
    public:

        IdealGasMix() : m_ok(false), m_r(0) {}

        IdealGasMix(std::string infile, std::string id="") : 
            m_ok(false), m_r(0) {
            
            m_r = Cantera::get_XML_File(infile); 
            if (id == "-") id = "";
            m_ok = Cantera::buildSolutionFromXML(*m_r, 
                id, "phase", this, this);
            if (!m_ok) throw Cantera::CanteraError("IdealGasMix",
                "Cantera::buildSolutionFromXML returned false");
        }


        IdealGasMix(Cantera::XML_Node& root, 
            std::string id) : m_ok(false), m_r(0) {
            m_ok = Cantera::buildSolutionFromXML(root, id, 
                "phase", this, this);
        }
        
        virtual ~IdealGasMix() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }
        friend std::ostream& operator<<(std::ostream& s, IdealGasMix& mix) {
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
