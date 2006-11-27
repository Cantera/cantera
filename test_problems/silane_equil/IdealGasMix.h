#ifndef CXX_IDEALGASMIX
#define CXX_IDEALGASMIX

#include <string>

#include "IdealGasPhase.h"
#include "GasKinetics.h"
#include "importCTML.h"

namespace Cantera {

    class IdealGasMix : 
        public IdealGasPhase, public GasKinetics
    {
    public:

        IdealGasMix() : m_ok(false), m_r(0) {}

        IdealGasMix(std::string infile, std::string id="") : m_ok(false), m_r(0) {
            
        m_r = get_XML_File(infile); 
        if (id == "-") id = "";
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, this);
        if (!m_ok) throw CanteraError("IdealGasMix",
            "buildSolutionFromXML returned false");
        }


        IdealGasMix(XML_Node& root, std::string id) : m_ok(false), m_r(0) {
            m_ok = buildSolutionFromXML(root, id, "phase", this, this);
        }
        
        virtual ~IdealGasMix() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }
        friend std::ostream& operator<<(std::ostream& s, IdealGasMix& mix) {
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
