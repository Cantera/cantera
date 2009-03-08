#ifndef CXX_IDEALGASMIX
#define CXX_IDEALGASMIX

#include <string>

#include "kernel/IdealGasPhase.h"
#include "kernel/GasKinetics.h"
#include "kernel/importCTML.h"

namespace Cantera {

    class IdealGasMix : 
        public IdealGasPhase, public GasKinetics
    {
    public:

        IdealGasMix() : m_ok(false), m_r(0) {}

        IdealGasMix(string infile, string id="") : m_ok(false), m_r(0) {
            
        m_r = get_XML_File(infile); 
        if (id == "-") id = "";
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, this);
        if (!m_ok) throw CanteraError("IdealGasMix",
            "buildSolutionFromXML returned false");
        }


        IdealGasMix(XML_Node& root, string id) : m_ok(false), m_r(0) {
            m_ok = buildSolutionFromXML(root, id, "phase", this, this);
        }
        
        virtual ~IdealGasMix() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }
        friend ostream& operator<<(ostream& s, IdealGasMix& mix) {
            string r = report(mix, true);
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
