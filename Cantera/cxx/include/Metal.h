#ifndef CXX_METAL
#define CXX_METAL

#include <string>

#include "kernel/MetalPhase.h"
#include "kernel/importKinetics.h"

namespace Cantera {

    class Metal : public MetalPhase
    {
    public:
        Metal(std::string infile, std::string id="") : m_ok(false), m_r(0) {
            
        m_r = get_XML_File(infile); 
        if (id == "-") id = "";
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("Metal",
            "buildSolutionFromXML returned false");
        }


        virtual ~Metal() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }

        //friend std::ostream& operator<<(std::ostream& s, IdealGasMix& mix) {
        //    std::string r = report(mix, true);
        //    s << r;
        //    return s;

    protected:
        bool m_ok;
        XML_Node* m_r;

    private:
    };
}


#endif
