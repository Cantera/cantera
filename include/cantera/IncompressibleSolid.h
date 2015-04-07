//! @file IncompressibleSolid.h
#ifndef CXX_INCOMPRESSIBLE
#define CXX_INCOMPRESSIBLE

#include "thermo/ConstDensityThermo.h"
#include "kinetics/importKinetics.h"

namespace Cantera
{

class IncompressibleSolid : public ConstDensityThermo
{
public:
    IncompressibleSolid(const std::string& infile,
                        std::string id="") : m_ok(false), m_r(0) {

        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }
        m_ok = buildSolutionFromXML(*m_r, id, "phase", this, 0);
        if (!m_ok) throw CanteraError("IncompressibleSolid",
                                          "buildSolutionFromXML returned false");
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }

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
