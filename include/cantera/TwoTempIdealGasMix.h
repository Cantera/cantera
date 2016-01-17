//! @file TwoTempGasMix.h
#ifndef CXX_TWOTEMPIDEALGASMIX
#define CXX_TWOTEMPIDEALGASMIX

#include "thermo/IdealGasPhase.h"
#include "kinetics/TwoTempGasKinetics.h"
#include "kinetics/importKinetics.h"
#include "base/stringUtils.h"

namespace Cantera
{

class TwoTempIdealGasMix :
    public IdealGasPhase,
    public TwoTempGasKinetics
{
public:

    TwoTempIdealGasMix() : m_ok(false), m_r(0) {}

    TwoTempIdealGasMix(const std::string& infile, std::string id_="") :
        m_ok(false), m_r(0) {

        m_r = get_XML_File(infile);
        m_id = id_;
        if (id_ == "-") {
            id_ = "";
        }
        m_ok = buildSolutionFromXML(*m_r,
                                    m_id, "phase", this, this);
        if (!m_ok) throw CanteraError("TwoTempIdealGasMix",
                                          "Cantera::buildSolutionFromXML returned false");
    }


    TwoTempIdealGasMix(XML_Node& root,
                std::string id_) : m_ok(false), m_r(&root), m_id(id_) {
        m_ok = buildSolutionFromXML(root, id_, "phase", this, this);
    }

    TwoTempIdealGasMix(const TwoTempIdealGasMix& other) : m_ok(false),
        m_r(other.m_r),
        m_id(other.m_id) {
        m_ok = buildSolutionFromXML(*m_r, m_id, "phase", this, this);
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }
    friend std::ostream& operator<<(std::ostream& s, TwoTempIdealGasMix& mix) {
        std::string r = mix.report(true);
        s << r;
        return s;
    }


protected:
    bool m_ok;
    XML_Node* m_r;
    std::string m_id;

private:
};
}


#endif
