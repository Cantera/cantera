//! @file Edge.h
#ifndef CXX_EDGE
#define CXX_EDGE

#include "thermo.h"
#include "kinetics/EdgeKinetics.h"
#include "kinetics/importKinetics.h"

namespace Cantera
{

class Edge :
    public EdgePhase, public EdgeKinetics
{
public:
    Edge(const std::string& infile, std::string id, std::vector<ThermoPhase*> phases)
        : m_ok(false), m_r(0) {

        m_r = get_XML_File(infile);
        if (id == "-") {
            id = "";
        }

        XML_Node* x = get_XML_Node("#"+id, m_r);
        if (!x) {
            throw CanteraError("Edge","error in get_XML_Node");
        }

        importPhase(*x, this);
        phases.push_back(this);
        importKinetics(*x, phases, this);
        m_ok = true;
    }

    bool operator!() {
        return !m_ok;
    }
    bool ready() const {
        return m_ok;
    }

protected:
    bool m_ok;
    XML_Node* m_r;

private:
};
}


#endif
