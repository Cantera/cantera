#ifndef CXX_INTERFACE
#define CXX_INTERFACE

#include <string>

#include "kernel/SurfPhase.h"
#include "kernel/InterfaceKinetics.h"
#include "kernel/importCTML.h"

namespace Cantera {

    class Interface : 
        public SurfPhase, public InterfaceKinetics
    {
    public:
        Interface(string infile, string id, vector<ThermoPhase*> phases) 
            : m_ok(false), m_r(0) {

            m_r = get_XML_File(infile); 
            if (id == "-") id = "";

            XML_Node* x = get_XML_Node("#"+id, m_r);
            if (!x)                 
                throw CanteraError("Interface","error in get_XML_Node");

            importPhase(*x, this);
            phases.push_back(this);
            importKinetics(*x, phases, this);
            m_ok = true;
        }

        virtual ~Interface() {}

        bool operator!() { return !m_ok;}
        bool ready() { return m_ok; }

    protected:
        bool m_ok;
        XML_Node* m_r;

    private:
    };

    inline Interface* importInterface(string infile, string id, 
        vector<ThermoPhase*> phases) {
        return new Interface(infile, id, phases);
    }

}


#endif
