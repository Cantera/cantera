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
            string path = findInputFile(infile);
            ifstream fin(path.c_str());
            if (!fin) {
                throw CanteraError("Interface","could not open "
                    +path+" for reading.");
            }
            
            m_r = new XML_Node("-");
            m_r->build(fin);

            XML_Node* x = find_XML("", m_r, id, "", "");
            if (!x)                 
                throw CanteraError("Interface","error in find_XML");

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
}


#endif
