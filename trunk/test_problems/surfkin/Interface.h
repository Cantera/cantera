/*
 * $Id: Interface.h,v 1.4 2006/04/30 21:43:56 hkmoffa Exp $
 */
#ifndef CXX_INTERFACE_H
#define CXX_INTERFACE_H

#include <string>
#ifdef SRCDIRTREE
#include "SurfPhase.h"
#include "InterfaceKinetics.h"
#include "importCTML.h"
#else
#include "kernel/SurfPhase.h"
#include "kernel/InterfaceKinetics.h"
#include "kernel/importCTML.h"
#endif

namespace Cantera {

    /**
     * The class interface inherits from both SurfPhase and
     * InterferFaceKinetics
     */
    class Interface : 
        public SurfPhase, public InterfaceKinetics
    {
    public:
	/**
	 *
	 * Constructor for the interface class:
	 *
	 *  infile = name of the file to get information about the
	 *           surface.
	 *  id = name of the surface phase.
	 *  phases = Pointer to the list of volume phases that participate
	 *           in the interface
	 *
	 */
        Interface(string infile, string id, vector<ThermoPhase*> phases) 
            : m_ok(false), m_r(0) {
            string path = findInputFile(infile);
            ifstream fin(path.c_str());
            if (!fin) {
                throw CanteraError("Interface","could not open "
                    +path+" for reading.");
            }
   

	    /*
	     * Create a top level xml node
	     */
            m_r = new XML_Node("-");
	    /*
	     * Fill the XML_Node with all of the information in the
	     * xml file
	     */
            m_r->build(fin);


	    /*
	     * Find the start of the surface phase data in the xml file
	     * Store a pointer to the position in the tree.
	     */
            XML_Node* x = get_XML_Node("#" + id, m_r);
            if (!x) {
	      throw CanteraError("Interface","error in find_XML");
	    }

	    /*
	     * Import the values of the surface into the current object.
	     * Note, since the current object inherits from SurfPhase
	     * object, it contains all of the surfphase fields of the 
	     * object. This operation fills all of those fields.
	     */
            importPhase(*x, this);
            
            phases.push_back(this);
            importKinetics(*x, phases, this);
            m_ok = true;
        }

	/**
	 * Destructor for the Interface class. This is a virtual function, meaning
	 * that after, it finishes with this class, it will call the destructor
	 * functions for the two classes that Interface inherits from.
	 */
        virtual ~Interface() {
	  /*
	   * We created the XML data tree structure for the Interface object
	   * in the constructor. Here, we delete the top level. of the structure.
	   * I believe this delete all of the daughter elements as well.
	   */
	  delete m_r;
	}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }

    protected:
        bool m_ok;
        XML_Node* m_r;

    private:
    };
}


#endif
