/**
 * @file Interface.h
 *   Declaration and Definition for the class Interface, part of 
 *   Cantera's Cantera_CXX namespace.
 */

/*
 * $Id: Interface.h,v 1.10 2008/01/09 13:43:26 dggoodwin Exp $
 */

#ifndef CXX_INTERFACE
#define CXX_INTERFACE

#include <string>
#include "Cantera.h"
#include "thermo.h"
#include "kinetics.h"
// #include "kernel/SurfPhase.h"
// #include "kernel/InterfaceKinetics.h"
// #include "kernel/importKinetics.h"

/**
 * This namespace is used for the Cantera C++ user interface.
 */
namespace Cantera_CXX {

    /**
     * An interface between multiple bulk phases. This class is
     * defined mostly for convenience. It inherits both from
     * Cantera::SurfPhase and Cantera::InterfaceKinetics. It therefore
     * represents a surface phase, and also acts as the kinetics
     * manager to manage reaction occurring on the surface, possibly
     * involving species from other phases.
     */
    class Interface : 
        public Cantera::SurfPhase,
        public Cantera::InterfaceKinetics
    {
    public:

        /**
         * Constructor. Construct an Interface instance from
         * a specification in an input file.
         * @param infile.  Cantera input file in CTI or CTML format.
         * @param id  Identification string to distinguish between
         * multiple definitions within one input file.
         * @param phases Neighboring phases that may participate in the
         * reactions on this interface.
         */
        Interface(std::string infile, std::string id, 
            std::vector<Cantera::ThermoPhase*> phases) 
            : m_ok(false), m_r(0) {
            m_r = Cantera::get_XML_File(infile); 
            if (id == "-") id = "";
            
            Cantera::XML_Node* x = Cantera::get_XML_Node("#"+id, m_r);
            if (!x)                 
                throw Cantera::CanteraError("Interface","error in get_XML_Node");

            Cantera::importPhase(*x, this);
            phases.push_back(this);
            Cantera::importKinetics(*x, phases, this);
            m_ok = true;
        }

        /// Destructor. Does nothing.
        virtual ~Interface() {}

        bool operator!() { return !m_ok;}
        bool ready() const { return m_ok; }

    protected:
        bool m_ok;
        Cantera::XML_Node* m_r;

    private:
    };

    /**
     * Import an instance of class Interface from a specification in an 
     * input file. This is the preferred method to create an Interface
     * instance.
     */
    inline Interface* importInterface(std::string infile, std::string id, 
        std::vector<Cantera::ThermoPhase*> phases) {
        return new Interface(infile, id, phases);
    }

}


#endif
