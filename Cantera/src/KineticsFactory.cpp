/**
 *  @file KineticsFactory.cpp
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "KineticsFactory.h"

#include "GasKinetics.h"
#include "GRI_30_Kinetics.h"
#include "InterfaceKinetics.h"
#include "EdgeKinetics.h"
#include "importCTML.h"

namespace Cantera {


    //---------- class Kinetics methods ------------------
    string Kinetics::kineticsSpeciesName(int k) const {
        int np = m_start.size();
        for (int n = np-1; n >= 0; n--) {
            if (k >= m_start[n]) {
                return thermo(n).speciesName(k - m_start[n]);
            }
        }
        return "<unknown>";
    }


    int Kinetics::kineticsSpeciesIndex(string nm, string ph) const {
        int np = static_cast<int>(m_thermo.size());
        int k;
        string id;
        for (int n = 0; n < np; n++) {
	    id = thermo(n).id();
	    if (ph == id) {
                k = thermo(n).speciesIndex(nm);
                if (k < 0) return -1;
                return k + m_start[n];
	    }
	    else if (ph == "<any>") {
                /*
                 * Call the speciesIndex() member function of the
                 * ThermoPhase object to find a match.
                 */
                k = thermo(n).speciesIndex(nm);
                if (k >= 0) return k + m_start[n];
	    }                    
        }
        return -2;
    }

    thermo_t& Kinetics::speciesPhase(string nm) {
        int np = static_cast<int>(m_thermo.size());
        int k;
        string id;
        for (int n = 0; n < np; n++) {
            k = thermo(n).speciesIndex(nm);
            if (k >= 0) return thermo(n);
        }
        throw CanteraError("speciesPhase", "unknown species "+nm);
    }

    int Kinetics::speciesPhaseIndex(int k) {
        int np = m_start.size();
        for (int n = np-1; n >= 0; n--) {
            if (k >= m_start[n]) {
                return n;
            }
        }
        throw CanteraError("speciesPhaseIndex", 
            "illegal species index: "+int2str(k));
    }

    void Kinetics::addPhase(thermo_t& thermo) {

            // if not the first thermo object, set the start position
            // to that of the last object added + the number of its species 
            if (m_thermo.size() > 0) {
                m_start.push_back(m_start.back() 
				  + m_thermo.back()->nSpecies());
            }
            // otherwise start at 0
            else {
                m_start.push_back(0);
            }
            // there should only be one surface phase
            int ptype = -100;
            if (type() == cEdgeKinetics) ptype = cEdge;
            else if (type() == cInterfaceKinetics) ptype = cSurf;
            if (thermo.eosType() == ptype) {
                if (m_surfphase >= 0) {
                    throw CanteraError("Kinetics::addPhase",
                        "cannot add more than one surface phase");
                }
                m_surfphase = nPhases();
            }
            m_thermo.push_back(&thermo);
            m_phaseindex[m_thermo.back()->id()] = nPhases();
        }

    void Kinetics::err(string m) const {
            throw CanteraError("Kinetics::" + m, 
			       "The default Base class method was called, when "
		               "the inherited class's method should have been called");
    }

    //-----------------------------------------------------

    KineticsFactory* KineticsFactory::__factory = 0;

    static int ntypes = 5;
    static string _types[] = {"none", "GasKinetics", "GRI30", "Interface", "Edge"};
    static int _itypes[]   = {0, cGasKinetics, cGRI30, cInterfaceKinetics, cEdgeKinetics};

    /**
     * Return a new kinetics manager that implements a reaction
     * mechanism specified in a CTML file. In other words, the
     * kinetics manager, given the rate constants and formulation of the
     * reactions that make up a kinetics mechanism, is responsible for
     * calculating the rates of progress of the reactions and for
     * calculating the source terms for species.
     * 
     * Input
     * ------
     *  phaseData = This is an XML_Node that contains the xml data
     *              describing the phase. Of particular note to ths
     *              routine is the child xml element called "kinetics".
     *              The element has one attribute called "model",
     *              with a string value. The value of this string
     *              is used to decide which kinetics manager is used
     *              to calculate the reacton mechanism. 
     *
     * Return
     * ---------
     *  Pointer to the new kinetics manager. 
     */

    Kinetics* KineticsFactory::
    newKinetics(XML_Node& phaseData, vector<ThermoPhase*> th) {
	/*
	 * Look for a child of the xml element phase called
	 * "kinetics". It has an attribute name "model".
	 * Store the value of that attribute in the variable kintype
	 */
        string kintype = phaseData.child("kinetics")["model"];
	/*
	 * look up the string kintype in the list of known
	 * kinetics managers (list is kept at the top of this file).
	 * Translate it to an integer value, ikin. 
	 */
        int ikin=-1;
        int n;
        for (n = 0; n < ntypes; n++) {
	  if (kintype == _types[n]) ikin = _itypes[n];
        }
	/*
	 * Assign the kinetics manager based on the value of ikin.
	 * Kinetics managers are classes derived from the base 
	 * Kinetics class. Unknown kinetics managers will throw a
	 * CanteraError here.
	 */
        Kinetics* k=0;
        switch (ikin) {
            
        case 0:
            k = new Kinetics;
            break;

        case cGasKinetics:
            k = new GasKinetics;
            break;

        case cGRI30:
            k = new GRI_30_Kinetics;
            break;

        case cInterfaceKinetics:
            k = new InterfaceKinetics;
            break;

        case cEdgeKinetics:
            k = new EdgeKinetics;
            break;

        default:
	    throw UnknownKineticsModel("KineticsFactory::newKinetics", 
				       kintype);
        }

        // Now that we have the kinetics manager, we can
        // import the reaction mechanism into it.
        importKinetics(phaseData, th, k);

        // Return the pointer to the kinetics manager
        return k;
    }


    /**
     * Return a new, empty kinetics manager.
     */
    Kinetics* KineticsFactory::newKinetics(string model) {

        int ikin = -1;
        int n;
        for (n = 0; n < ntypes; n++) {
            if (model == _types[n]) ikin = _itypes[n];
        }
        Kinetics* k=0;
        switch (ikin) {

        case cGasKinetics:
            k = new GasKinetics;
            break;

        case cGRI30:
            k = new GRI_30_Kinetics;
            break;

        case cInterfaceKinetics:
            k = new InterfaceKinetics;
            break;

        default:
	    throw UnknownKineticsModel("KineticsFactory::newKinetics", 
				       model);
        }
        return k;
    }

}
