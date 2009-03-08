/**
 *
 *  @file ThermoPhase.cpp
 */

/*
 *  $Author: hkmoffa $
 *  $Date: 2006/10/20 21:19:35 $
 *  $Revision: 1.15 $
 *
 *  Copyright 2002 California Institute of Technology
 *
 */

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"


namespace Cantera {
    /**
     * Copy Constructor for the ThermoPhase object. 
     *
     * Currently, this is implemented, but not tested. If called it will
     * throw an exception until fully tested.
     */
    ThermoPhase::ThermoPhase(const ThermoPhase &right)  :
	Phase(),
	m_spthermo(0), 
	m_speciesData(0),
	m_index(-1),
	m_phi(0.0), 
	m_hasElementPotentials(false)
    {
	/*
	 * Call the assignment operator
	 */
	*this = operator=(right);
    }

    /*
     * operator=()
     *
     *  Note this stuff will not work until the underlying phase
     *  has a working assignment operator
     */
    ThermoPhase& ThermoPhase::
    operator=(const ThermoPhase &right) {
	/*
         * Check for self assignment.
         */
        if (this == &right) return *this;

	(void)Phase::operator=(right);

	/*
	 * Pointer to the species thermodynamic property manager
	 * We own this, so we need to do a deep copy
	 */
	if (m_spthermo) {
	  delete m_spthermo;
	}
        //m_spthermo = (right.m_spthermo)->duplMyselfAsSpeciesThermo();
	throw CanteraError("ThermoPhase assignment", "not implemented");

        /// Pointer to  the XML tree containing the species
        /// data for this phase. This is used to access data needed to
        /// construct the transport manager and other properties
        /// later in the initialization process.
        m_speciesData = right.m_speciesData;

      
	m_index = right.m_index;
        m_phi = right.m_phi;
	m_lambda = right.m_lambda;
	m_hasElementPotentials = right.m_hasElementPotentials;

	return *this;
    }

    /*
     * Duplication routine for objects which inherit from 
     * ThermoPhase.
     *
     *  This virtual routine can be used to duplicate thermophase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     * 
     *  Currently, this is not fully implemented. If called, an
     *  exception will be called by the ThermoPhase copy constructor.
     */
    ThermoPhase *ThermoPhase::duplMyselfAsThermoPhase() {
	ThermoPhase* tp = new ThermoPhase(*this);
	return tp;
    }

    int ThermoPhase::activityConvention() const {
	return cAC_CONVENTION_MOLAR;
    }

    void ThermoPhase::getActivities(doublereal* a) {
        getActivityConcentrations(a);
        int nsp = nSpecies();
        int k;
        for (k = 0; k < nsp; k++) a[k] /= standardConcentration(k);
    }

    void ThermoPhase::setState_TPX(doublereal t, doublereal p, 
        const doublereal* x) {
        setMoleFractions(x); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TPX(doublereal t, doublereal p, 
        compositionMap& x) {
        setMoleFractionsByName(x); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TPX(doublereal t, doublereal p, 
        const string& x) {
        compositionMap xx;
        int kk = nSpecies();
        for (int k = 0; k < kk; k++) xx[speciesName(k)] = -1.0;
        try {
            parseCompString(x, xx);
        }
        catch (CanteraError) {
            throw CanteraError("setState_TPX",
                "Unknown species in composition map: "+ x);
        }
        setMoleFractionsByName(xx); setTemperature(t); setPressure(p);
    }        

    void ThermoPhase::setState_TPY(doublereal t, doublereal p, 
        const doublereal* y) {
        setMassFractions(y); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TPY(doublereal t, doublereal p, 
        compositionMap& y) {
        setMassFractionsByName(y); setTemperature(t); setPressure(p);
    }
        
    void ThermoPhase::setState_TPY(doublereal t, doublereal p, 
        const string& y) {
        compositionMap yy;
        int kk = nSpecies();
        for (int k = 0; k < kk; k++) yy[speciesName(k)] = -1.0;
        try {
            parseCompString(y, yy);
        }
        catch (CanteraError) {
            throw CanteraError("setState_TPY",
                "Unknown species in composition map: "+ y);
        }
        setMassFractionsByName(yy); setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_TP(doublereal t, doublereal p) {
        setTemperature(t); setPressure(p);
    }

    void ThermoPhase::setState_PX(doublereal p, doublereal* x) {
        setMoleFractions(x); setPressure(p);
    }

    void ThermoPhase::setState_PY(doublereal p, doublereal* y) {
        setMassFractions(y); setPressure(p);
    }

    void ThermoPhase::setState_HP(doublereal h, doublereal p, 
        doublereal tol) {
        doublereal dt;
        setPressure(p);

        // Newton iteration
        for (int n = 0; n < 50; n++) {
            double h0 = enthalpy_mass();
            dt = (h - h0)/cp_mass();
            // limit step size to 100 K
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            setState_TP(temperature() + dt, p);
            if (fabs(dt) < tol) {
                return;
            }
        }
        throw CanteraError("setState_HP","No convergence. dt = " + fp2str(dt));
    }

    void ThermoPhase::setState_UV(doublereal u, doublereal v, 
        doublereal tol) {
        doublereal dt;
        setDensity(1.0/v);
        for (int n = 0; n < 50; n++) {
            dt = (u - intEnergy_mass())/cv_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0;
            if (fabs(dt) < tol) {
                setTemperature(temperature() + dt);
                return;
            }
            setTemperature(temperature() + 0.5*dt);
        }
        throw CanteraError("setState_UV",
            "no convergence. dt = " + fp2str(dt)+"\n"
            +"tol = "+fp2str(tol)+"\n"
            +"u = "+fp2str(u)+" v = "+fp2str(v)+"\n");
    }

    void ThermoPhase::setState_SP(doublereal s, doublereal p, 
        doublereal tol) {
        doublereal dt;
        setPressure(p);
        for (int n = 0; n < 50; n++) {
            dt = (s - entropy_mass())*temperature()/cp_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            if (fabs(dt) < tol) {
                setState_TP(temperature() + dt, p);
                return;
            }
            setState_TP(temperature() + 0.5*dt, p);
        }
        throw CanteraError("setState_SP","no convergence. dt = " + fp2str(dt));
    }

    void ThermoPhase::setState_SV(doublereal s, doublereal v, 
        doublereal tol) {
        doublereal dt;
        setDensity(1.0/v);
        for (int n = 0; n < 50; n++) {
            dt = (s - entropy_mass())*temperature()/cv_mass();
            if (dt > 100.0) dt = 100.0;
            else if (dt < -100.0) dt = -100.0; 
            if (fabs(dt) < tol) {
                setTemperature(temperature() + dt);
                return;
            }
            setTemperature(temperature() + 0.5*dt);
        }
        throw CanteraError("setState_SV","no convergence. dt = " + fp2str(dt));
    }

    doublereal ThermoPhase::err(string msg) const {
            throw CanteraError("ThermoPhase","Base class method "
                +msg+" called. Equation of state type: "+int2str(eosType()));
            return 0;
    }

	/**
	 * Returns the units of the standard and general concentrations
	 * Note they have the same units, as their divisor is 
	 * defined to be equal to the activity of the kth species
	 * in the solution, which is unitless.
	 *
	 * This routine is used in print out applications where the
	 * units are needed. Usually, MKS units are assumed throughout
	 * the program and in the XML input files. 
	 *
	 * On return uA contains the powers of the units (MKS assumed)
	 * of the standard concentrations and generalized concentrations
	 * for the kth species.
	 *
	 *  uA[0] = kmol units - default  = 1
	 *  uA[1] = m    units - default  = -nDim(), the number of spatial
	 *                                dimensions in the Phase class.
	 *  uA[2] = kg   units - default  = 0;
	 *  uA[3] = Pa(pressure) units - default = 0;
	 *  uA[4] = Temperature units - default = 0;
	 *  uA[5] = time units - default = 0
	 */
    void ThermoPhase::getUnitsStandardConc(double *uA, int k, int sizeUA) {
	for (int i = 0; i < sizeUA; i++) {
	  if (i == 0) uA[0] = 1.0;
	  if (i == 1) uA[1] = -nDim();
	  if (i == 2) uA[2] = 0.0;
	  if (i == 3) uA[3] = 0.0;
	  if (i == 4) uA[4] = 0.0;
	  if (i == 5) uA[5] = 0.0;
	}
    }

    /*
     * initThermoFile():
     *
     * Initialization of a Debye-Huckel phase using an
     * xml file.
     *
     * This routine is a precursor to initThermoXML(XML_Node*)
     * routine, which does most of the work. 
     *
     * @param infile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void ThermoPhase::initThermoFile(string inputFile, string id) {

	if (inputFile.size() == 0) {
	  throw CanteraError("ThermoPhase::initThermoFile",
			     "input file is null");
	}
	string path = findInputFile(inputFile);
	ifstream fin(path.c_str());
	if (!fin) {
	  throw CanteraError("initThermoFile","could not open "
			     +path+" for reading.");
	}
	/*
	 * The phase object automatically constructs an XML object.
	 * Use this object to store information.
	 */
	XML_Node &phaseNode_XML = xml();
        XML_Node *fxml = new XML_Node();
	fxml->build(fin);
	XML_Node *fxml_phase = findXMLPhase(fxml, id);
	if (!fxml_phase) {
	  throw CanteraError("ThermoPhase::initThermo",
			     "ERROR: Can not find phase named " +
			     id + " in file named " + inputFile);
	}
	fxml_phase->copy(&phaseNode_XML);	
	initThermoXML(*fxml_phase, id);
	delete fxml;
    }

    /*
     *   Import and initialize a ThermoPhase
     *   object
     *
     * @param phaseNode This object must be the phase node of a
     *             complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id. 
     */
    void ThermoPhase::initThermoXML(XML_Node& phaseNode, string id) {
	/*
	 * The default implementation just calls initThermo();
	 */
	initThermo();
	/*
	 * and sets the state
	 */
	if (phaseNode.hasChild("state")) {
	  XML_Node& stateNode = phaseNode.child("state");
	  setStateFromXML(stateNode);
	}
    }

    /*
     * Initialize. 
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase.
     *
     * @see importCTML.cpp
     */
    void ThermoPhase::initThermo() {

    }
    
    /**
     * Set the thermodynamic state.
     */
    void ThermoPhase::setStateFromXML(const XML_Node& state) {

        string comp = getString(state,"moleFractions");
        if (comp != "") 
            setMoleFractionsByName(comp);
        else {
            comp = getString(state,"massFractions");
            if (comp != "") 
                setMassFractionsByName(comp);
        }
        if (state.hasChild("temperature")) {
            double t = getFloat(state, "temperature", "temperature");
            setTemperature(t);
        }
        if (state.hasChild("pressure")) {
            double p = getFloat(state, "pressure", "pressure");
            setPressure(p);
        }
        if (state.hasChild("density")) {
            double rho = getFloat(state, "density", "density");
            setDensity(rho);
        }
    }

}
