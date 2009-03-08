/**
 *  @file WaterTP.cpp
 *
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: WaterTP.cpp,v 1.1 2006/07/04 00:01:54 hkmoffa Exp $
 */

#include "xml.h"
#include "WaterTP.h"
#include "WaterPropsIAPWS.h"
#include "importCTML.h"

namespace Cantera {
    /**
     * Basic list of constructors and duplicators
     */

    WaterTP::WaterTP() :
	ThermoPhase(),
	m_sub(0),
	m_subflag(0),
	m_mw(0.0),
	EW_Offset(0.0),
	SW_Offset(0.0),
	m_verbose(0),
	m_allowGasPhase(false)
    {
	constructPhase();
    }


    WaterTP::WaterTP(string inputFile, string id) :
	ThermoPhase(),
	m_sub(0),
	m_subflag(0),
	m_mw(0.0),
	EW_Offset(0.0),
	SW_Offset(0.0),
	m_verbose(0),
	m_allowGasPhase(false)
    {
	constructPhaseFile(inputFile, id);
    }


    WaterTP::WaterTP(XML_Node& phaseRoot, string id) :
	ThermoPhase(),
	m_sub(0),
	m_subflag(0),
	m_mw(0.0),
	EW_Offset(0.0),
	SW_Offset(0.0),
	m_verbose(0),
	m_allowGasPhase(false)
    {
	constructPhaseXML(phaseRoot, id) ;
    }



   WaterTP::WaterTP(const WaterTP &b) :
       ThermoPhase(b),
       m_sub(0),
       m_subflag(b.m_subflag),
       m_mw(b.m_mw),
       EW_Offset(b.EW_Offset),
       SW_Offset(b.SW_Offset),
       m_verbose(b.m_verbose),
       m_allowGasPhase(b.m_allowGasPhase)
    {
        m_sub = new WaterPropsIAPWS(*(b.m_sub));  
	/*
         * Use the assignment operator to do the brunt
         * of the work for the copy construtor.
         */
        *this = b;
    }

    /**
     * Assignment operator
     */
    WaterTP& WaterTP::operator=(const WaterTP&b) {
	if (&b == this) return *this;
	m_sub->operator=(*(b.m_sub));
	m_subflag = b.m_subflag;
	m_mw = b.m_mw;
	m_verbose = b.m_verbose;
	m_allowGasPhase = b.m_allowGasPhase;
	return *this;
    }


    ThermoPhase *WaterTP::duplMyselfAsThermoPhase() {
	WaterTP* wtp = new WaterTP(*this);
	return (ThermoPhase *) wtp;
    }

    WaterTP::~WaterTP() { 
	delete m_sub; 
    }


  
    void WaterTP::constructPhase() {
	throw CanteraError("constructPhaseXML", "unimplemented");

    }

   
    /**
     * constructPhase:
     *
     * Initialization of a Debye-Huckel phase using an
     * xml file.
     *
     * This routine is a precursor to initThermo(XML_Node*)
     * routine, which does most of the work.
     *
     * @param infile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void WaterTP::constructPhaseXML(XML_Node& phaseNode, string id) {

	/*
         * Call the Cantera importPhase() function. This will import
         * all of the species into the phase. This will also handle
         * all of the solvent and solute standard states.
         */
        bool m_ok = importPhase(phaseNode, this);
        if (!m_ok) {
          throw CanteraError("initThermo","importPhase failed ");
        }
	
    }



   
    /**
     * initThermo():
     *
     * Initialization of a Debye-Huckel phase using an
     * xml file.
     *
     * This routine is a precursor to initThermo(XML_Node*)
     * routine, which does most of the work.
     *
     * @param infile XML file containing the description of the
     *        phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void WaterTP::constructPhaseFile(string inputFile, string id) {

	if (inputFile.size() == 0) {
	  throw CanteraError("WaterTp::initThermo",
			     "input file is null");
	}
	string path = findInputFile(inputFile);
	ifstream fin(path.c_str());
	if (!fin) {
	  throw CanteraError("WaterTP::initThermo","could not open "
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
	  throw CanteraError("WaterTP::initThermo",
			     "ERROR: Can not find phase named " +
			     id + " in file named " + inputFile);
	}
	fxml_phase->copy(&phaseNode_XML);	
	constructPhaseXML(*fxml_phase, id);
	delete fxml;
    }



    void WaterTP::initThermo() {




    }

    void WaterTP::
    initThermoXML(XML_Node& phaseNode, string id) {
        if (m_sub) delete m_sub;
        m_sub = new WaterPropsIAPWS();
        if (m_sub == 0) {
            throw CanteraError("WaterTP::initThermo",
                "could not create new substance object.");
        }
	/*
	 * Calculate the molecular weight. Note while there may
	 * be a very good calculated weight in the steam table
	 * class, using this weight may lead to codes exhibiting
	 * mass loss issues. We need to grab the elemental
	 * atomic weights used in the Element class and calculate
	 * a consistent H2O molecular weight based on that.
	 */
	int nH = elementIndex("H");
	if (nH < 0) {
	  throw CanteraError("WaterTP::initThermo",
			     "H not an element");
	}
	double mw_H = atomicWeight(nH);
	int nO = elementIndex("O");
	if (nO < 0) {
	  throw CanteraError("WaterTP::initThermo",
			     "O not an element");
	}
	double mw_O = atomicWeight(nO);
        m_mw = 2.0 * mw_H + mw_O;
        m_weight[0] = m_mw;
        setMolecularWeight(0,m_mw);
        double one = 1.0;
        setMoleFractions(&one);

     	/*
	 * Set the baseline 
	 */
	doublereal T = 298.15;

	doublereal presLow = 1.0E-2;
	doublereal oneBar = 1.0E5;
	doublereal dens = density();
	doublereal dd = m_sub->density(T, presLow, WATER_GAS, dens);
	setTemperature(T);
	setDensity(dd);
	SW_Offset = 0.0;
	doublereal s = entropy_mole();
	s -=  GasConstant * log(oneBar/presLow);
	if (s != 188.835E3) {
	  SW_Offset = 188.835E3 - s;
	}
	s = entropy_mole();
	s -=  GasConstant * log(oneBar/presLow);
	printf("s = %g\n", s);

	doublereal h = enthalpy_mole();
	if (h != -241.826E6) {
	  EW_Offset = -241.826E6 - h;
	}
	h = enthalpy_mole();

	printf("h = %g\n", h);


	/*
	 * Set the initial state of the system to 298.15 K and 
	 * 1 bar.
	 */
        setTemperature(298.15);
	double rho0 = m_sub->density(298.15, OneAtm, WATER_LIQUID);
	setDensity(rho0);

	/*
	 * We have to do something with the thermo function here.
	 */
	if (m_spthermo) {
	  delete m_spthermo;
	  m_spthermo = 0;
	}
    }

    void WaterTP::
    setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","PureFluid");
        m_subflag = atoi(eosdata["fluid_type"].c_str());
        if (m_subflag < 0) 
            throw CanteraError("WaterTP::setParametersFromXML",
                "missing or negative substance flag");
    }

    /**
     * Return the molar enthalpy in units of J kmol-1
     */
    doublereal WaterTP::
    enthalpy_mole() const {
	double T = temperature();
	double dens = density();
        doublereal h = m_sub->enthalpy(T, dens);
        return (h + EW_Offset);
    }

    /**
     * Calculate the internal energy in mks units of
     * J kmol-1 
     */
    doublereal WaterTP::
    intEnergy_mole() const {
	double T = temperature();
	double dens = density();
        doublereal u = m_sub->intEnergy(T, dens);
        return (u + EW_Offset);            
    }

    /**
     * Calculate the entropy in mks units of 
     * J kmol-1 K-1
     */
    doublereal WaterTP::
    entropy_mole() const {
	double T = temperature();
	double dens = density();
        doublereal s = m_sub->entropy(T, dens);
        return (s + SW_Offset); 
    }

    /**
     * Calculate the Gibbs free energy in mks units of
     * J kmol-1 K-1.
     */
    doublereal WaterTP::
    gibbs_mole() const {
	double T = temperature();
	double dens = density();
        doublereal g = m_sub->Gibbs(T, dens);
        return (g + EW_Offset - SW_Offset*T);
    }

    /**
     * Calculate the constant pressure heat capacity
     * in mks units of J kmol-1 K-1
     */
    doublereal WaterTP::
    cp_mole() const {
	double T = temperature();
	double dens = density();
        doublereal cp = m_sub->cp(T, dens);
        return cp;            
    }

    /**
     * Calculate the constant volume heat capacity
     * in mks units of J kmol-1 K-1
     */
    doublereal WaterTP::
    cv_mole() const {
	double T = temperature();
	double dens = density();
        doublereal cv = m_sub->cv(T, dens);
        return cv;
    }

    /**
     * Calculate the pressure (Pascals), given the temperature and density
     *  Temperature: kelvin
     *  rho: density in kg m-3
     */
    doublereal WaterTP::
    pressure() const {
	double T = temperature();
	double dens = density();
        doublereal p = m_sub->pressure(T, dens);
        return p;
    }
        
    void WaterTP::
    setPressure(doublereal p) {
	double T = temperature();
	double dens = density();
	int waterState = WATER_GAS;
	double rc = m_sub->Rhocrit();
	if (dens > rc) {
	  waterState = WATER_LIQUID;
	}
	doublereal dd = m_sub->density(T, p, waterState, dens);
	if (dd <= 0.0) {
	  throw CanteraError("setPressure", "error");
	}
	setDensity(dd);
    }
 

    /// critical temperature 
    doublereal WaterTP::critTemperature() const { return m_sub->Tcrit(); }
        
    /// critical pressure
    doublereal WaterTP::critPressure() const { return m_sub->Pcrit(); }
        
    /// critical density
    doublereal WaterTP::critDensity() const { return m_sub->Rhocrit(); }
        
        

    void WaterTP::setTemperature(double temp) {
	State::setTemperature(temp);
	doublereal dd = density();
	m_sub->setState(temp, dd);
    }

     

    /// saturation pressure
    doublereal WaterTP::satPressure(doublereal t){
        doublereal pp = m_sub->psat(t);
	double dens = density();
	setTemperature(t);
	setDensity(dens);
	return pp;
    }
        


}
