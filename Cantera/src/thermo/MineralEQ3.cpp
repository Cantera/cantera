/**
 *  @file MineralEQ3.cpp
 *   Definition file for the MineralEQ3 class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and 
 * class \link Cantera::MineralEQ3 MineralEQ3\endlink)
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 *
 * Copyright 2001 California Institute of Technology
 */

/*
 * $Id$ 
 */

#include "ct_defs.h"
#include "mix_defs.h"
#include "MineralEQ3.h"
#include "SpeciesThermo.h"

#include "ThermoFactory.h"
#include "MineralEQ3.h"

#include <string>

using namespace std;

namespace Cantera {

  /*
   * ----  Constructors -------
   */

  /*
   * Default Constructor for the MineralEQ3 class
   */
  MineralEQ3::MineralEQ3():
    StoichSubstanceSSTP ()
  {
  }

  // Create and initialize a MineralEQ3 ThermoPhase object 
  // from an asci input file
  /*
   * @param infile name of the input file
   * @param id     name of the phase id in the file.
   *               If this is blank, the first phase in the file is used.
   */
  MineralEQ3::MineralEQ3(std::string infile, std::string id) :
    StoichSubstanceSSTP()
  {
    XML_Node* root = get_XML_File(infile);
    if (id == "-") id = "";
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
      throw CanteraError("MineralEQ3::MineralEQ3",
			  "Couldn't find phase name in file:" + id);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "MineralEQ3") {
      throw CanteraError("MineralEQ3::MineralEQ3",
			 "thermo model attribute must be StoichSubstance");
    }
    importPhase(*xphase, this);
  }

  // Full Constructor.
  /*
   *  @param phaseRef XML node pointing to a MineralEQ3 description
   *  @param id       Id of the phase. 
   */
  MineralEQ3::MineralEQ3(XML_Node& xmlphase, std::string id) :
    StoichSubstanceSSTP()
  {
    if (id != "") {
      std::string idxml = xmlphase["id"];
      if (id != idxml) {
	throw CanteraError("MineralEQ3::MineralEQ3",
			   "id's don't match");
      }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "MineralEQ3") {
      throw CanteraError("MineralEQ3::MineralEQ3",
			 "thermo model attribute must be StoichSubstance");
    }
    importPhase(xmlphase, this);
  }

  //! Copy constructor
  /*!
   * @param right Object to be copied
   */
  MineralEQ3::MineralEQ3(const MineralEQ3  &right) :
    StoichSubstanceSSTP()
  {
    *this = operator=(right);
  }
  
  //! Assignment operator
  /*!
   * @param right Object to be copied
   */
  MineralEQ3 & 
  MineralEQ3::operator=(const MineralEQ3 & right) {
    if (&right == this) {
      return *this;
    }
    StoichSubstanceSSTP::operator=(right);
    m_Mu0_tr_pr = right.m_Mu0_tr_pr;
    m_Entrop_tr_pr = right.m_Entrop_tr_pr;
    m_deltaG_formation_tr_pr = right.m_deltaG_formation_tr_pr;
    m_deltaH_formation_tr_pr = right.m_deltaH_formation_tr_pr;
    m_V0_tr_pr               = right.m_V0_tr_pr;
    m_a                      = right.m_a;
    m_b                      = right.m_b;
    m_c                      = right.m_c;

    return *this;
  }

  /*
   * Destructor for the routine (virtual)
   *        
   */
  MineralEQ3::~MineralEQ3() 
  {
  }

  // Duplication function
  /*
   * This virtual function is used to create a duplicate of the
   * current phase. It's used to duplicate the phase when given
   * a ThermoPhase pointer to the phase.
   *
   * @return It returns a ThermoPhase pointer.
   */
  ThermoPhase *MineralEQ3::duplMyselfAsThermoPhase() const {
    MineralEQ3 *stp = new MineralEQ3(*this);
    return (ThermoPhase *) stp;
  }


  /*
   * ---- Utilities -----
   */

  /*
   * Equation of state flag. Returns the value cStoichSubstance,
   * defined in mix_defs.h.
   */
  int MineralEQ3::eosType() const {
    return cStoichSubstance; 
  }

  /*
   * ---- Molar Thermodynamic properties of the solution ----
   */

  /**
   * ----- Mechanical Equation of State ------
   */

  /*
   * Pressure. Units: Pa.
   * For an incompressible substance, the density is independent
   * of pressure. This method simply returns the stored
   * pressure value.
   */ 
  doublereal MineralEQ3::pressure() const {
    return m_press;
  }
    
  /*
   * Set the pressure at constant temperature. Units: Pa.
   * For an incompressible substance, the density is 
   * independent of pressure. Therefore, this method only 
   * stores the specified pressure value. It does not 
   * modify the density.
   */
  void MineralEQ3::setPressure(doublereal p) {
    m_press = p;
  }

  /*
   * The isothermal compressibility. Units: 1/Pa.
   * The isothermal compressibility is defined as
   * \f[
   * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
   * \f]
   *
   *  It's equal to zero for this model, since the molar volume
   *  doesn't change with pressure or temperature.
   */
  doublereal MineralEQ3::isothermalCompressibility() const {
    return 0.0;
  } 
    
  /*
   * The thermal expansion coefficient. Units: 1/K.
   * The thermal expansion coefficient is defined as
   *
   * \f[
   * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
   * \f]
   *
   *  It's equal to zero for this model, since the molar volume
   *  doesn't change with pressure or temperature.
   */
  doublereal MineralEQ3::thermalExpansionCoeff() const {
    return 0.0;
  }
    
  /*
   * ---- Chemical Potentials and Activities ----
   */

  /*
   * This method returns the array of generalized
   * concentrations.  For a stoichiomeetric substance, there is
   * only one species, and the generalized concentration is 1.0.
   */
  void MineralEQ3::
  getActivityConcentrations(doublereal* c) const {
    c[0] = 1.0;
  }

  /*
   * The standard concentration. This is defined as the concentration 
   * by which the generalized concentration is normalized to produce 
   * the activity. 
   */ 
  doublereal MineralEQ3::standardConcentration(int k) const {
    return 1.0;
  }

  /*
   * Returns the natural logarithm of the standard 
   * concentration of the kth species
   */
  doublereal MineralEQ3::logStandardConc(int k) const {
    return 0.0;
  }

  /*
   * Returns the units of the standard and generalized
   * concentrations Note they have the same units, as their
   * ratio is defined to be equal to the activity of the kth
   * species in the solution, which is unitless.
   *
   * This routine is used in print out applications where the
   * units are needed. Usually, MKS units are assumed throughout
   * the program and in the XML input files.
   *
   *  uA[0] = kmol units - default  = 1
   *  uA[1] = m    units - default  = -nDim(), the number of spatial
   *                                dimensions in the Phase class.
   *  uA[2] = kg   units - default  = 0;
   *  uA[3] = Pa(pressure) units - default = 0;
   *  uA[4] = Temperature units - default = 0;
   *  uA[5] = time units - default = 0
   */
  void MineralEQ3::
  getUnitsStandardConc(doublereal *uA, int k, int sizeUA) const {
    for (int i = 0; i < 6; i++) {
      uA[i] = 0;
    }
  }

  /*
   *  ---- Partial Molar Properties of the Solution ----
   */



  /*
   * ---- Properties of the Standard State of the Species in the Solution
   * ----
   */

  /*
   * Get the array of chemical potentials at unit activity 
   * \f$ \mu^0_k \f$.
   *
   * For a stoichiometric substance, there is no activity term in 
   * the chemical potential expression, and therefore the
   * standard chemical potential and the chemical potential
   * are both equal to the molar Gibbs function.
   */
  void MineralEQ3::
  getStandardChemPotentials(doublereal* mu0) const {
    getGibbs_RT(mu0);
    mu0[0] *= GasConstant * temperature();
  }
    
  /*
   * Get the nondimensional Enthalpy functions for the species
   * at their standard states at the current 
   * <I>T</I> and <I>P</I> of the solution.
   * Molar enthalpy. Units: J/kmol.  For an incompressible,
   * stoichiometric substance, the internal energy is
   * independent of pressure, and therefore the molar enthalpy
   * is \f[ \hat h(T, P) = \hat u(T) + P \hat v \f], where the
   * molar specific volume is constant.
   */
  void MineralEQ3::getEnthalpy_RT(doublereal* hrt) const {
    getEnthalpy_RT_ref(hrt);
    doublereal RT = GasConstant * temperature();
    doublereal presCorrect = (m_press - m_p0) /  molarDensity();
    hrt[0] += presCorrect / RT;
  }

  /*
   * Get the array of nondimensional Entropy functions for the
   * standard state species
   * at the current <I>T</I> and <I>P</I> of the solution.
   */
  void MineralEQ3::getEntropy_R(doublereal* sr) const {
    getEntropy_R_ref(sr);
  }

  /*
   * Get the nondimensional Gibbs functions for the species
   * at their standard states of solution at the current T and P
   * of the solution
   */
  void MineralEQ3::getGibbs_RT(doublereal* grt) const {
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R[0];
  }

  /*
   * Get the nondimensional Gibbs functions for the standard
   * state of the species at the current T and P.
   */
  void MineralEQ3::getCp_R(doublereal* cpr) const {
    _updateThermo();
    cpr[0] = m_cp0_R[0];
  }
   
  /*
   * Molar internal energy (J/kmol).
   * For an incompressible,
   * stoichiometric substance, the molar internal energy is
   * independent of pressure. Since the thermodynamic properties
   * are specified by giving the standard-state enthalpy, the
   * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
   * enthalpy to compute the molar internal energy.
   */
  void MineralEQ3::getIntEnergy_RT(doublereal* urt) const {
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    doublereal PV = m_p0 / molarDensity();
    urt[0] = m_h0_RT[0] - PV / RT;
  }

  /*
   * ---- Thermodynamic Values for the Species Reference States ----
   */
  /*
   * Molar internal energy or the reference state at the current
   * temperature, T  (J/kmol).  
   * For an incompressible,
   * stoichiometric substance, the molar internal energy is
   * independent of pressure. Since the thermodynamic properties
   * are specified by giving the standard-state enthalpy, the
   * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
   * enthalpy to compute the molar internal energy.
   *
   * Note, this is equal to the standard state internal energy
   * evaluated at the reference pressure.
   */
  void MineralEQ3::getIntEnergy_RT_ref(doublereal* urt) const {
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    doublereal PV = m_p0 / molarDensity();
    urt[0] = m_h0_RT[0] - PV / RT;
  }
              
  /*
   * ---- Saturation Properties
   */
    
 

  /*
   * ---- Initialization and Internal functions
   */

  /**
   * @internal Initialize. This method is provided to allow
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
  void MineralEQ3::initThermo() {
  
    /*
     * Call the base class thermo initializer
     */
    StoichSubstanceSSTP::initThermo();
  }

  /**
   * setParameters:
   *
   *   Generic routine that is used to set the parameters used
   *   by this model.
   *        C[0] = density of phase [ kg/m3 ]
   */
  void MineralEQ3::setParameters(int n, doublereal * const c) {
    doublereal rho = c[0];
    setDensity(rho);
  }

  /**
   * getParameters:
   *
   *   Generic routine that is used to get the parameters used
   *   by this model.
   *        n = 1
   *        C[0] = density of phase [ kg/m3 ]
   */
  void MineralEQ3::getParameters(int &n, doublereal * const c) const {
    doublereal rho = density();
    n = 1;
    c[0] = rho;
  }

  // Initialize the phase parameters from an XML file.
  /*
   * initThermoXML()                 (virtual from ThermoPhase)
   *
   *  This gets called from importPhase(). It processes the XML file
   *  after the species are set up. This is the main routine for
   *  reading in activity coefficient parameters.
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
  void MineralEQ3::initThermoXML(XML_Node& phaseNode, std::string id) {
    /*
     * Find the Thermo XML node 
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("HMWSoln::initThermoXML",
			 "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");


    
  }

  /*
   * Reads an xml data block for the parameters needed by this
   * routine. eosdata is a reference to the xml thermo block, and looks
   * like this:
   * 
   *   <phase id="stoichsolid" >
   *     <thermo model="StoichSubstance">
   *         <density units="g/cm3">3.52</density>
   *     </thermo>
   *   </phase>
   */
  void MineralEQ3::setParametersFromXML(const XML_Node& eosdata) {
    std::string model = eosdata["model"];
    if (model != "StoichSubstance" && model != "MineralEQ3") {
      throw CanteraError("MineralEQ3::MineralEQ3",
			 "thermo model attribute must be StoichSubstance");
    }
    doublereal rho = getFloat(eosdata, "density", "toSI");
    setDensity(rho);
  }

 doublereal MineralEQ3::LookupGe(const std::string& elemName) {
#ifdef OLDWAY
    int num = sizeof(geDataTable) / sizeof(struct GeData);
    string s3 = elemName.substr(0,3);
    for (int i = 0; i < num; i++) {
      //if (!std::strncmp(elemName.c_str(), aWTable[i].name, 3)) {
      if (s3 == geDataTable[i].name) {
        return (geDataTable[i].GeValue);
      }
    }
    throw CanteraError("LookupGe", "element " + s + " not found");
    return -1.0;
#else
    int iE = elementIndex(elemName);
    if (iE < 0) {
      throw CanteraError("PDSS_HKFT::LookupGe", "element " + elemName + " not found");
    }
    doublereal geValue = entropyElement298(iE);
    if (geValue == ENTROPY298_UNKNOWN) {
      throw CanteraError("PDSS_HKFT::LookupGe",
                         "element " + elemName + " doesn not have a supplied entropy298");
    }
    geValue *= (-298.15);
    return geValue;
#endif
  }

 void MineralEQ3::convertDGFormation() {
    /*
     * Ok let's get the element compositions and conversion factors.
     */
    int ne = nElements();
    doublereal na;
    doublereal ge;
    string ename;

    doublereal totalSum = 0.0;
    for (int m = 0; m < ne; m++) {
      na = nAtoms(0, m);
      if (na > 0.0) {
	ename = elementName(m);
	ge = LookupGe(ename);
	totalSum += na * ge;
      }
    }
    // Add in the charge
    // if (m_charge_j != 0.0) {
    // ename = "H";
      // ge = LookupGe(ename);
    // totalSum -= m_charge_j * ge;
    //}
    // Ok, now do the calculation. Convert to joules kmol-1
    doublereal dg = m_deltaG_formation_tr_pr * 4.184 * 1.0E3;
    //! Store the result into an internal variable.
    m_Mu0_tr_pr = dg + totalSum;
  }

}


