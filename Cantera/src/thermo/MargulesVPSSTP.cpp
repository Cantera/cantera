/**
 *  @file MargulesVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations related to Margules
 *   expansions (see \ref thermoprops 
 *    and class \link Cantera::MargulesVPSSTP MargulesVPSSTP\endlink).
 *
 */
/*
 * Copywrite (2009) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 *  $Date$
 *  $Revision$
 */


#include "MargulesVPSSTP.h"
#include "ThermoFactory.h"


using namespace std;

namespace Cantera {

 static  const double xxSmall = 1.0E-150;
  /*
   * Default constructor.
   *
   */
  MargulesVPSSTP::MargulesVPSSTP() :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
  {
  }

  /*
   * Working constructors
   *
   *  The two constructors below are the normal way
   *  the phase initializes itself. They are shells that call
   *  the routine initThermo(), with a reference to the
   *  XML database to get the info for the phase.

   */
  MargulesVPSSTP::MargulesVPSSTP(std::string inputFile, std::string id) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
  {
    constructPhaseFile(inputFile, id);
  }

  MargulesVPSSTP::MargulesVPSSTP(XML_Node& phaseRoot, std::string id) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
  {
    constructPhaseXML(phaseRoot, id);
  }


  /*
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  MargulesVPSSTP::MargulesVPSSTP(const MargulesVPSSTP &b) :
    GibbsExcessVPSSTP()
  {
    MargulesVPSSTP::operator=(b);
  }

  /*
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  MargulesVPSSTP& MargulesVPSSTP::
  operator=(const MargulesVPSSTP &b) {
    if (&b == this) {
      return *this;
    }
   
    GibbsExcessVPSSTP::operator=(b);
    
    numBinaryInteractions_      = b.numBinaryInteractions_ ;
    m_HE_b_ij                   = b.m_HE_b_ij;
    m_HE_c_ij                   = b.m_HE_c_ij;
    m_HE_d_ij                   = b.m_HE_d_ij;
    m_SE_b_ij                   = b.m_SE_b_ij;
    m_SE_c_ij                   = b.m_SE_c_ij;
    m_SE_d_ij                   = b.m_SE_d_ij;
    m_pSpecies_A_ij             = b.m_pSpecies_A_ij;
    m_pSpecies_B_ij             = b.m_pSpecies_B_ij;
    formMargules_               = b.formMargules_;
    formTempModel_              = b.formTempModel_;
    
    return *this;
  }

  /**
   *
   * ~MargulesVPSSTP():   (virtual)
   *
   * Destructor: does nothing:
   *
   */
  MargulesVPSSTP::~MargulesVPSSTP() {
  }

  /*
   * This routine duplicates the current object and returns
   * a pointer to ThermoPhase.
   */
  ThermoPhase* 
  MargulesVPSSTP::duplMyselfAsThermoPhase() const {
    MargulesVPSSTP* mtp = new MargulesVPSSTP(*this);
    return (ThermoPhase *) mtp;
  }

  // Special constructor for a hard-coded problem
  /*
   *
   *   LiKCl treating the PseudoBinary layer as passthrough.
   *   -> test to predict the eutectic and liquidus correctly.
   *
   */
  MargulesVPSSTP::MargulesVPSSTP(int testProb)  :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
  {
  

    constructPhaseFile("LiKCl_liquid.xml", "");


    numBinaryInteractions_ = 1;

    m_HE_b_ij.resize(1);
    m_HE_c_ij.resize(1);
    m_HE_d_ij.resize(1);

    m_SE_b_ij.resize(1);
    m_SE_c_ij.resize(1);
    m_SE_d_ij.resize(1);
    
    m_pSpecies_A_ij.resize(1);
    m_pSpecies_B_ij.resize(1);


    m_HE_b_ij[0] = -17570E3;
    m_HE_c_ij[0] = -377.0E3;
    m_HE_d_ij[0] = 0.0;
      
    m_SE_b_ij[0] = -7.627E3;
    m_SE_c_ij[0] =  4.958E3;
    m_SE_d_ij[0] =  0.0;

    int iLiCl = speciesIndex("LiCl(L)");
    if (iLiCl < 0) {
      throw CanteraError("MargulesVPSSTP test1 constructor",
			 "Unable to find LiCl(L)");
    }
    m_pSpecies_B_ij[0] = iLiCl;


    int iKCl = speciesIndex("KCl(L)");
    if (iKCl < 0) {
      throw CanteraError("MargulesVPSSTP test1 constructor",
			 "Unable to find KCl(L)");
    }
    m_pSpecies_A_ij[0] = iKCl;
  }


  /*
   *  -------------- Utilities -------------------------------
   */

 
  // Equation of state type flag.
  /*
   * The ThermoPhase base class returns
   * zero. Subclasses should define this to return a unique
   * non-zero value. Known constants defined for this purpose are
   * listed in mix_defs.h. The MargulesVPSSTP class also returns
   * zero, as it is a non-complete class.
   */
  int MargulesVPSSTP::eosType() const { 
    return 0;
  }

  /*
   *   Import, construct, and initialize a phase
   *   specification from an XML tree into the current object.
   *
   * This routine is a precursor to constructPhaseXML(XML_Node*)
   * routine, which does most of the work.
   *
   * @param infile XML file containing the description of the
   *        phase
   *
   * @param id  Optional parameter identifying the name of the
   *            phase. If none is given, the first XML
   *            phase element will be used.
   */
  void MargulesVPSSTP::constructPhaseFile(std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("MargulesVPSSTP:constructPhaseFile",
                         "input file is null");
    }
    string path = findInputFile(inputFile);
    std::ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("MargulesVPSSTP:constructPhaseFile","could not open "
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
      throw CanteraError("MargulesVPSSTP:constructPhaseFile",
                         "ERROR: Can not find phase named " +
                         id + " in file named " + inputFile);
    }
    fxml_phase->copy(&phaseNode_XML);
    constructPhaseXML(*fxml_phase, id);
    delete fxml;
  }

  /*
   *   Import, construct, and initialize a HMWSoln phase 
   *   specification from an XML tree into the current object.
   *
   *   Most of the work is carried out by the cantera base
   *   routine, importPhase(). That routine imports all of the
   *   species and element data, including the standard states
   *   of the species.
   *
   *   Then, In this routine, we read the information 
   *   particular to the specification of the activity 
   *   coefficient model for the Pitzer parameterization.
   *
   *   We also read information about the molar volumes of the
   *   standard states if present in the XML file.
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
  void MargulesVPSSTP::constructPhaseXML(XML_Node& phaseNode, std::string id) {
    string stemp;
    if (id.size() > 0) {
      string idp = phaseNode.id();
      if (idp != id) {
	throw CanteraError("MargulesVPSSTP::constructPhaseXML", 
			   "phasenode and Id are incompatible");
      }
    }

    /*
     * Find the Thermo XML node 
     */
    if (!phaseNode.hasChild("thermo")) {
      throw CanteraError("MargulesVPSSTP::constructPhaseXML",
			 "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Make sure that the thermo model is Margules
     */ 
    stemp = thermoNode.attrib("model");
    string formString = lowercase(stemp);
    if (formString != "margules") {
      throw CanteraError("MargulesVPSSTP::constructPhaseXML",
			 "model name isn't Margules: " + formString);
    
    }

    /*
     * Call the Cantera importPhase() function. This will import
     * all of the species into the phase. This will also handle
     * all of the solvent and solute standard states
     */
    bool m_ok = importPhase(phaseNode, this);
    if (!m_ok) {
      throw CanteraError("MargulesVPSSTP::constructPhaseXML","importPhase failed "); 
    }
    
  }


  /*
   * ------------ Molar Thermodynamic Properties ----------------------
   */


  /*
   * - Activities, Standard States, Activity Concentrations -----------
   */


  doublereal MargulesVPSSTP::standardConcentration(int k) const {
    err("standardConcentration");
    return -1.0;
  }

  doublereal MargulesVPSSTP::logStandardConc(int k) const {
    err("logStandardConc");
    return -1.0;
  }

  // Get the array of non-dimensional molar-based activity coefficients at
  // the current solution temperature, pressure, and solution concentration.
  /*
   * @param ac Output vector of activity coefficients. Length: m_kk.
   */
  void MargulesVPSSTP::getActivityCoefficients(doublereal* ac) const {
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();

    /*
     * take the exp of the internally storred coefficients.
     */
    for (int k = 0; k < m_kk; k++) {
      ac[k] = exp(lnActCoeff_Scaled_[k]);      
    }
  }

  /*
   * ------------ Partial Molar Properties of the Solution ------------
   */



  void MargulesVPSSTP::getElectrochemPotentials(doublereal* mu) const {
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (int k = 0; k < m_kk; k++) {
      mu[k] += ve*charge(k);
    }
  }


  void MargulesVPSSTP::getChemPotentials(doublereal* mu) const {
    doublereal xx;
    /*
     * First get the standard chemical potentials in
     * molar form.
     *  -> this requires updates of standard state as a function
     *     of T and P
     */
    getStandardChemPotentials(mu);
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();
    /*
     *
     */
    doublereal RT = GasConstant * temperature();
    for (int k = 0; k < m_kk; k++) {
      xx = fmaxx(moleFractions_[k], xxSmall);
      mu[k] += RT * (log(xx) + lnActCoeff_Scaled_[k]);      
    }
  }

 
  // Returns an array of partial molar enthalpies for the species
  // in the mixture.
  /*
   * Units (J/kmol)
   *
   * For this phase, the partial molar enthalpies are equal to the
   * standard state enthalpies modified by the derivative of the
   * molality-based activity coefficent wrt temperature
   *
   *  \f[
   * \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
   * \f]
   *
   */
  void MargulesVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const {
    /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * dimensionalize it.
     */
    double T = temperature();
    double RT = GasConstant * T;
    for (int k = 0; k < m_kk; k++) {
      hbar[k] *= RT;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally storred molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    double RTT = RT * T;
    for (int k = 0; k < m_kk; k++) {
      hbar[k] -= RTT * dlnActCoeffdT_Scaled_[k];
    }
  }

  // Returns an array of partial molar entropies for the species
  // in the mixture.
  /*
   * Units (J/kmol)
   *
   * For this phase, the partial molar enthalpies are equal to the
   * standard state enthalpies modified by the derivative of the
   * activity coefficent wrt temperature
   *
   *  \f[
   * \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
   * \f]
   *
   */
  void MargulesVPSSTP::getPartialMolarEntropies(doublereal* sbar) const {
    double xx;
    /*
     * Get the nondimensional standard state entropies
     */
    getEntropy_R(sbar);
    double T = temperature();
    /*
     * Update the activity coefficients, This also update the
     * internally storred molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (int k = 0; k < m_kk; k++) {
      xx = fmaxx(moleFractions_[k], xxSmall);
      sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }  
    /*
     * dimensionalize it.
     */
   for (int k = 0; k < m_kk; k++) {
      sbar[k] *= GasConstant;
    }
  }
  
  

  doublereal MargulesVPSSTP::err(std::string msg) const {
    throw CanteraError("MargulesVPSSTP","Base class method "
		       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
  }


  /*
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
  void MargulesVPSSTP::initThermo() {
    initLengths();
    GibbsExcessVPSSTP::initThermo();
  }


  //   Initialize lengths of local variables after all species have
  //   been identified.
  void  MargulesVPSSTP::initLengths() {
    m_kk = nSpecies();
 
  }

  /*
   * initThermoXML()                (virtual from ThermoPhase)
   *   Import and initialize a ThermoPhase object
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
  void MargulesVPSSTP::initThermoXML(XML_Node& phaseNode, std::string id) {
    string subname = "MargulesVPSSTP::initThermoXML";
    string stemp;

    /*
     * Check on the thermo field. Must have:
     * <thermo model="IdealSolidSolution" />
     */
  
    XML_Node& thermoNode = phaseNode.child("thermo");
    string mStringa = thermoNode.attrib("model");
    string mString = lowercase(mStringa);
    if (mString != "margules") {
      throw CanteraError(subname.c_str(),
			 "Unknown thermo model: " + mStringa);
    }
 

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node *acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
      XML_Node& acNode = thermoNode.child("activityCoefficients");
      acNodePtr = &acNode;
      string mStringa = thermoNode.attrib("model");
      string mString = lowercase(mStringa);
      if (mString != "margules") {
	throw CanteraError(subname.c_str(),
			   "Unknown activity coefficient model: " + mStringa);
      }
      int n = acNodePtr->nChildren();
      for (int i = 0; i < n; i++) {
	XML_Node &xmlACChild = acNodePtr->child(i);
	stemp = xmlACChild.name();
	string nodeName = lowercase(stemp);
	/*
	 * Process a binary salt field, or any of the other XML fields
	 * that make up the Pitzer Database. Entries will be ignored
	 * if any of the species in the entry isn't in the solution.
	 */
	if (nodeName == "binaryneutralspeciesparameters") {
	  readXMLBinarySpecies(xmlACChild);

	}
      }
    }

    /*
     * Go down the chain
     */
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);


  }

  // Update the activity coefficients
  /*
   * This function will be called to update the internally storred
   * natural logarithm of the activity coefficients
   *
   *   he = X_A X_B(B + C(X_A - X_B))
   */
  void MargulesVPSSTP::s_update_lnActCoeff() const {
    int iA, iB;
    double XA, XB, g0 , g1;
    double T = temperature();

    fvo_zero_dbl_1(lnActCoeff_Scaled_, m_kk);

    double RT = GasConstant * temperature();
    for (int i = 0; i <  numBinaryInteractions_; i++) {
      iA =  m_pSpecies_A_ij[i];    
      iB =  m_pSpecies_B_ij[i];

      XA = moleFractions_[iA];
      XB = moleFractions_[iB];
      
      g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT ;
      g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;
      
      lnActCoeff_Scaled_[iA] += XB * XB * (g0 + g1 * (XB - XA));
      lnActCoeff_Scaled_[iB] += XA * XA * g0 +  XA * XB * g1 * (2 * XA);
    }
  }

  // Update the derivative of the log of the activity coefficients wrt T
  /*
   * This function will be called to update the internally storred
   * natural logarithm of the activity coefficients
   *
   *   he = X_A X_B(B + C(X_A - X_B))
   */
  void MargulesVPSSTP::s_update_dlnActCoeff_dT() const {
    int iA, iB;
    doublereal XA, XB, h0 , h1;
    doublereal T = temperature();

    fvo_zero_dbl_1(dlnActCoeffdT_Scaled_, m_kk);

    doublereal RTT = GasConstant * T * T;
    for (int i = 0; i <  numBinaryInteractions_; i++) {
      iA =  m_pSpecies_A_ij[i];    
      iB =  m_pSpecies_B_ij[i];

      XA = moleFractions_[iA];
      XB = moleFractions_[iB];
      
      h0 = m_HE_b_ij[i];
      h1 = m_HE_c_ij[i];
      
      dlnActCoeffdT_Scaled_[iA] += -(XB * XB * (h0 + h1 * (XB - XA))) / RTT;
      dlnActCoeffdT_Scaled_[iB] += -(XA * XA * h0 +  XA * XB * h1 * (2 * XA))/RTT;
    }
  }

  void MargulesVPSSTP::getdlnActCoeffdT(doublereal *dlnActCoeffdT) const {
    s_update_dlnActCoeff_dT();
    for (int k = 0; k < m_kk; k++) {
      dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
  }

  // Update the derivative of the log of the activity coefficients wrt ln(X)
  /*
   * This function will be called to update the internally stored gradients of the 
   * logarithm of the activity coefficients.  These are used in the determination 
   * of the diffusion coefficients.
   *
   *   he = X_A X_B(B + C(X_A - X_B))
   */
  void MargulesVPSSTP::s_update_dlnActCoeff_dlnC() const {
    int iA, iB;
    doublereal XA, XB, g0 , g1;
    doublereal T = temperature();

    fvo_zero_dbl_1(dlnActCoeffdlnC_Scaled_, m_kk);

    doublereal RT = GasConstant * T;
    for (int i = 0; i <  numBinaryInteractions_; i++) {
      iA =  m_pSpecies_A_ij[i];    
      iB =  m_pSpecies_B_ij[i];

      XA = moleFractions_[iA];
      XB = moleFractions_[iB];
      
      g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT ;
      g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;
      
      dlnActCoeffdlnC_Scaled_[iA] += XA * ( ( - 2.0 + 2.0 * XA ) * g0
					    + ( - 4.0 + 10.0 * XA - 6.0 * XA*XA ) * g1 ) ;
      dlnActCoeffdlnC_Scaled_[iB] += XB * ( ( - 2.0 + 2.0 * XB ) * g0
					    + (   2.0 -  8.0 * XB + 6.0 * XB*XB ) * g1 ) ;
    }
  }

  void MargulesVPSSTP::getdlnActCoeffdlnC(doublereal *dlnActCoeffdlnC) const {
    s_update_dlnActCoeff_dlnC();
    for (int k = 0; k < m_kk; k++) {
      dlnActCoeffdlnC[k] = dlnActCoeffdlnC_Scaled_[k];
    }
  }


  void MargulesVPSSTP::resizeNumInteractions(const int num) {
    numBinaryInteractions_ = num;
    m_HE_b_ij.resize(num, 0.0);
    m_HE_c_ij.resize(num, 0.0);
    m_HE_d_ij.resize(num, 0.0);
    m_SE_b_ij.resize(num, 0.0);
    m_SE_c_ij.resize(num, 0.0);
    m_SE_d_ij.resize(num, 0.0);

    m_pSpecies_A_ij.resize(num, -1);
    m_pSpecies_B_ij.resize(num, -1);

  }


  /*
   * Process an XML node called "binaryNeutralSpeciesParameters"
   * This node contains all of the parameters necessary to describe
   * the Margules Interaction for a single binary interaction
   * This function reads the XML file and writes the coefficients
   * it finds to an internal data structures.
   */
  void MargulesVPSSTP::readXMLBinarySpecies(XML_Node &xmLBinarySpecies) {
    string xname = xmLBinarySpecies.name();
    if (xname != "binaryNeutralSpeciesParameters") {
      throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies",
                         "Incorrect name for processing this routine: " + xname);
    }
    double *charge = DATA_PTR(m_speciesCharge);
    string stemp;
    int nParamsFound;
    vector_fp vParams;
    string iName = xmLBinarySpecies.attrib("speciesA");
    if (iName == "") {
      throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "no speciesA attrib");
    }
    string jName = xmLBinarySpecies.attrib("speciesB");
    if (jName == "") {
      throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "no speciesB attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    int iSpecies = speciesIndex(iName);
    if (iSpecies < 0) {
      return;
    }
    string ispName = speciesName(iSpecies);
    if (charge[iSpecies] != 0) {
      throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "speciesA charge problem");
    }
    int jSpecies = speciesIndex(jName);
    if (jSpecies < 0) {
      return;
    }
    string jspName = speciesName(jSpecies);
    if (charge[jSpecies] != 0) {
      throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "speciesB charge problem");
    }

    resizeNumInteractions(numBinaryInteractions_ + 1);
    int iSpot = numBinaryInteractions_ - 1;
    m_pSpecies_A_ij[iSpot] = iSpecies;
    m_pSpecies_B_ij[iSpot] = jSpecies;
  
    int num = xmLBinarySpecies.nChildren();
    for (int iChild = 0; iChild < num; iChild++) {
      XML_Node &xmlChild = xmLBinarySpecies.child(iChild);
      stemp = xmlChild.name();
      string nodeName = lowercase(stemp);
      /*
       * Process the binary species interaction child elements
       */
      if (nodeName == "excessenthalpy") {
        /*
         * Get the string containing all of the values
         */
        getFloatArray(xmlChild, vParams, true, "", "excessEnthalpy");
        nParamsFound = vParams.size();
       
	if (nParamsFound != 2) {
	  throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies::excessEnthalpy for " + ispName
			     + "::" + jspName,
			     "wrong number of params found");
	}
	m_HE_b_ij[iSpot] = vParams[0];
	m_HE_c_ij[iSpot] = vParams[1];
      }

      if (nodeName == "excessentropy") {
        /*
         * Get the string containing all of the values
         */
        getFloatArray(xmlChild, vParams, true, "", "excessEntropy");
        nParamsFound = vParams.size();
       
	if (nParamsFound != 2) {
	  throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies::excessEntropy for " + ispName
			     + "::" + jspName,
			     "wrong number of params found");
	}
	m_SE_b_ij[iSpot] = vParams[0];
	m_SE_c_ij[iSpot] = vParams[1];
      }


    } 
    
  }

  /**
   * Format a summary of the mixture state for output.
   */           
  std::string MargulesVPSSTP::report(bool show_thermo) const {


    char p[800];
    string s = "";
    try {
      if (name() != "") {
	sprintf(p, " \n  %s:\n", name().c_str());
	s += p;
      }
      sprintf(p, " \n       temperature    %12.6g  K\n", temperature());
      s += p;
      sprintf(p, "          pressure    %12.6g  Pa\n", pressure());
      s += p;
      sprintf(p, "           density    %12.6g  kg/m^3\n", density());
      s += p;
      sprintf(p, "  mean mol. weight    %12.6g  amu\n", meanMolecularWeight());
      s += p;

      doublereal phi = electricPotential();
      sprintf(p, "         potential    %12.6g  V\n", phi);
      s += p;

      int kk = nSpecies();
      array_fp x(kk);
      array_fp molal(kk);
      array_fp mu(kk);
      array_fp muss(kk);
      array_fp acMolal(kk);
      array_fp actMolal(kk);
      getMoleFractions(&x[0]);
   
      getChemPotentials(&mu[0]);
      getStandardChemPotentials(&muss[0]);
      getActivities(&actMolal[0]);
 

      if (show_thermo) {
        sprintf(p, " \n");
        s += p;
        sprintf(p, "                          1 kg            1 kmol\n");
        s += p;
        sprintf(p, "                       -----------      ------------\n");
        s += p;
        sprintf(p, "          enthalpy    %12.6g     %12.4g     J\n", 
		enthalpy_mass(), enthalpy_mole());
        s += p;
        sprintf(p, "   internal energy    %12.6g     %12.4g     J\n", 
		intEnergy_mass(), intEnergy_mole());
        s += p;
        sprintf(p, "           entropy    %12.6g     %12.4g     J/K\n", 
		entropy_mass(), entropy_mole());
        s += p;
        sprintf(p, "    Gibbs function    %12.6g     %12.4g     J\n", 
		gibbs_mass(), gibbs_mole());
        s += p;
        sprintf(p, " heat capacity c_p    %12.6g     %12.4g     J/K\n", 
		cp_mass(), cp_mole());
        s += p;
        try {
	  sprintf(p, " heat capacity c_v    %12.6g     %12.4g     J/K\n", 
		  cv_mass(), cv_mole());
	  s += p;
        }
        catch(CanteraError) {
	  sprintf(p, " heat capacity c_v    <not implemented>       \n");
	  s += p;
        }
      }
  
    } catch (CanteraError) {
      ;
    }
    return s;
  }

 
}

