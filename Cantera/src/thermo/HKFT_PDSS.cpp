/*
 * $Id$
 */
#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "HKFT_PDSS.h"
#include "WaterPDSS.h"
#include "WaterProps.h"

#include "ThermoPhase.h"

using namespace std;

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */





  HKFT_PDSS::HKFT_PDSS(ThermoPhase *tp, int spindex) :
    PDSS(tp, spindex)
  {
  }


  HKFT_PDSS::HKFT_PDSS(ThermoPhase *tp, int spindex, std::string inputFile, std::string id) :
    PDSS(tp, spindex, inputFile, id)
  {
  }

  HKFT_PDSS::HKFT_PDSS(ThermoPhase *tp, int spindex, XML_Node& phaseRoot, std::string id) :
    PDSS(tp, spindex, phaseRoot, id)
  {
  }


  HKFT_PDSS::HKFT_PDSS(const HKFT_PDSS &b) :
    PDSS(b)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * Assignment operator
   */
  HKFT_PDSS& HKFT_PDSS::operator=(const HKFT_PDSS&b) {
    if (&b == this) return *this;
    m_tp = b.m_tp;
    m_spindex = b.m_spindex;
    m_temp = b.m_temp;
    m_dens = b.m_dens;
    m_mw = b.m_mw;
    return *this;
  }
  
  /**
   * Destructor for the HKFT_PDSS class
   */
  HKFT_PDSS::~HKFT_PDSS() { 
  }
  
  void HKFT_PDSS::constructHKFT_PDSS(ThermoPhase *tp, int spindex) {
    initThermo();
  }
  
   
  /**
   * constructHKFT_PDSSXML:
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
  void HKFT_PDSS::constructHKFT_PDSSXML(ThermoPhase *tp, int spindex, 
					XML_Node& phaseNode, std::string id) {
    initThermo();
  }

   
  /**
   * constructHKFT_PDSSFile():
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
  void HKFT_PDSS::constructHKFT_PDSSFile(ThermoPhase *tp, int spindex,
					 std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("WaterTp::initThermo",
			 "input file is null");
    }
    string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("HKFT_PDSS::initThermo","could not open "
			 +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */

    XML_Node *fxml = new XML_Node();
    fxml->build(fin);
    XML_Node *fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
      throw CanteraError("HKFT_PDSS::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }	
    constructHKFT_PDSSXML(tp, spindex, *fxml_phase, id);
    delete fxml;
  }

  void HKFT_PDSS::
  initThermoXML(XML_Node& phaseNode, std::string id) {
    initThermo();
  }

  void HKFT_PDSS::initThermo() {
    
    /*
     *  Section to initialize  m_Z_pr_tr and   m_Y_pr_tr
     */
    double temp = 273.15 + 25.;
    double pres = OneAtm;
    double relepsilon = m_wprops->relEpsilon(temp, pres, 0);
 
    m_Z_pr_tr = -1.0 / relepsilon;
    //double m_Z_pr_tr = -0.0127803;
    //printf("m_Z_pr_tr = %20.10g\n",  m_Z_pr_tr ); 
    double drelepsilondT = m_wprops->relEpsilon(temp, pres, 1);
    //double m_Y_pr_tr = -5.799E-5;
    m_Y_pr_tr = drelepsilondT / (relepsilon * relepsilon);
    //printf("m_Y_pr_tr = %20.10g\n",  m_Y_pr_tr ); 
  }

  void HKFT_PDSS::
  setParametersFromXML(const XML_Node& eosdata) {
  }

  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal HKFT_PDSS::
  enthalpy_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal HKFT_PDSS::
  intEnergy_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  entropy_mole() const {

    throw CanteraError("HKFT_PDSS::entropy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1 K-1.
   */
  doublereal HKFT_PDSS::gibbs_mole() const {
    double val = deltaG();
    return (m_Mu0_tr_pr + val);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cp_mole() const {
    throw CanteraError("HKFT_PDSS::cp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cv_mole() const {
    throw CanteraError("HKFT_PDSS::cv_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Return the difference in enthalpy between current p
   * and ref p0, in mks units of
   * in units of J kmol-1
   */
  doublereal HKFT_PDSS::
  enthalpyDelp_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpy_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate difference in the internal energy between current p
   * and ref p0, in mks units of
   * J kmol-1 
   */
  doublereal HKFT_PDSS::
  intEnergyDelp_mole() const {
    throw CanteraError("HKFT_PDSS::enthalpyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   *  Return the difference in entropy between current p
   * and ref p0, in mks units of
   * J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  entropyDelp_mole() const {

    throw CanteraError("HKFT_PDSS::entropyDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in Gibbs free energy between current p and
   * the ref p0, in mks units of
   * J kmol-1 K-1.
   */
  doublereal HKFT_PDSS::
  gibbsDelp_mole() const {
    throw CanteraError("HKFT_PDSS::gibbsDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in the constant pressure heat capacity
   * between the current p and the ref p0,
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cpDelp_mole() const {
    throw CanteraError("HKFT_PDSS::cpDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the difference in constant volume heat capacity
   * between the current p and the ref p0
   * in mks units of J kmol-1 K-1
   */
  doublereal HKFT_PDSS::
  cvDelp_mole() const {
    throw CanteraError("HKFT_PDSS::cvDelp_mole()", "unimplemented");
    return (0.0);
  }

  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal HKFT_PDSS::
  pressure() const {
    throw CanteraError("HKFT_PDSS::pressure()", "unimplemented");
    return (0.0);
  }
        
  void HKFT_PDSS::
  setPressure(doublereal p) {
    throw CanteraError("HKFT_PDSS::pressure()", "unimplemented");
  }
 

  /// critical temperature 
  doublereal HKFT_PDSS::critTemperature() const { 
    throw CanteraError("HKFT_PDSS::critTemperature()", "unimplemented");
    return (0.0);
  }
        
  /// critical pressure
  doublereal HKFT_PDSS::critPressure() const {
    throw CanteraError("HKFT_PDSS::critPressure()", "unimplemented");
    return (0.0);
  }
        
  /// critical density
  doublereal HKFT_PDSS::critDensity() const {
    throw CanteraError("HKFT_PDSS::critDensity()", "unimplemented");
    return (0.0);
  }
        
  void HKFT_PDSS::setDensity(double dens) {
    m_dens = dens;
  }

  double HKFT_PDSS::density() const {
    return m_dens;
  }

  double HKFT_PDSS::temperature() const {
    return m_temp;
  }
 
  void HKFT_PDSS::setTemperature(double temp) {
    m_temp = temp;
  }

  doublereal HKFT_PDSS::molecularWeight() const {
    return m_mw;
  }
  void HKFT_PDSS::setMolecularWeight(double mw) {
    m_mw = mw;
  }

  void HKFT_PDSS::setState_TP(double temp, double pres) {
    throw CanteraError("HKFT_PDSS::setState_TP()", "unimplemented");
  }

  /// saturation pressure
  doublereal HKFT_PDSS::satPressure(doublereal t){
    throw CanteraError("HKFT_PDSS::satPressure()", "unimplemented");
    return (0.0);
  }


  double HKFT_PDSS::deltaG() const {
    
    double pbar = m_pres * 1.0E-5;
    double m_presR_bar = OneAtm * 1.0E-5;

    double sterm = -  m_Entrop_tr_pr * (m_temp - 298.15);

    double c1term = -m_c1 * (m_temp * log(m_temp/298.15) - (m_temp - 298.15));
    double a1term = m_a1 * (pbar - m_presR_bar);

    double a2term = m_a2 * log((2600. + pbar)/(2600. + m_presR_bar));

    double c2term = -m_c2 * (( 1.0/(m_temp - 228.) - 1.0/(298.15 - 228.) ) * (228. - m_temp)/228.
			     - m_temp / (228.*228.) * log( (298.15*(m_temp-228.)) / (m_temp*(298.15-228.)) ));
    
    double a3term = m_a3 / (m_temp - 228.) * (pbar - m_presR_bar);

    double a4term = m_a4 / (m_temp - 228.) * log((2600. + pbar)/(2600. + m_presR_bar));

    double nu = 166027;
    double r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
  
    double gval = gstar(m_temp, m_pres, 0);

    double r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;

    double omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval)  );

    double relepsilon = m_wprops->relEpsilon(m_temp, m_pres, 0);

    double Z = -1.0 / relepsilon;

    double wterm = - omega_j * (Z + 1.0);

    double wrterm = m_omega_pr_tr * (m_Z_pr_tr + 1.0);

    double yterm = m_omega_pr_tr * m_Y_pr_tr * (m_temp - 298.15);

    double deltaG_calgmol = sterm + c1term + a1term + a2term + c2term + a3term + a4term + wterm + wrterm + yterm;

    // Convert to Joules / kmol
    double deltaG = deltaG_calgmol * 1.0E3 * 4.184;
    return deltaG;
  }
   
  double HKFT_PDSS::electrostatic_radii_calc() {
    return 0.0;
  }

  //! Internal formula for the calculation of a_g()
  /*
   * The output of this is in units of Angstroms
   */
  double HKFT_PDSS::ag(const double temp, const int ifunc) const {
    static double ag_coeff[3] = { -2.037662,  5.747000E-3,  -6.557892E-6};
    if (ifunc == 0) {
      double t2 = temp * temp;
      double val = ag_coeff[0] + ag_coeff[1] * temp + ag_coeff[2] * t2;
      return val;
    } else if (ifunc == 1) {
      return  ag_coeff[1] + ag_coeff[2] * 2.0 * temp;
    }
    if (ifunc != 2) {
      return 0.0;
    }
    return ag_coeff[2] * 2.0;;
  }

  //! Internal formula for the calculation of b_g()
  /*
   * the output of this is unitless
   */
  double HKFT_PDSS::bg(const double temp, const int ifunc) const {
    static double bg_coeff[3] = { 6.107361, -1.074377E-2,  1.268348E-5};
    if (ifunc == 0) {
      double t2 = temp * temp;
      double val = bg_coeff[0] + bg_coeff[1] * temp + bg_coeff[2] * t2;
      return val;
    }   else if (ifunc == 1) {
      return bg_coeff[1] + bg_coeff[2] * 2.0 * temp;
    }
    if (ifunc != 2) {
      return 0.0;
    }
    return bg_coeff[2] * 2.0;
  }

  double HKFT_PDSS::f(const double temp, const double pres, const int ifunc) const {
    
    static double af_coeff[3] = { 3.666666E1, -0.1504956E-9, 0.5107997E-13};
    double TC = temp - 273.15;
    double presBar = pres / 1.0E5;

    if (TC < 155.0) return 0.0;
    if (TC > 355.0) TC = 355.0;
    if (presBar > 1000.) return 0.0;
    

    double T1 = (TC-155.0)/300.;
    double fac1;

    double p2 = presBar * presBar;
    double p3 = presBar * p2;
    double p4 = p2 * p2;
    double fac2 = af_coeff[1] * p3 +   af_coeff[2] * p4;
    if (ifunc == 0) {
      fac1 = pow(T1,4.8) + af_coeff[0] * pow(T1, 16.0);
      return fac1 * fac2;
    } else if (ifunc == 1) {
      fac1 =  (4.8 * pow(T1,3.8) + 16.0 * af_coeff[0] * pow(T1, 15.0)) / 300.;
      return fac1 * fac2;
    } else if (ifunc == 2) {
      fac1 =  (4.8 * 3.8 * pow(T1,2.8) + 16.0 * 15.0 * af_coeff[0] * pow(T1, 14.0)) / (300. * 300.);
      return fac1 * fac2;
    } else if (ifunc == 3) {
      fac1 = pow(T1,4.8) + af_coeff[0] * pow(T1, 16.0);
      fac2 = (3.0 * af_coeff[1] * p2 +   4.0 * af_coeff[2] * p3 )/ 1.0E5;
      return fac1 * fac2;
    } else {
      throw CanteraError("HKFT_PDSS::gg", "unimplemented");
    }
    return 0.0;
  }

  double HKFT_PDSS::g(const double temp, const double pres, const int ifunc) const {
    double afunc = ag(temp, 0);
    double bfunc = bg(temp, 0);
    m_waterSS->setState_TP(temp, pres);
    m_densWaterSS = m_waterSS->density();
    // density in gm cm-3
    double dens = m_densWaterSS * 1.0E-3;
    double gval = afunc * pow((1.0-dens), bfunc);
    if (dens >= 1.0) {
      return 0.0;
    }
    if (ifunc == 0) {  
      return gval;

    } else if (ifunc == 1 || ifunc == 2) {
      double afuncdT = ag(temp, 1);
      double bfuncdT = bg(temp, 1);
      double alpha   = m_waterSS->thermalExpansionCoeff();

      double fac1 = afuncdT * gval / afunc;
      double fac2 = bfuncdT * gval * log(1.0 - dens);
      double fac3 = gval * alpha * bfunc * dens / (1.0 - dens);

      double dgdt = fac1 + fac2 + fac3;
      if (ifunc == 1) {
	return dgdt;
      }

      double afuncdT2 = ag(temp, 2);
      double bfuncdT2 = bg(temp, 2);

      double dfac1dT = dgdt * afuncdT / afunc + afuncdT2 * gval / afunc 
	-  afuncdT * afuncdT * gval / (afunc * afunc);

      double ddensdT = - alpha * dens;
      double dfac2dT =  bfuncdT2 * gval * log(1.0 - dens) 
	+ bfuncdT * dgdt * log(1.0 - dens) 
	- bfuncdT * gval /(1.0 - dens) * ddensdT;

      double dalphadT = m_waterSS->dthermalExpansionCoeffdT();
      
      double dfac3dT = dgdt * alpha * bfunc * dens / (1.0 - dens)
	+ gval * dalphadT * bfunc * dens / (1.0 - dens)
	+ gval * alpha * bfuncdT * dens / (1.0 - dens)
	+ gval * alpha * bfunc * ddensdT / (1.0 - dens)
	- gval * alpha * bfunc * dens / ((1.0 - dens) * (1.0 - dens)) * ddensdT;

      return dfac1dT + dfac2dT + dfac3dT;

    } else {
      throw CanteraError("HKFT_PDSS::gg", "unimplemented");
    }
    return 0.0;
  }

  double HKFT_PDSS::gstar(const double temp, const double pres, const int ifunc) const {
    double gval = g(temp, pres, ifunc);
    double fval = f(temp, pres, ifunc);
    return gval - fval;
  }

  /* awData structure */
  /**
   * Database for atomic molecular weights
   *
   *  Values are taken from the 1989 Standard Atomic Weights, CRC
   *
   *  awTable[] is a static function with scope limited to this file.
   *  It can only be referenced via the static Elements class function,
   *  LookupWtElements().
   *
   *  units = kg / kg-mol (or equivalently gm / gm-mol)
   *
   * (note: this structure was picked because it's simple, compact,
   *          and extensible).
   *
   */
  struct GeData {
    char name[4];     ///< Null Terminated name, First letter capitalized
    double GeValue;   /// < Gibbs free energies of elements J kmol-1
  };
  

  //! Values of G_elements(T=298.15,1atm) 
  /*!
   *  all units are Joules kmol-1
   */
  static struct GeData geDataTable[] = {
    {"H",   -19.48112E6}, // NIST Webbook - Cox, Wagman 1984
    {"Na",  -15.29509E6}, // NIST Webbook - Cox, Wagman 1984
    {"O",   -30.58303E6}, // NIST Webbook - Cox, Wagman 1984
    {"Cl",  -33.25580E6}, // NIST Webbook - Cox, Wagman 1984
    {"Si",   -5.61118E6}, // Janaf
    {"C",    -1.71138E6}, // barin, Knack, NBS Bulletin 1971
    {"S",    -9.55690E6}, // Yellow - webbook
    {"Al",   -8.42870E6}, // Webbook polynomial
    {"K",   -19.26943E6}  // Webbook
  };

  //!  Static function to look up Element Free Energies
  /*!
   *
   *   This static function looks up the argument string in the
   *   database above and returns the associated Gibbs Free energies.
   
   *
   *  @param  ElemName  String. Only the first 3 characters are significant
   *
   *  @return
   *    Return value contains the Gibbs free energy for that element
   *
   *  @exception CanteraError
   *    If a match is not found, a CanteraError is thrown as well
   */
  double HKFT_PDSS::LookupGe(const std::string& s) {
    int num = sizeof(geDataTable) / sizeof(struct GeData);
    string s3 = s.substr(0,3);
    for (int i = 0; i < num; i++) {
      //if (!std::strncmp(s.c_str(), aWTable[i].name, 3)) {
      if (s3 == geDataTable[i].name) {
        return (geDataTable[i].GeValue);
      }
    }
    throw CanteraError("LookupGe", "element not found");
    return -1.0;
  }

  void HKFT_PDSS::convertDGFormation() {
    /*
     * Ok let's get the element compositions and conversion factors.
     */
    int ne = m_tp->nElements();
    double na;
    double ge;
    string ename;

    double totalSum = 0.0;
    for (int m = 0; m < ne; m++) {
      na = m_tp->nAtoms(m_spindex, m);
      if (na > 0.0) {
	ename = m_tp->elementName(m);
	ge = LookupGe(ename);
	totalSum += na * ge;
      }
    }
    // Add in the charge
    if (m_charge_j != 0.0) {
      ename = "H";
      ge = LookupGe(ename);
      totalSum -= m_charge_j * ge;
    }
    // Ok, now do the calculation. Convert to joules kmol-1
    double dg = m_deltaG_formation_tr_pr * 4.184 * 1.0E3;
    //! Store the result into an internal variable.
    m_Mu0_tr_pr = dg + totalSum;
  }

}
