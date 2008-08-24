/*
 * $Id$
 */
#include "ct_defs.h"
#include "xml.h"
#include "ctml.h"
#include "PDSS_HKFT.h"
#include "WaterProps.h"
#include "PDSS_Water.h"

#include "VPStandardStateTP.h"

using namespace std;

namespace Cantera {
  /**
   * Basic list of constructors and duplicators
   */



  PDSS_HKFT::PDSS_HKFT(VPStandardStateTP *tp, int spindex) :
    PDSS(tp, spindex),
    m_waterSS(0),
    m_pres(OneAtm),
    m_densWaterSS(-1.0),
    m_waterProps(0),
    m_born_coeff_j(-1.0),
    m_r_e_j(-1.0),
    m_deltaG_formation_tr_pr(0.0),
    m_deltaH_formation_tr_pr(0.0),
    m_Mu0_tr_pr(0.0),
    m_Entrop_tr_pr(0.0),
    m_a1(0.0),
    m_a2(0.0),
    m_a3(0.0),
    m_a4(0.0),
    m_c1(0.0),
    m_c2(0.0),
    m_omega_pr_tr(0.0),
    m_Y_pr_tr(0.0),
    m_Z_pr_tr(0.0),
    m_presR_bar(0.0),
    m_charge_j(0.0)
  {
    m_pdssType = cPDSS_MOLAL_HKFT;
  }


  PDSS_HKFT::PDSS_HKFT(VPStandardStateTP *tp, int spindex, std::string inputFile, std::string id) :
    PDSS(tp, spindex),
    m_waterSS(0),
    m_pres(OneAtm),
    m_densWaterSS(-1.0),
    m_waterProps(0),
    m_born_coeff_j(-1.0),
    m_r_e_j(-1.0),
    m_deltaG_formation_tr_pr(0.0),
    m_deltaH_formation_tr_pr(0.0),
    m_Mu0_tr_pr(0.0),
    m_Entrop_tr_pr(0.0),
    m_a1(0.0),
    m_a2(0.0),
    m_a3(0.0),
    m_a4(0.0),
    m_c1(0.0),
    m_c2(0.0),
    m_omega_pr_tr(0.0),
    m_Y_pr_tr(0.0),
    m_Z_pr_tr(0.0),
    m_presR_bar(0.0),
    m_charge_j(0.0)
  {
    m_pdssType = cPDSS_MOLAL_HKFT;
    constructPDSSFile(tp, spindex, inputFile, id);
  }

  PDSS_HKFT::PDSS_HKFT(VPStandardStateTP *tp, int spindex, const XML_Node& speciesNode,
                       const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex),
    m_waterSS(0),
    m_pres(OneAtm),
    m_densWaterSS(-1.0),
    m_waterProps(0),
    m_born_coeff_j(-1.0),
    m_r_e_j(-1.0),
    m_deltaG_formation_tr_pr(0.0),
    m_deltaH_formation_tr_pr(0.0),
    m_Mu0_tr_pr(0.0),
    m_Entrop_tr_pr(0.0),
    m_a1(0.0),
    m_a2(0.0),
    m_a3(0.0),
    m_a4(0.0),
    m_c1(0.0),
    m_c2(0.0),
    m_omega_pr_tr(0.0),
    m_Y_pr_tr(0.0),
    m_Z_pr_tr(0.0),
    m_presR_bar(0.0),
    m_charge_j(0.0)
  {
    m_pdssType = cPDSS_MOLAL_HKFT;
    // We have to read the info from here
    constructPDSSXML(tp, spindex, speciesNode, phaseRoot, spInstalled); 
  }

  PDSS_HKFT::PDSS_HKFT(const PDSS_HKFT &b) :
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
  PDSS_HKFT& PDSS_HKFT::operator=(const PDSS_HKFT&b) {
    if (&b == this) return *this;
    m_tp = b.m_tp;
    m_spindex = b.m_spindex;
    m_temp = b.m_temp;
    m_pres = b.m_pres;
    m_mw = b.m_mw;
    return *this;
  }
  
  /**
   * Destructor for the PDSS_HKFT class
   */
  PDSS_HKFT::~PDSS_HKFT() { 
    delete m_waterProps;
  }
  
  //! Duplicator
  PDSS* PDSS_HKFT::duplMyselfAsPDSS() const {
    PDSS_HKFT * idg = new PDSS_HKFT(*this);
    return (PDSS *) idg;
  }

  /**
   * Return the molar enthalpy in units of J kmol-1
   */
  doublereal 
  PDSS_HKFT::enthalpy_mole() const {
    // Ok we may change this evaluation method in the future.
    double GG = gibbs_mole();
    double SS = entropy_mole();
    double h = GG - m_temp * SS;
    return h;
  }

  doublereal 
  PDSS_HKFT::enthalpy_RT() const {
    double hh = enthalpy_mole();
    double RT = GasConstant * m_temp;
    return hh / RT;
  }

  /**
   * Calculate the internal energy in mks units of
   * J kmol-1 
   */
  doublereal 
  PDSS_HKFT::intEnergy_mole() const {
    double hh = enthalpy_RT();
    double mv = molarVolume();
    return (hh - mv * m_pres);
  }

  /**
   * Calculate the entropy in mks units of 
   * J kmol-1 K-1
   */
  doublereal
  PDSS_HKFT::entropy_mole() const {
    double delS = deltaS();
    return (m_Entrop_tr_pr * 1.0E3 * 4.184 + delS);
  }

  /**
   * Calculate the Gibbs free energy in mks units of
   * J kmol-1
   */
  doublereal PDSS_HKFT::gibbs_mole() const {
    double delG = deltaG();
    return (m_Mu0_tr_pr + delG);
  }

  /**
   * Calculate the constant pressure heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal 
  PDSS_HKFT::cp_mole() const {


    double pbar = m_pres * 1.0E-5;
    double m_presR_bar = OneAtm * 1.0E-5;

    double c1term = m_c1;

    double c2term = m_c2 / (m_temp - 228.) / (m_temp - 228.);
    
    double a3term = m_a3 / (m_temp - 228.) / (m_temp - 228.) / (m_temp - 228.) * 2.0 * m_temp * (m_pres - OneAtm);

    double a4term = m_a4 / (m_temp - 228.) / (m_temp - 228.) / (m_temp - 228.) * 2.0 * m_temp 
      * log((2600. + pbar)/(2600. + m_presR_bar));

    double nu = 166027;
    double r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
  
    double gval = gstar(m_temp, m_pres, 0);

    double dgvaldT = gstar(m_temp, m_pres, 1);
    double d2gvaldT2 = gstar(m_temp, m_pres, 2);

    double r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
    double dr_e_jdT = fabs(m_charge_j) * dgvaldT;

    double omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval)  );

    double domega_jdT = - 2.0 * nu * (m_charge_j * m_charge_j * m_charge_j * m_charge_j / (r_e_j * r_e_j* r_e_j) 
				      - m_charge_j / (3.082 + gval) / (3.082 + gval)  / (3.082 + gval)) * dgvaldT * dgvaldT
      - nu * (m_charge_j * m_charge_j * fabs(m_charge_j) / (r_e_j * r_e_j) 
	      - m_charge_j / (3.082 + gval) / (3.082 + gval)) * d2gvaldT2;


    double d2omega_jdT2 =   nu * (m_charge_j * m_charge_j / (r_e_j * r_e_j) * dr_e_jdT)
      + nu * m_charge_j / (3.082 + gval) / (3.082 + gval) * dgvaldT;
    
    double relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    double drelepsilondT = m_waterProps->relEpsilon(m_temp, m_pres, 1);

    double Y = drelepsilondT / (relepsilon * relepsilon);

    double d2relepsilondT2 = m_waterProps->relEpsilon(m_temp, m_pres, 2);

    double X = d2relepsilondT2 / (relepsilon* relepsilon) - 2.0 * relepsilon * Y * Y;

    double Z = -1.0 / relepsilon;

    double yterm = 2.0 * m_temp * Y * domega_jdT;

    double xterm = omega_j * m_temp * X;

    double otterm = m_temp * d2omega_jdT2 * (Z + 1.0);
 
    double Cp_calgmol = c1term + c2term + a3term + a4term + yterm + xterm + otterm;

    // Convert to Joules / kmol
    double Cp = Cp_calgmol * 1.0E3 * 4.184;
    return Cp;

   
  }

  /**
   * Calculate the constant volume heat capacity
   * in mks units of J kmol-1 K-1
   */
  doublereal 
  PDSS_HKFT::cv_mole() const {
    throw CanteraError("PDSS_HKFT::cv_mole()", "unimplemented");
    return (0.0);
  }

  doublereal 
  PDSS_HKFT::molarVolume() const {
   
    double a1term = m_a1;

    double a2term = m_a2 / (2600.E5 + m_pres);

    double a3term = m_a3 / (m_temp - 228.);

    double a4term = m_a4 / (m_temp - 228.) / (2600.E5 + m_pres);

    double nu = 166027;
    double r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
  
    double gval = gstar(m_temp, m_pres, 0);
    double dgvaldP = gstar(m_temp, m_pres, 3);

    double r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;

    double omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval)  );

    double relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);

    double dr_e_jdP = fabs(m_charge_j) * dgvaldP;

    double domega_jdP = -  nu * (m_charge_j * m_charge_j / (r_e_j * r_e_j) * dr_e_jdP)
      + nu * m_charge_j / (3.082 + gval) / (3.082 + gval) * dgvaldP;
    
    double drelepsilondP = m_waterProps->relEpsilon(m_temp, m_pres, 3);

    double Q = drelepsilondP / (relepsilon * relepsilon);

    double Z = -1.0 / relepsilon;

    double wterm = - domega_jdP * (Z + 1.0);

    double qterm = - omega_j * Q;

    double molVol_calgmolbar = a1term + a2term +  a3term + a4term + wterm + qterm;

    // Convert to m**3 / kmol
    double molVol = molVol_calgmolbar * 4.184 / 100.;
    return molVol;
  }

  doublereal 
  PDSS_HKFT::gibbs_RT_ref() const {
    double m_psave = m_pres;
    m_pres = OneAtm;
    double ee = gibbs_RT();
    m_pres = m_psave;
    return ee;
  }

  doublereal 
  PDSS_HKFT::enthalpy_RT_ref() const {
    double m_psave = m_pres;
    m_pres = OneAtm;
    double hh = enthalpy_RT();
    m_pres = m_psave;
    return hh;
  }

  doublereal 
  PDSS_HKFT::entropy_R_ref() const {
    double m_psave = m_pres;
    m_pres = OneAtm;
    double ee = entropy_R();
    m_pres = m_psave;
    return ee;
  }

  doublereal 
  PDSS_HKFT::cp_R_ref() const {
    double m_psave = m_pres;
    m_pres = OneAtm;
    double ee = cp_R();
    m_pres = m_psave;
    return ee;
  }
  
  doublereal
  PDSS_HKFT::molarVolume_ref() const {
    double m_psave = m_pres;
    m_pres = OneAtm;
    double ee = molarVolume();
    m_pres = m_psave;
    return ee;
  }
  
  /**
   * Calculate the pressure (Pascals), given the temperature and density
   *  Temperature: kelvin
   *  rho: density in kg m-3
   */
  doublereal
  PDSS_HKFT::pressure() const {
    return m_pres;
  }
        
  void 
  PDSS_HKFT::setPressure(doublereal p) {
    m_pres = p;
  }
 
  void PDSS_HKFT::setTemperature(doublereal temp) {
    m_temp = temp;
  }

  doublereal PDSS_HKFT::temperature() const {
    return m_temp;
  }

  void PDSS_HKFT::setState_TP(doublereal temp, doublereal pres) {
    setTemperature(temp);
    setPressure(pres);
  }

  /// critical temperature 
  doublereal 
  PDSS_HKFT::critTemperature() const { 
    throw CanteraError("PDSS_HKFT::critTemperature()", "unimplemented");
    return (0.0);
  }
        
  /// critical pressure
  doublereal PDSS_HKFT::critPressure() const {
    throw CanteraError("PDSS_HKFT::critPressure()", "unimplemented");
    return (0.0);
  }
        
  /// critical density
  doublereal PDSS_HKFT::critDensity() const {
    throw CanteraError("PDSS_HKFT::critDensity()", "unimplemented");
    return (0.0);
  }
        

  void PDSS_HKFT::initThermo() {
    PDSS::initThermo();
    SpeciesThermo &sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);

    m_waterSS = (PDSS_Water *) m_tp->providePDSS(0);
    /*
     *  Section to initialize  m_Z_pr_tr and   m_Y_pr_tr
     */
    double temp = 273.15 + 25.;
    double pres = OneAtm;
    double relepsilon = m_waterProps->relEpsilon(temp, pres, 0);

    m_waterSS->setState_TP(temp, pres);
    m_densWaterSS = m_waterSS->density();
    m_Z_pr_tr = -1.0 / relepsilon;
    //double m_Z_pr_tr = -0.0127803;
    //printf("m_Z_pr_tr = %20.10g\n",  m_Z_pr_tr ); 
    double drelepsilondT = m_waterProps->relEpsilon(temp, pres, 1);
    //double m_Y_pr_tr = -5.799E-5;
    m_Y_pr_tr = drelepsilondT / (relepsilon * relepsilon);
    //printf("m_Y_pr_tr = %20.10g\n",  m_Y_pr_tr );

    m_presR_bar = OneAtm / 1.0E5;
    m_charge_j = m_tp->charge(m_spindex);
    convertDGFormation();

    m_waterProps = new WaterProps(m_waterSS);


  }


  void PDSS_HKFT::initThermoXML(const XML_Node& phaseNode, std::string id) {
    PDSS::initThermoXML(phaseNode, id);
  }

  void PDSS_HKFT::constructPDSSXML(VPStandardStateTP *tp, int spindex,
				   const XML_Node& speciesNode, 
				   const XML_Node& phaseNode, bool spInstalled) {
    PDSS::initThermo();
    SpeciesThermo &sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);

    if (!spInstalled) {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", "spInstalled false not handled");
    }

    const XML_Node *tn = speciesNode.findByName("thermo");
    if (!tn) {
      throw CanteraError("PDSS_HKFT::constructPDSSXML",
			 "no thermo Node for species " + speciesNode.name());
    }
    std::string model = lowercase((*tn)["model"]);
    if (model != "hkft") {
      throw CanteraError("PDSS_HKFT::initThermoXML",
			 "thermo model for species isn't hkft: "
			 + speciesNode.name());
    }
    const XML_Node *hh = tn->findByName("HKFT");
    if (!hh) {
      throw CanteraError("PDSS_HKFT::constructPDSSXML",
			 "no Thermo::HKFT Node for species " + speciesNode.name());
    }

    if (hh->hasChild("DG0_f_Pr_Tr")) {
      double val = getFloat(*hh, "DG0_f_Pr_Tr");
      m_deltaG_formation_tr_pr = val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing DG0_f_Pr_Tr field");
    }

    if (hh->hasChild("DH0_f_Pr_Tr")) {
      double val = getFloat(*hh, "DH0_f_Pr_Tr");
      m_deltaH_formation_tr_pr = val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing DH0_f_Pr_Tr field");
    }

    if (hh->hasChild("S0_Pr_Tr")) {
      double val = getFloat(*hh, "S0_Pr_Tr");
      m_Entrop_tr_pr= val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing S0_Pr_Tr field");
    }

    const XML_Node *ss = speciesNode.findByName("standardState");
    if (!ss) {
      throw CanteraError("PDSS_HKFT::constructPDSSXML",
			 "no standardState Node for species " + speciesNode.name());
    }
    model = lowercase((*ss)["model"]);
    if (model != "hkft") {
      throw CanteraError("PDSS_HKFT::initThermoXML",
			 "standardState model for species isn't hkft: "
			 + speciesNode.name());
    }
    if (ss->hasChild("a1")) {
      double val = getFloat(*ss, "a1");
      m_a1 = val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing a1 field");
    }
    if (ss->hasChild("a2")) {
      double val = getFloat(*ss, "a2");
      m_a2 = val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing a2 field");
    }
    if (ss->hasChild("a3")) {
      double val = getFloat(*ss, "a3");
      m_a3 = val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing a3 field");
    }
    if (ss->hasChild("a4")) {
      double val = getFloat(*ss, "a4");
      m_a4 = val;
    } else {
      throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing a4 field");
    }
    
   if (ss->hasChild("c1")) {
      double val = getFloat(*ss, "c1");
      m_c1 = val;
   } else {
     throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing c1 field");
   }
   if (ss->hasChild("c2")) {
      double val = getFloat(*ss, "c2");
      m_c2 = val;
   } else {
     throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing c2 field");
   }
   if (ss->hasChild("omega_Pr_Tr")) {
      double val = getFloat(*ss, "omega_Pr_Tr");
      m_omega_pr_tr = val;
   } else {
     throw CanteraError("PDSS_HKFT::constructPDSSXML", " missing omega_Pr_Tr field");
   }
    
   
    std::string id = "";
    initThermoXML(phaseNode, id);
  }

  void PDSS_HKFT::constructPDSSFile(VPStandardStateTP *tp, int spindex,
				    std::string inputFile, std::string id) {

    if (inputFile.size() == 0) {
      throw CanteraError("PDSS_HKFT::initThermo",
			 "input file is null");
    }
    std::string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
      throw CanteraError("PDSS_HKFT::initThermo","could not open "
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
      throw CanteraError("PDSS_HKFT::initThermo",
			 "ERROR: Can not find phase named " +
			 id + " in file named " + inputFile);
    }

    XML_Node& speciesList = fxml_phase->child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
					 &(fxml_phase->root()));
    const vector<string>&sss = tp->speciesNames();
    const XML_Node* s =  speciesDB->findByAttr("name", sss[spindex]);

    constructPDSSXML(tp, spindex, *s, *fxml_phase, true);
    delete fxml;
  }








  double PDSS_HKFT::deltaG() const {
    
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

    double relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);

    double Z = -1.0 / relepsilon;

    double wterm = - omega_j * (Z + 1.0);

    double wrterm = m_omega_pr_tr * (m_Z_pr_tr + 1.0);

    double yterm = m_omega_pr_tr * m_Y_pr_tr * (m_temp - 298.15);

    double deltaG_calgmol = sterm + c1term + a1term + a2term + c2term + a3term + a4term + wterm + wrterm + yterm;

    // Convert to Joules / kmol
    double deltaG = deltaG_calgmol * 1.0E3 * 4.184;
    return deltaG;
  }


  double PDSS_HKFT::deltaS() const {
    
    double pbar = m_pres * 1.0E-5;
    double m_presR_bar = OneAtm * 1.0E-5;

    double c1term = m_c1 * log(m_temp/298.15);

    double c2term = -m_c2 / 228. * (( 1.0/(m_temp - 228.) - 1.0/(298.15 - 228.) )
				    + 1.0 / 228. * log( (298.15*(m_temp-228.)) / (m_temp*(298.15-228.)) ));
    
    double a3term = m_a3 / (m_temp - 228.) / (m_temp - 228.) * (pbar - m_presR_bar);

    double a4term = m_a4 / (m_temp - 228.) / (m_temp - 228.) * log((2600. + pbar)/(2600. + m_presR_bar));

    double nu = 166027;
    double r_e_j_pr_tr = m_charge_j * m_charge_j / (m_omega_pr_tr/nu + m_charge_j/3.082);
  
    double gval = gstar(m_temp, m_pres, 0);

    double dgvaldT = gstar(m_temp, m_pres, 1);

    double r_e_j = r_e_j_pr_tr + fabs(m_charge_j) * gval;
    double dr_e_jdT = fabs(m_charge_j) * dgvaldT;

    double omega_j = nu * (m_charge_j * m_charge_j / r_e_j - m_charge_j / (3.082 + gval)  );

    double domega_jdT = -  nu * (m_charge_j * m_charge_j / (r_e_j * r_e_j) * dr_e_jdT)
      + nu * m_charge_j / (3.082 + gval) / (3.082 + gval) * dgvaldT;
    
    double relepsilon = m_waterProps->relEpsilon(m_temp, m_pres, 0);
    double drelepsilondT = m_waterProps->relEpsilon(m_temp, m_pres, 1);

    double Y = drelepsilondT / (relepsilon * relepsilon);

    double Z = -1.0 / relepsilon;

    double wterm = omega_j * Y;

    double wrterm = - m_omega_pr_tr * m_Y_pr_tr;

    double otterm = domega_jdT * (Z + 1.0);

    double domega_jdT_prtr = 0.0;
    double otrterm = - domega_jdT_prtr * m_Z_pr_tr + 1.0;

   
    double deltaS_calgmol = c1term + c2term + a3term + a4term + wterm + wrterm  + otterm + otrterm;

    // Convert to Joules / kmol
    double deltaS = deltaS_calgmol * 1.0E3 * 4.184;
    return deltaS;
  }
   
  double PDSS_HKFT::electrostatic_radii_calc() {
    return 0.0;
  }



  //! Internal formula for the calculation of a_g()
  /*
   * The output of this is in units of Angstroms
   */
  double PDSS_HKFT::ag(const double temp, const int ifunc) const {
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
  double PDSS_HKFT::bg(const double temp, const int ifunc) const {
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


  double PDSS_HKFT::f(const double temp, const double pres, const int ifunc) const {
    
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


  double PDSS_HKFT::g(const double temp, const double pres, const int ifunc) const {
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


 double PDSS_HKFT::gstar(const double temp, const double pres, const int ifunc) const {
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
    {"K",   -19.26943E6}, // Webbook
    {"E",    0.0}         // Don't overcount
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
  double PDSS_HKFT::LookupGe(const std::string& s) {
    int num = sizeof(geDataTable) / sizeof(struct GeData);
    string s3 = s.substr(0,3);
    for (int i = 0; i < num; i++) {
      //if (!std::strncmp(s.c_str(), aWTable[i].name, 3)) {
      if (s3 == geDataTable[i].name) {
        return (geDataTable[i].GeValue);
      }
    }
    throw CanteraError("LookupGe", "element " + s + " not found");
    return -1.0;
  }

 void PDSS_HKFT::convertDGFormation() {
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
