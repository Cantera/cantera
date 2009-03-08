/**
 *  @file HMWSoln.cpp
 *
 * Member functions of Pitzer activity coefficient implementation.
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
/*
 * $Id: HMWSoln.cpp,v 1.5 2006/07/13 20:05:11 hkmoffa Exp $
 */

#ifndef MAX
#define MAX(x,y)    (( (x) > (y) ) ? (x) : (y))
#endif

#include "HMWSoln.h"
#include "importCTML.h"
#include "WaterProps.h"
#include "WaterPDSS.h"

namespace Cantera {

  /**
   * Default constructor
   */
  HMWSoln::HMWSoln() :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_Pcurrent(OneAtm),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_debugCalc(0)
  {
    for (int i = 0; i < 17; i++) {
      elambda[i] = 0.0;
      elambda1[i] = 0.0;
    }
  }
  /**
   * Working constructors
   *
   *  The two constructors below are the normal way
   *  the phase initializes itself. They are shells that call
   *  the routine initThermo(), with a reference to the
   *  XML database to get the info for the phase.
   */
  HMWSoln::HMWSoln(string inputFile, string id) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_Pcurrent(OneAtm),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_debugCalc(0)
  {
    for (int i = 0; i < 17; i++) {
      elambda[i] = 0.0;
      elambda1[i] = 0.0;
    }
    constructPhaseFile(inputFile, id);
  }

  HMWSoln::HMWSoln(XML_Node& phaseRoot, string id) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_Pcurrent(OneAtm),
    m_IionicMolality(0.0),
    m_maxIionicStrength(100.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_debugCalc(0)
  {
    for (int i = 0; i < 17; i++) {
      elambda[i] = 0.0;
      elambda1[i] = 0.0;
    }
    constructPhaseXML(phaseRoot, id);
  }
 
  /**
   * Copy Constructor:
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working copy constructor
   */
  HMWSoln::HMWSoln(const HMWSoln &b) :
    MolalityVPSSTP(b)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = b;
  }

  /**
   * operator=()
   *
   *  Note this stuff will not work until the underlying phase
   *  has a working assignment operator
   */
  HMWSoln& HMWSoln::
  operator=(const HMWSoln &b) {
    if (&b != this) {
      MolalityVPSSTP::operator=(b);
      m_formPitzer          = b.m_formPitzer;
      m_formPitzerTemp      = b.m_formPitzerTemp;
      m_formGC              = b.m_formGC;
      m_Pcurrent            = b.m_Pcurrent;
      m_Aionic              = b.m_Aionic;
      m_IionicMolality      = b.m_IionicMolality;
      m_maxIionicStrength   = b.m_maxIionicStrength;
      m_TempPitzerRef       = b.m_TempPitzerRef;
      m_IionicMolalityStoich= b.m_IionicMolalityStoich;
      m_form_A_Debye        = b.m_form_A_Debye;
      m_A_Debye             = b.m_A_Debye;
      if (m_waterSS) {
	delete m_waterSS;
	m_waterSS = 0;
      }
      if (b.m_waterSS) {
	m_waterSS = new WaterPDSS(*(b.m_waterSS));
      }
      m_densWaterSS         = b.m_densWaterSS;
      if (m_waterProps) {
	delete m_waterProps;
	m_waterProps = 0;
      }
      if (b.m_waterProps) {
	m_waterProps = new WaterProps(*(b.m_waterProps));
      }
      m_expg0_RT            = b.m_expg0_RT;
      m_pe                  = b.m_pe;
      m_pp                  = b.m_pp;
      m_tmpV                = b.m_tmpV;
      m_speciesCharge_Stoich= b.m_speciesCharge_Stoich;
      m_Beta0MX_ij          = b.m_Beta0MX_ij;
      m_Beta0MX_ij_L        = b.m_Beta0MX_ij_L;
      m_Beta0MX_ij_LL       = b.m_Beta0MX_ij_LL;
      m_Beta0MX_ij_P        = b.m_Beta0MX_ij_P;
      m_Beta0MX_ij_coeff    = b.m_Beta0MX_ij_coeff;
      m_Beta1MX_ij          = b.m_Beta1MX_ij;
      m_Beta1MX_ij_L        = b.m_Beta1MX_ij_L;
      m_Beta1MX_ij_LL       = b.m_Beta1MX_ij_LL;
      m_Beta1MX_ij_P        = b.m_Beta1MX_ij_P;
      m_Beta1MX_ij_coeff    = b.m_Beta1MX_ij_coeff;
      m_Beta2MX_ij          = b.m_Beta2MX_ij;
      m_Beta2MX_ij_L        = b.m_Beta2MX_ij_L;
      m_Beta2MX_ij_LL       = b.m_Beta2MX_ij_LL;
      m_Beta2MX_ij_P        = b.m_Beta2MX_ij_P;
      m_Alpha1MX_ij         = b.m_Alpha1MX_ij;
      m_CphiMX_ij           = b.m_CphiMX_ij;
      m_CphiMX_ij_L         = b.m_CphiMX_ij_L;
      m_CphiMX_ij_LL        = b.m_CphiMX_ij_LL;
      m_CphiMX_ij_P         = b.m_CphiMX_ij_P;
      m_CphiMX_ij_coeff     = b.m_CphiMX_ij_coeff;
      m_Theta_ij            = b.m_Theta_ij;
      m_Theta_ij_L          = b.m_Theta_ij_L;
      m_Theta_ij_LL         = b.m_Theta_ij_LL;
      m_Theta_ij_P          = b.m_Theta_ij_P;
      m_Psi_ijk             = b.m_Psi_ijk;
      m_Psi_ijk_L           = b.m_Psi_ijk_L;
      m_Psi_ijk_LL          = b.m_Psi_ijk_LL;
      m_Psi_ijk_P           = b.m_Psi_ijk_P;
      m_Lambda_ij           = b.m_Lambda_ij;
      m_Lambda_ij_L         = b.m_Lambda_ij_L;
      m_Lambda_ij_LL        = b.m_Lambda_ij_LL;
      m_Lambda_ij_P         = b.m_Lambda_ij_P;
      m_lnActCoeffMolal     = b.m_lnActCoeffMolal;
      m_dlnActCoeffMolaldT  = b.m_dlnActCoeffMolaldT;
      m_d2lnActCoeffMolaldT2= b.m_d2lnActCoeffMolaldT2;
      m_dlnActCoeffMolaldP  = b.m_dlnActCoeffMolaldP;

      m_gfunc_IJ            = b.m_gfunc_IJ;
      m_hfunc_IJ            = b.m_hfunc_IJ;
      m_BMX_IJ              = b.m_BMX_IJ;
      m_BMX_IJ_L            = b.m_BMX_IJ_L;
      m_BMX_IJ_LL           = b.m_BMX_IJ_LL;
      m_BMX_IJ_P            = b.m_BMX_IJ_P;
      m_BprimeMX_IJ         = b.m_BprimeMX_IJ;
      m_BprimeMX_IJ_L       = b.m_BprimeMX_IJ_L;
      m_BprimeMX_IJ_LL      = b.m_BprimeMX_IJ_LL;
      m_BprimeMX_IJ_P       = b.m_BprimeMX_IJ_P;
      m_BphiMX_IJ           = b.m_BphiMX_IJ;
      m_BphiMX_IJ_L         = b.m_BphiMX_IJ_L;
      m_BphiMX_IJ_LL        = b.m_BphiMX_IJ_LL;
      m_BphiMX_IJ_P         = b.m_BphiMX_IJ_P;
      m_Phi_IJ              = b.m_Phi_IJ;
      m_Phi_IJ_L            = b.m_Phi_IJ_L;
      m_Phi_IJ_LL           = b.m_Phi_IJ_LL;
      m_Phi_IJ_P            = b.m_Phi_IJ_P;
      m_Phiprime_IJ         = b.m_Phiprime_IJ;
      m_PhiPhi_IJ           = b.m_PhiPhi_IJ;
      m_PhiPhi_IJ_L         = b.m_PhiPhi_IJ_L;
      m_PhiPhi_IJ_LL        = b.m_PhiPhi_IJ_LL;
      m_PhiPhi_IJ_P         = b.m_PhiPhi_IJ_P;
      m_CMX_IJ              = b.m_CMX_IJ;
      m_CMX_IJ_L            = b.m_CMX_IJ_L;
      m_CMX_IJ_LL           = b.m_CMX_IJ_LL;
      m_CMX_IJ_P            = b.m_CMX_IJ_P;
      m_gamma               = b.m_gamma;

      m_CounterIJ           = b.m_CounterIJ;
      m_debugCalc           = b.m_debugCalc;
    }
    return *this;
  }


  /**
   * Test matrix for this object
   *
   *
   *  test problems:
   *  1 = NaCl problem - 5 species -
   *   the thermo is read in from an XML file
   *
   * speci   molality                        charge
   *  Cl-     6.0954          6.0997E+00      -1
   *  H+      1.0000E-08      2.1628E-09      1
   *  Na+     6.0954E+00      6.0997E+00      1
   *  OH-     7.5982E-07      1.3977E-06     -1
   *  HMW_params____beta0MX__beta1MX__beta2MX__CphiMX_____alphaMX__thetaij
   * 10
   * 1  2          0.1775  0.2945   0.0      0.00080    2.0      0.0
   * 1  3          0.0765  0.2664   0.0      0.00127    2.0      0.0
   * 1  4          0.0     0.0      0.0      0.0        0.0     -0.050
   * 2  3          0.0     0.0      0.0      0.0        0.0      0.036
   * 2  4          0.0     0.0      0.0      0.0        0.0      0.0
   * 3  4          0.0864  0.253    0.0      0.0044     2.0      0.0
   * Triplet_interaction_parameters_psiaa'_or_psicc'
   * 2
   * 1  2  3   -0.004
   * 1  3  4   -0.006
   */
  HMWSoln::HMWSoln(int testProb) :
    MolalityVPSSTP(),
    m_formPitzer(PITZERFORM_BASE),
    m_formPitzerTemp(PITZER_TEMP_CONSTANT),
    m_formGC(2),
    m_Pcurrent(OneAtm),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_TempPitzerRef(298.15),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_WATER),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0),
    m_debugCalc(0)
  {
    if (testProb != 1) {
      printf("unknown test problem\n");
      exit(-1);
    }

    constructPhaseFile("HMW_NaCl.xml", "");

    int i = speciesIndex("Cl-");
    int j = speciesIndex("H+");
    int n = i * m_kk + j;
    int ct = m_CounterIJ[n];
    m_Beta0MX_ij[ct] = 0.1775;
    m_Beta1MX_ij[ct] = 0.2945;
    m_CphiMX_ij[ct]  = 0.0008;
    m_Alpha1MX_ij[ct]= 2.000;


    i = speciesIndex("Cl-");
    j = speciesIndex("Na+");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Beta0MX_ij[ct] = 0.0765;
    m_Beta1MX_ij[ct] = 0.2664;
    m_CphiMX_ij[ct]  = 0.00127;
    m_Alpha1MX_ij[ct]= 2.000;


    i = speciesIndex("Cl-");
    j = speciesIndex("OH-");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Theta_ij[ct] = -0.05;

    i = speciesIndex("H+");
    j = speciesIndex("Na+");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Theta_ij[ct] = 0.036;

    i = speciesIndex("Na+");
    j = speciesIndex("OH-");
    n = i * m_kk + j;
    ct = m_CounterIJ[n];
    m_Beta0MX_ij[ct] = 0.0864;
    m_Beta1MX_ij[ct] = 0.253;
    m_CphiMX_ij[ct]  = 0.0044;
    m_Alpha1MX_ij[ct]= 2.000;

    i = speciesIndex("Cl-");
    j = speciesIndex("H+");
    int k = speciesIndex("Na+");
    double param = -0.004;
    n = i * m_kk *m_kk + j * m_kk + k ;
    m_Psi_ijk[n] = param;
    n = i * m_kk *m_kk + k * m_kk + j ;
    m_Psi_ijk[n] = param;
    n = j * m_kk *m_kk + i * m_kk + k ;
    m_Psi_ijk[n] = param;
    n = j * m_kk *m_kk + k * m_kk + i ;
    m_Psi_ijk[n] = param;
    n = k * m_kk *m_kk + j * m_kk + i ;
    m_Psi_ijk[n] = param;
    n = k * m_kk *m_kk + i * m_kk + j ;
    m_Psi_ijk[n] = param;

    i = speciesIndex("Cl-");
    j = speciesIndex("Na+");
    k = speciesIndex("OH-");
    param = -0.006;
    n = i * m_kk *m_kk + j * m_kk + k ;
    m_Psi_ijk[n] = param;
    n = i * m_kk *m_kk + k * m_kk + j ;
    m_Psi_ijk[n] = param;
    n = j * m_kk *m_kk + i * m_kk + k ;
    m_Psi_ijk[n] = param;
    n = j * m_kk *m_kk + k * m_kk + i ;
    m_Psi_ijk[n] = param;
    n = k * m_kk *m_kk + j * m_kk + i ;
    m_Psi_ijk[n] = param;
    n = k * m_kk *m_kk + i * m_kk + j ;
    m_Psi_ijk[n] = param;

    printCoeffs();
  }

  /**
   * ~HMWSoln():   (virtual)
   *
   *     Destructor: does nothing:
   */
  HMWSoln::~HMWSoln() {
    if (m_waterProps) {
      delete m_waterProps; m_waterProps = 0;
    }
    if (m_waterSS) {
      delete m_waterSS; m_waterSS = 0;
    } 
  }

  /**
   *  duplMyselfAsThermoPhase():
   *
   *  This routine operates at the ThermoPhase level to 
   *  duplicate the current object. It uses the copy constructor
   *  defined above.
   */
  ThermoPhase* HMWSoln::duplMyselfAsThermoPhase() {
    HMWSoln* mtp = new HMWSoln(*this);
    return (ThermoPhase *) mtp;
  }

  /** 
   * Equation of state type flag. The base class returns
   * zero. Subclasses should define this to return a unique
   * non-zero value. Constants defined for this purpose are
   * listed in mix_defs.h.
   */
  int HMWSoln::eosType() const {
    int res;
    switch (m_formGC) {
    case 0:
      res = cHMWSoln0;
      break;
    case 1:
      res = cHMWSoln1;
      break;
    case 2:
      res = cHMWSoln2;
      break;
    default:
      throw CanteraError("eosType", "Unknown type");
      break;
    }
    return res;
  }

  //
  // -------- Molar Thermodynamic Properties of the Solution --------------- 
  //
  /**
   * Molar enthalpy of the solution. Units: J/kmol.
   */
  doublereal HMWSoln::enthalpy_mole() const {
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    getMoleFractions(DATA_PTR(m_pp));
    double val = mean_X(DATA_PTR(m_tmpV));
#ifdef DEBUG_HKM
    double val0 = 0.0;
    for (int k = 0; k < m_kk; k++) {
      val0 += m_tmpV[k] * m_pp[k];
    }
    //if (val != val0) {
    // printf("ERROR\n");
    //}
#endif
    return val;
  }

  doublereal HMWSoln::relative_enthalpy() const {
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    double hbar = mean_X(DATA_PTR(m_tmpV));
    getEnthalpy_RT(DATA_PTR(m_gamma));
    double RT = GasConstant * temperature();
    for (int k = 0; k < m_kk; k++) {
      m_gamma[k] *= RT;
    }
    double h0bar = mean_X(DATA_PTR(m_gamma));
    return (hbar - h0bar);
  }



  doublereal HMWSoln::relative_molal_enthalpy() const {
    double L = relative_enthalpy();
    getMoleFractions(DATA_PTR(m_tmpV));
    double xanion = 0.0;
    int kcation = -1;
    double xcation = 0.0;
    int kanion = -1;
    const double *charge =  DATA_PTR(m_speciesCharge);
    for (int k = 0; k < m_kk; k++) {
      if (charge[k] > 0.0) {
	if (m_tmpV[k] > xanion) {
	  xanion = m_tmpV[k];
	  kanion = k;
	}
      } else if (charge[k] < 0.0) {
	if (m_tmpV[k] > xcation) {
	  xcation = m_tmpV[k];
	  kcation = k;
	}
      }
    }
    if (kcation < 0 || kanion < 0) {
      return L;
    }
    double xuse = xcation;
    int kuse = kcation;
    double factor = 1;
    if (xanion < xcation) {
      xuse = xanion;
      kuse = kanion;
      if (charge[kcation] != 1.0) {
	factor = charge[kcation];
      }
    } else {
      if (charge[kanion] != 1.0) {
	factor = charge[kanion];
      }
    }
    xuse = xuse / factor;
    L = L / xuse;
    return L;
  }

  /**
   * Molar internal energy of the solution. Units: J/kmol.
   *
   * This is calculated from the soln enthalpy and then
   * subtracting pV.
   */
  doublereal HMWSoln::intEnergy_mole() const {
    double hh = enthalpy_mole();
    double pres = pressure();
    double molarV = 1.0/molarDensity();
    double uu = hh - pres * molarV;
    return uu;
  }

  /**
   *  Molar soln entropy at constant pressure. Units: J/kmol/K. 
   *
   *  This is calculated from the partial molar entropies.
   */
  doublereal HMWSoln::entropy_mole() const {
    getPartialMolarEntropies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
  }

  /// Molar Gibbs function. Units: J/kmol. 
  doublereal HMWSoln::gibbs_mole() const {
    getChemPotentials(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
  }

  /** Molar heat capacity at constant pressure. Units: J/kmol/K.
   *
   * Returns the solution heat capacition at constant pressure. 
   * This is calculated from the partial molar heat capacities.
   */ 
  doublereal HMWSoln::cp_mole() const {
    getPartialMolarCp(DATA_PTR(m_tmpV));
    double val = mean_X(DATA_PTR(m_tmpV));
    return val;
  }

  /// Molar heat capacity at constant volume. Units: J/kmol/K. 
  doublereal HMWSoln::cv_mole() const {
    //getPartialMolarCv(m_tmpV.begin());
    //return mean_X(m_tmpV.begin());
    err("not implemented");
    return 0.0;
  }

  //
  // ------- Mechanical Equation of State Properties ------------------------
  //

  /**
   * Pressure. Units: Pa.
   * For this incompressible system, we return the internally storred
   * independent value of the pressure.
   */
  doublereal HMWSoln::pressure() const {
    return m_Pcurrent;
  }

  /**
   * Set the pressure at constant temperature. Units: Pa.
   * This method sets a constant within the object.
   * The mass density is not a function of pressure.
   */
  void HMWSoln::setPressure(doublereal p) {
#ifdef DEBUG_MODE
    //printf("setPressure: %g\n", p);
#endif
    double temp = temperature();
    /*
     * Call the water SS and set it's internal state
     */
    m_waterSS->setTempPressure(temp, p);

    /*
     * Store the internal density of the water SS.
     * Note, we would have to do this for all other
     * species if they had pressure dependent properties.
     */
    m_densWaterSS = m_waterSS->density();
    /*
     * Store the current pressure
     */
    m_Pcurrent = p;
    /*
     * Calculate all of the other standard volumes
     * -> note these are constant for now
     */
    /*
     * Get the partial molar volumes of all of the
     * species. -> note this is a lookup for 
     * water, here since it was done above.
     */
    double *vbar = &m_pp[0];
    getPartialMolarVolumes(vbar);

    /*
     * Get mole fractions of all species.
     */
    double *x = &m_tmpV[0];
    getMoleFractions(x);
	
    /*
     * Calculate the solution molar volume and the 
     * solution density.
     */
    doublereal vtotal = 0.0;
    for (int i = 0; i < m_kk; i++) {
      vtotal += vbar[i] * x[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;

    /*
     * Now, update the State class with the results. This
     * store the denisty.
     */
    State::setDensity(dd);
  }

  /**
   * The isothermal compressibility. Units: 1/Pa.
   * The isothermal compressibility is defined as
   * \f[
   * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
   * \f]
   *
   *  It's equal to zero for this model, since the molar volume
   *  doesn't change with pressure or temperature.
   */
  doublereal HMWSoln::isothermalCompressibility() const {
    throw CanteraError("HMWSoln::isothermalCompressibility",
		       "unimplemented");
    return 0.0;
  }

  /**
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
  doublereal HMWSoln::thermalExpansionCoeff() const {
    throw CanteraError("HMWSoln::thermalExpansionCoeff",
		       "unimplemented");
    return 0.0;
  }
    
  /**
   * Overwritten setDensity() function is necessary because the
   * density is not an indendent variable.
   *
   * This function will now throw an error condition
   *
   * Note, in general, setting the phase density is now a nonlinear
   * calculation. P and T are the fundamental variables. This
   * routine should be revamped to do the nonlinear problem
   *
   * @internal May have to adjust the strategy here to make
   * the eos for these materials slightly compressible, in order
   * to create a condition where the density is a function of
   * the pressure.
   *
   * This function will now throw an error condition.
   *
   *  NOTE: This is an overwritten function from the State.h
   *        class
   */
  void HMWSoln::setDensity(doublereal rho) {
    double dens_old = density();

    if (rho != dens_old) {
      throw CanteraError("HMWSoln::setDensity",
			 "Density is not an independent variable");
    }
  }	

  /**
   * Overwritten setMolarDensity() function is necessary because the
   * density is not an indendent variable.
   *
   * This function will now throw an error condition.
   *
   *  NOTE: This is an overwritten function from the State.h
   *        class
   */
  void HMWSoln::setMolarDensity(doublereal rho) {
    throw CanteraError("HMWSoln::setMolarDensity",
		       "Density is not an independent variable");
  }

  /**
   * Overwritten setTemperature(double) from State.h. This
   * function sets the temperature, and makes sure that
   * the value propagates to underlying objects.
   */
  void HMWSoln::setTemperature(double temp) {
    m_waterSS->setTemperature(temp);
    State::setTemperature(temp);
  }

  //
  // ------- Activities and Activity Concentrations
  //

  /**
   * This method returns an array of generalized concentrations
   * \f$ C_k\f$ that are defined such that 
   * \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$ 
   * is a standard concentration
   * defined below.  These generalized concentrations are used
   * by kinetics manager classes to compute the forward and
   * reverse rates of elementary reactions. 
   *
   * @param c Array of generalized concentrations. The 
   *           units depend upon the implementation of the
   *           reaction rate expressions within the phase.
   */
  void HMWSoln::getActivityConcentrations(doublereal* c) const {
    double c_solvent = standardConcentration();
    getActivities(c);
    for (int k = 0; k < m_kk; k++) {
      c[k] *= c_solvent;
    }
  }

  /**
   * The standard concentration \f$ C^0_k \f$ used to normalize
   * the generalized concentration. In many cases, this quantity
   * will be the same for all species in a phase - for example,
   * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
   * reason, this method returns a single value, instead of an
   * array.  However, for phases in which the standard
   * concentration is species-specific (e.g. surface species of
   * different sizes), this method may be called with an
   * optional parameter indicating the species.
   *
   * For the time being we will use the concentration of pure
   * solvent for the the standard concentration of all species.
   * This has the effect of making reaction rates
   * based on the molality of species proportional to the
   * molality of the species.
   */
  doublereal HMWSoln::standardConcentration(int k) const {
    double mvSolvent = m_speciesSize[m_indexSolvent];
    return 1.0 / mvSolvent;
  }
    
  /**
   * Returns the natural logarithm of the standard 
   * concentration of the kth species
   */
  doublereal HMWSoln::logStandardConc(int k) const {
    double c_solvent = standardConcentration(k);
    return log(c_solvent);
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
  void HMWSoln::getUnitsStandardConc(double *uA, int k, int sizeUA) {
    for (int i = 0; i < sizeUA; i++) {
      if (i == 0) uA[0] = 1.0;
      if (i == 1) uA[1] = -nDim();
      if (i == 2) uA[2] = 0.0;
      if (i == 3) uA[3] = 0.0;
      if (i == 4) uA[4] = 0.0;
      if (i == 5) uA[5] = 0.0;
    }
  }    


  /**
   * Get the array of non-dimensional activities at
   * the current solution temperature, pressure, and
   * solution concentration.
   * (note solvent activity coefficient is on the molar scale).
   *
   */
  void HMWSoln::getActivities(doublereal* ac) const {
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    s_update_lnMolalityActCoeff();
    /*
     * Now calculate the array of activities.
     */
    for (int k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) {
	ac[k] = m_molalities[k] * exp(m_lnActCoeffMolal[k]);
      }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    ac[m_indexSolvent] =
      exp(m_lnActCoeffMolal[m_indexSolvent]) * xmolSolvent;
  }

  /**
   * getMolalityActivityCoefficients()             (virtual, const)
   *
   * Get the array of non-dimensional Molality based
   * activity coefficients at
   * the current solution temperature, pressure, and
   * solution concentration.
   * (note solvent activity coefficient is on the molar scale).
   *
   *  Note, most of the work is done in an internal private routine
   */
  void HMWSoln::
  getMolalityActivityCoefficients(doublereal* acMolality) const {

    A_Debye_TP(-1.0, -1.0);
    s_update_lnMolalityActCoeff();
    copy(m_lnActCoeffMolal.begin(), m_lnActCoeffMolal.end(), acMolality);
    for (int k = 0; k < m_kk; k++) {
      acMolality[k] = exp(acMolality[k]);
    }
  }

  //
  // ------ Partial Molar Properties of the Solution -----------------
  //
  /**
   * Get the species chemical potentials. Units: J/kmol.
   *
   * This function returns a vector of chemical potentials of the 
   * species in solution.
   *
   * \f[
   *    \mu_k = \mu^{o}_k(T,P) + R T ln(m_k)
   * \f]
   *
   * \f[
   *    \mu_solvent = \mu^{o}_solvent(T,P) +
   *            R T ((X_solvent - 1.0) / X_solvent)
   * \f]
   */
  void HMWSoln::getChemPotentials(doublereal* mu) const{
    double xx;
    const double xxSmall = 1.0E-150; 
    /*
     * First get the standard chemical potentials in
     * molar form.
     *  -> this requires updates of standard state as a function
     *     of T and P
     */
    getStandardChemPotentials(mu);
    /*
     * Update the activity coefficients
     * This also updates the internal molality array.
     */
    s_update_lnMolalityActCoeff();
    /*
     *   
     */
    doublereal RT = GasConstant * temperature();
    double xmolSolvent = moleFraction(m_indexSolvent);
    for (int k = 0; k < m_kk; k++) {
      if (m_indexSolvent != k) {
	xx = MAX(m_molalities[k], xxSmall);
	mu[k] += RT * (log(xx) + m_lnActCoeffMolal[k]);
      }
    }
    xx = MAX(xmolSolvent, xxSmall);
    mu[m_indexSolvent] +=  
      RT * (log(xx) + m_lnActCoeffMolal[m_indexSolvent]);
  }


  /**
   * Returns an array of partial molar enthalpies for the species
   * in the mixture.
   * Units (J/kmol)
   *
   * We calculate this quantity partially from the relation and
   * partially by calling the standard state enthalpy function.
   *
   *     hbar_i = - T**2 * d(chemPot_i/T)/dT 
   *
   * We calculate 
   */
  void HMWSoln::getPartialMolarEnthalpies(doublereal* hbar) const {
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
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dT();
    double RTT = RT * T;
    for (int k = 0; k < m_kk; k++) {
      hbar[k] -= RTT * m_dlnActCoeffMolaldT[k];
    }
  }

  /**
   *
   * getPartialMolarEntropies()        (virtual, const)
   *
   * Returns an array of partial molar entropies of the species in the
   * solution. Units: J/kmol.
   *
   * Maxwell's equations provide an insight in how to calculate this
   * (p.215 Smith and Van Ness)
   *
   *      d(chemPot_i)/dT = -sbar_i
   *   
   * Combining this with the expression H = G + TS yields:
   *
   *  \f[
   * \bar s_k(T,P) =  \hat s^0_k(T) - R log(M0 * molality[k] ac[k])
   *                      - R T^2 d log(ac[k]) / dT
   * \f]
   *
   *
   * The reference-state pure-species entropies,\f$ \hat s^0_k(T) \f$,
   * at the reference pressure, \f$ P_{ref} \f$,  are computed by the
   * species thermodynamic
   * property manager. They are polynomial functions of temperature.
   * @see SpeciesThermo
   */
  void HMWSoln::
  getPartialMolarEntropies(doublereal* sbar) const {
    int k;
    /*
     * Get the standard state entropies at the temperature
     * and pressure of the solution.
     */
    getEntropy_R(sbar);
    /*
     * Dimensionalize the entropies
     */
    doublereal R = GasConstant;
    for (k = 0; k < m_kk; k++) {
      sbar[k] *= R;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnMolalityActCoeff();
    /*
     * First we will add in the obvious dependence on the T
     * term out front of the log activity term
     */
    doublereal mm;
    for (k = 0; k < m_kk; k++) {
      if (k != m_indexSolvent) {
	mm = fmaxx(SmallNumber, m_molalities[k]);
	sbar[k] -= R * (log(mm) + m_lnActCoeffMolal[k]);
      }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    mm = fmaxx(SmallNumber, xmolSolvent);
    sbar[m_indexSolvent] -= R *(log(mm) + m_lnActCoeffMolal[m_indexSolvent]);
    /*
     * Check to see whether activity coefficients are temperature
     * dependent. If they are, then calculate the their temperature
     * derivatives and add them into the result.
     */
    s_update_dlnMolalityActCoeff_dT();
    double RT = R * temperature();
    for (k = 0; k < m_kk; k++) {
      sbar[k] -= RT * m_dlnActCoeffMolaldT[k];
    }
  }

  /**
   * getPartialMolarVolumes()                (virtual, const)
   *
   * Returns an array of partial molar volumes of the species
   * in the solution. Units: m^3 kmol-1.
   *
   * For this solution, the partial molar volumes are a
   * complex function of pressure.
   *
   * The general relation is 
   *
   *       vbar_i = d(chemPot_i)/dP at const T, n
   *
   *              = V0_i + d(Gex)/dP)_T,M
   *
   *              = V0_i + RT d(lnActCoeffi)dP _T,M
   *
   */
  void HMWSoln::getPartialMolarVolumes(doublereal* vbar) const {
    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);
    /*
     * Update the derivatives wrt the activity coefficients.
     */
    s_update_lnMolalityActCoeff();
    s_Pitzer_dlnMolalityActCoeff_dP();
    double T = temperature();
    double RT = GasConstant * T;
    for (int k = 0; k < m_kk; k++) {
      vbar[k] += RT * m_dlnActCoeffMolaldP[k];
    }
  }

  /*
   * Partial molar heat capacity of the solution:
   *   The kth partial molar heat capacity is  equal to 
   *   the temperature derivative of the partial molar
   *   enthalpy of the kth species in the solution at constant
   *   P and composition (p. 220 Smith and Van Ness).
   *
   *     Cp = -T d2(chemPot_i)/dT2
   */
  void HMWSoln::getPartialMolarCp(doublereal* cpbar) const {
    /*
     * Get the nondimensional gibbs standard state of the
     * species at the T and P of the solution.
     */
    getCp_R(cpbar);
	
    for (int k = 0; k < m_kk; k++) {
      cpbar[k] *= GasConstant;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally storred molalities.
     */
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dT();
    s_update_d2lnMolalityActCoeff_dT2();
    double T = temperature();
    double RT = GasConstant * T;
    double RTT = RT * T;
    for (int k = 0; k < m_kk; k++) {
      cpbar[k] -= (2.0 * RT * m_dlnActCoeffMolaldT[k] +
		   RTT * m_d2lnActCoeffMolaldT2[k]);
    }
  }

  /*
   * -------- Properties of the Standard State of the Species
   *           in the Solution ------------------
   */

  /**
   *  getStandardChemPotentials()      (virtual, const)
   *
   *
   *  Get the standard state chemical potentials of the species.
   *  This is the array of chemical potentials at unit activity 
   *  (Mole fraction scale)
   *  \f$ \mu^0_k(T,P) \f$.
   *  We define these here as the chemical potentials of the pure
   *  species at the temperature and pressure of the solution.
   *  This function is used in the evaluation of the 
   *  equilibrium constant Kc. Therefore, Kc will also depend
   *  on T and P. This is the norm for liquid and solid systems.
   *
   *  units = J / kmol
   */
  void HMWSoln::getStandardChemPotentials(doublereal* mu) const {
    getGibbs_ref(mu);
    doublereal pref;
    doublereal delta_p;
    for (int k = 1; k < m_kk; k++) {
      pref = m_spthermo->refPressure(k);
      delta_p = m_Pcurrent - pref;
      mu[k] += delta_p * m_speciesSize[k];
    }
    mu[0] = m_waterSS->gibbs_mole();
  }
    
  /**
   * Get the nondimensional gibbs function for the species
   * standard states at the current T and P of the solution.
   *
   *  \f[
   *  \mu^0_k(T,P) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
   * \f]
   * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
   * \f$ \mu^{ref}_k(T)\f$ is the chemical potential of pure
   * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
   *
   * @param grt Vector of length m_kk, which on return sr[k]
   *           will contain the nondimensional 
   *           standard state gibbs function for species k. 
   */
  void HMWSoln::getGibbs_RT(doublereal* grt) const {
    getStandardChemPotentials(grt);
    doublereal invRT = 1.0 / _RT();
    for (int k = 0; k < m_kk; k++) {
      grt[k] *= invRT;
    }
  }
    
  /**
   *
   * getPureGibbs()
   *
   * Get the Gibbs functions for the pure species
   * at the current <I>T</I> and <I>P</I> of the solution.
   * We assume an incompressible constant partial molar
   * volume here:
   * \f[
   *  \mu^0_k(T,p) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
   * \f]
   * where \f$V_k\f$ is the molar volume of pure species <I>k<\I>.
   * \f$ u^{ref}_k(T)\f$ is the chemical potential of pure
   * species <I>k<\I> at the reference pressure, \f$P_{ref}\f$.
   */
  void HMWSoln::getPureGibbs(doublereal* gpure) const {
    getStandardChemPotentials(gpure);
  }

  /**
   *
   * getEnthalpy_RT()        (virtual, const)
   *
   * Get the array of nondimensional Enthalpy functions for the ss
   * species at the current <I>T</I> and <I>P</I> of the solution.
   * We assume an incompressible constant partial molar
   * volume here:
   * \f[
   *  h^0_k(T,P) = h^{ref}_k(T) + (P - P_{ref}) * V_k
   * \f]
   * where \f$V_k\f$ is the molar volume of SS species <I>k<\I>.
   * \f$ h^{ref}_k(T)\f$ is the enthalpy of the SS
   * species <I>k<\I> at the reference pressure, \f$P_{ref}\f$.
   */
  void HMWSoln::
  getEnthalpy_RT(doublereal* hrt) const {
    getEnthalpy_RT_ref(hrt);
    doublereal pref;
    doublereal delta_p;
    double RT = _RT();
    for (int k = 1; k < m_kk; k++) {
      pref = m_spthermo->refPressure(k);
      delta_p = m_Pcurrent - pref;
      hrt[k] += delta_p/ RT * m_speciesSize[k];
    }
    hrt[0] = m_waterSS->enthalpy_mole();
    hrt[0] /= RT;
  }
    
  /**
   *  getEntropy_R()         (virtual, const)
   * 
   * Get the nondimensional Entropies for the species
   * standard states at the current T and P of the solution.
   *
   * Note, this is equal to the reference state entropies
   * due to the zero volume expansivity:
   * i.e., (dS/dp)_T = (dV/dT)_P = 0.0
   *
   * @param sr Vector of length m_kk, which on return sr[k]
   *           will contain the nondimensional
   *           standard state entropy of species k.
   */
  void HMWSoln::
  getEntropy_R(doublereal* sr) const {
    getEntropy_R_ref(sr);
    sr[0] = m_waterSS->entropy_mole();
    sr[0] /= GasConstant;
  }

  /**
   * Get the nondimensional heat capacity at constant pressure
   * function for the species
   * standard states at the current T and P of the solution.
   * \f[
   *  Cp^0_k(T,P) = Cp^{ref}_k(T)
   * \f]
   * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
   * \f$ Cp^{ref}_k(T)\f$ is the constant pressure heat capacity
   * of species <I>k</I> at the reference pressure, \f$p_{ref}\f$.
   *
   * @param cpr Vector of length m_kk, which on return cpr[k]
   *           will contain the nondimensional 
   *           constant pressure heat capacity for species k. 
   */
  void HMWSoln::getCp_R(doublereal* cpr) const {
    getCp_R_ref(cpr); 
    cpr[0] = m_waterSS->cp_mole();
    cpr[0] /= GasConstant;
  }
    
  /**
   * Get the molar volumes of each species in their standard
   * states at the current
   * <I>T</I> and <I>P</I> of the solution.
   * units = m^3 / kmol
   *
   * The water calculation is done separately.
   */
  void HMWSoln::getStandardVolumes(doublereal *vol) const {
    copy(m_speciesSize.begin(),
	 m_speciesSize.end(), vol);
    double dd = m_waterSS->density();
    vol[0] = molecularWeight(0)/dd;
  }

  /*
   * ------ Thermodynamic Values for the Species Reference States ---
   */

  // -> This is handled by VPStandardStatesTP

  /*
   *  -------------- Utilities -------------------------------
   */

  /**
   * @internal
   * Set equation of state parameters. The number and meaning of
   * these depends on the subclass. 
   * @param n number of parameters
   * @param c array of \i n coefficients
   * 
   */
  void HMWSoln::setParameters(int n, doublereal* c) {
  }
  void HMWSoln::getParameters(int &n, doublereal * const c) {
  }
  /**
   * Set equation of state parameter values from XML
   * entries. This method is called by function importPhase in
   * file importCTML.cpp when processing a phase definition in
   * an input file. It should be overloaded in subclasses to set
   * any parameters that are specific to that particular phase
   * model.
   *
   * @param eosdata An XML_Node object corresponding to
   * the "thermo" entry for this phase in the input file.
   *
   * HKM -> Right now, the parameters are set elsewhere (initThermoXML)
   *        It just didn't seem to fit.
   */
  void HMWSoln::setParametersFromXML(const XML_Node& eosdata) {
  }
    
  /**
   * Get the saturation pressure for a given temperature. 
   * Note the limitations of this function. Stability considerations
   * concernting multiphase equilibrium are ignored in this 
   * calculation. Therefore, the call is made directly to the SS of 
   * water underneath. The object is put back into its original
   * state at the end of the call.
   */
  doublereal HMWSoln::satPressure(doublereal t) const {
    double p_old = pressure();
    double t_old = temperature();
    double pres = m_waterSS->satPressure(t);
    /*
     * Set the underlying object back to its original state.
     */
    m_waterSS->setState_TP(t_old, p_old);
    return pres;
  }

  /**
   * Report the molar volume of species k
   *
   * units - \f$ m^3 kmol^-1 \f$
   */
  double HMWSoln::speciesMolarVolume(int k) const {
    double vol = m_speciesSize[k];
    if (k == 0) {
      double dd = m_waterSS->density();
      vol = molecularWeight(0)/dd;
    }
    return vol;
  }
 
  /**
   * A_Debye_TP()                              (virtual)
   *
   *   Returns the A_Debye parameter as a function of temperature
   *  and pressure. This function also sets the internal value
   *  of the parameter within the object, if it is changeable. 
   *
   *  The default is to assume that it is constant, given
   *  in the initialization process and storred in the
   *  member double, m_A_Debye
   *
   *            A_Debye = (1/(8 Pi)) sqrt(2 Na dw /1000) 
   *                          (e e/(epsilon R T))^3/2
   *
   *                    where epsilon = e_rel * e_naught
   *
   * Note, this is si units. Frequently, gaussian units are
   * used in Pitzer's papers where D is used, D = epsilon/(4 Pi)
   * units = A_Debye has units of sqrt(gmol kg-1).
   */
  double HMWSoln::A_Debye_TP(double tempArg, double presArg) const {
    double T = temperature(); 
    double A;
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
      P = presArg;
    }

    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
      A = m_A_Debye;
      break;
    case A_DEBYE_WATER:
      A = m_waterProps->ADebye(T, P, 0);
      m_A_Debye = A;
      break;
    default:
      printf("shouldn't be here\n");
      exit(-1);
    }
    return A;
  }

  /**
   * dA_DebyedT_TP()                              (virtual)
   *
   *  Returns the derivative of the A_Debye parameter with
   *  respect to temperature as a function of temperature
   *  and pressure. 
   *
   * units = A_Debye has units of sqrt(gmol kg-1).
   *         Temp has units of Kelvin.
   */
  double HMWSoln::dA_DebyedT_TP(double tempArg, double presArg) const {
    double T = temperature();
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
      P = presArg;
    }
    double dAdT;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
      dAdT = 0.0;
      break;
    case A_DEBYE_WATER:
      dAdT = m_waterProps->ADebye(T, P, 1);
      //dAdT = WaterProps::ADebye(T, P, 1);
      break;
    default:
      printf("shouldn't be here\n");
      exit(-1);
    }
    return dAdT;
  }

  /**
   * dA_DebyedP_TP()                              (virtual)
   *
   *  Returns the derivative of the A_Debye parameter with
   *  respect to pressure, as a function of temperature
   *  and pressure. 
   *
   * units = A_Debye has units of sqrt(gmol kg-1).
   *         Pressure has units of pascals.
   */
  double HMWSoln::dA_DebyedP_TP(double tempArg, double presArg) const {
    double T = temperature();
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
      P = presArg;
    }
    double dAdP;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
      dAdP = 0.0;
      break;
    case A_DEBYE_WATER:
      dAdP = m_waterProps->ADebye(T, P, 3);
      break;
    default:
      printf("shouldn't be here\n");
      exit(-1);
    }
    return dAdP;
  }


  /**
   *  Calculate the DH Parameter used for the Enthalpy calcalations
   *
   *      ADebye_L = 4 R T**2 d(Aphi) / dT
   *
   *   where   Aphi = A_Debye/3
   *
   *   units -> J / (kmolK) * sqrt( kg/gmol)
   * 
   */
  double HMWSoln::ADebye_L(double tempArg, double presArg) const {
    double dAdT = dA_DebyedT_TP();
    double dAphidT = dAdT /3.0;
    double T = temperature();
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double retn = dAphidT * (4.0 * GasConstant * T * T);
    return retn;
  }

  /**
   *  Calculate the DH Parameter used for the Volume calcalations
   *
   *      ADebye_V = - 4 R T d(Aphi) / dP
   *
   *   where   Aphi = A_Debye/3
   *
   *   units -> J / (kmolK) * sqrt( kg/gmol)
   * 
   */
  double HMWSoln::ADebye_V(double tempArg, double presArg) const {
    double dAdP = dA_DebyedP_TP();
    double dAphidP = dAdP /3.0;
    double T = temperature();
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double retn = - dAphidP * (4.0 * GasConstant * T);
    return retn;
  }

  /**
   * Return Pitzer's definition of A_J. This is basically the
   * temperature derivative of A_L, and the second derivative
   * of Aphi
   * It's the DH parameter used in heat capacity calculations
   *  
   *  A_J = 2 A_L/T + 4 * R * T * T * d2(A_phi)/dT2
   *
   *    Units = sqrt(kg/gmol) (R)
   *  
   *   where
   *      ADebye_L = 4 R T**2 d(Aphi) / dT
   *
   *   where   Aphi = A_Debye/3
   *
   *   units -> J / (kmolK) * sqrt( kg/gmol)
   * 
   */
  double HMWSoln::ADebye_J(double tempArg, double presArg) const {
    double T = temperature();
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double A_L = ADebye_L(T, presArg);
    double d2 = d2A_DebyedT2_TP(T, presArg);
    double d2Aphi = d2 / 3.0;
    double retn = 2.0 * A_L / T + 4.0 * GasConstant * T * T *d2Aphi;
    return retn;
  }

  /**
   * d2A_DebyedT2_TP()                              (virtual)
   *
   *  Returns the 2nd derivative of the A_Debye parameter with
   *  respect to temperature as a function of temperature
   *  and pressure. 
   *
   * units = A_Debye has units of sqrt(gmol kg-1).
   *         Temp has units of Kelvin.
   */
  double HMWSoln::d2A_DebyedT2_TP(double tempArg, double presArg) const {
    double T = temperature();
    if (tempArg != -1.0) {
      T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
      P = presArg;
    }
    double d2AdT2;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
      d2AdT2 = 0.0;
      break;
    case A_DEBYE_WATER:
      d2AdT2 = m_waterProps->ADebye(T, P, 2);
      break;
    default:
      printf("shouldn't be here\n");
      exit(-1);
    }
    return d2AdT2;
  }

  /*
   * ----------- Critical State Properties --------------------------
   */

  /*
   * ---------- Other Property Functions
   */
  double HMWSoln::AionicRadius(int k) const {
    return m_Aionic[k];
  }

  /*
   * ------------ Private and Restricted Functions ------------------
   */

  /**
   * Bail out of functions with an error exit if they are not
   * implemented.
   */
  doublereal HMWSoln::err(string msg) const {
    throw CanteraError("HMWSoln",
		       "Unfinished func called: " + msg );
    return 0.0;
  }


  /**
   * initLengths():
   *
   * This internal function adjusts the lengths of arrays based on
   * the number of species. This is done before these arrays are
   * populated with parameter values.
   */
  void HMWSoln::initLengths() {
    m_kk = nSpecies();
    MolalityVPSSTP::initThermo();
 
    /*
     * Resize lengths equal to the number of species in
     * the phase.
     */
    int leng = m_kk;
    m_electrolyteSpeciesType.resize(m_kk, cEST_polarNeutral);
    m_speciesSize.resize(leng);
    m_Aionic.resize(leng, 0.0);

    m_expg0_RT.resize(leng, 0.0);
    m_pe.resize(leng, 0.0);
    m_pp.resize(leng, 0.0);
    m_tmpV.resize(leng, 0.0);


    int maxCounterIJlen = 1 + (leng-1) * (leng-2) / 2;

    /*
     * Figure out the size of the temperature coefficient
     * arrays
     */
    int TCoeffLength = 1;
    if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
      TCoeffLength = 2;
    } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
      TCoeffLength = 5;
    }

    m_Beta0MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta0MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);
	
    m_Beta1MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_P.resize(maxCounterIJlen, 0.0);
    m_Beta1MX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Beta2MX_ij.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_L.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Beta2MX_ij_P.resize(maxCounterIJlen, 0.0);

    m_CphiMX_ij.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_L.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_LL.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_P.resize(maxCounterIJlen, 0.0);
    m_CphiMX_ij_coeff.resize(TCoeffLength, maxCounterIJlen, 0.0);

    m_Alpha1MX_ij.resize(maxCounterIJlen, 0.0);
    m_Theta_ij.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_L.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_LL.resize(maxCounterIJlen, 0.0);
    m_Theta_ij_P.resize(maxCounterIJlen, 0.0);

    m_Psi_ijk.resize(m_kk*m_kk*m_kk, 0.0);
    m_Psi_ijk_L.resize(m_kk*m_kk*m_kk, 0.0);
    m_Psi_ijk_LL.resize(m_kk*m_kk*m_kk, 0.0);
    m_Psi_ijk_P.resize(m_kk*m_kk*m_kk, 0.0);

    m_Lambda_ij.resize(leng, leng, 0.0);
    m_Lambda_ij_L.resize(leng, leng, 0.0);
    m_Lambda_ij_LL.resize(leng, leng, 0.0);
    m_Lambda_ij_P.resize(leng, leng, 0.0);

    m_lnActCoeffMolal.resize(leng, 0.0);
    m_dlnActCoeffMolaldT.resize(leng, 0.0);
    m_d2lnActCoeffMolaldT2.resize(leng, 0.0);
    m_dlnActCoeffMolaldP.resize(leng, 0.0);

    m_CounterIJ.resize(m_kk*m_kk, 0);

    m_gfunc_IJ.resize(maxCounterIJlen, 0.0);
    m_hfunc_IJ.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BprimeMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_BphiMX_IJ_P.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_L.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_Phi_IJ_P.resize(maxCounterIJlen, 0.0);
    m_Phiprime_IJ.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_L.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_PhiPhi_IJ_P.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_L.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_LL.resize(maxCounterIJlen, 0.0);
    m_CMX_IJ_P.resize(maxCounterIJlen, 0.0);

    m_gamma.resize(leng, 0.0);

    counterIJ_setup();
  }

  /**
   * Calcuate the natural log of the molality-based
   * activity coefficients.
   *
   */
  void HMWSoln::s_update_lnMolalityActCoeff() const {

    /*
     * Calculate the molalities. Currently, the molalities
     * may not be current with respect to the contents of the
     * State objects' data.
     */
    calcMolalities();
    /*
     * Calculate the stoichiometric ionic charge. This isn't used in the
     * Pitzer formulation.
     */
    m_IionicMolalityStoich = 0.0;
    for (int k = 0; k < m_kk; k++) {
      double z_k = m_speciesCharge[k];
      double zs_k1 =  m_speciesCharge_Stoich[k];
      if (z_k == zs_k1) {
	m_IionicMolalityStoich += m_molalities[k] * z_k * z_k;
      } else {
	double zs_k2 = z_k - zs_k1;
	m_IionicMolalityStoich
	  += m_molalities[k] * (zs_k1 * zs_k1 + zs_k2 * zs_k2);
      }
    }
    m_IionicMolalityStoich /= 2.0;                       
    if (m_IionicMolalityStoich > m_maxIionicStrength) {
      m_IionicMolalityStoich = m_maxIionicStrength;
    }

    /*
     * Update the temperature dependence of the pitzer coefficients
     * and their derivatives
     */
    s_updatePitzerCoeffWRTemp();

    /*
     * Now do the main calculation.
     */
    s_updatePitzerSublnMolalityActCoeff();
  }

  /*
   * Set up a counter variable for keeping track of symmetric binary
   * interactactions amongst the solute species.
   *
   * n = m_kk*i + j 
   * m_Counter[n] = counter
   */
  void HMWSoln::counterIJ_setup(void) const {
    int n, nc, i, j;
    m_CounterIJ.resize(m_kk * m_kk);
    int counter = 0;
    for (i = 0; i < m_kk; i++) {
      n = i;
      nc = m_kk * i;
      m_CounterIJ[n] = 0;
      m_CounterIJ[nc] = 0;
    }
    for (i = 1; i < (m_kk - 1); i++) {
      n = m_kk * i + i;
      m_CounterIJ[n] = 0;
      for (j = (i+1); j < m_kk; j++) {
	n = m_kk * j + i;
	nc = m_kk * i + j;
	counter++;
	m_CounterIJ[n] = counter;
	m_CounterIJ[nc] = counter;
      }
    }
  }

  /**
   * Calculates the Pitzer coefficients' dependence on the
   * temperature. It will also calculate the temperature
   * derivatives of the coefficients, as they are important
   * in the calculation of the latent heats and the
   * heat capacities of the mixtures.
   *
   * @param doDerivs If >= 1, then the routine will calculate
   *                 the first derivative. If >= 2, the 
   *                 routine will calculate the first and second
   *                 temperature derivative.
   *                 default = 2
   */
  void HMWSoln::s_updatePitzerCoeffWRTemp(int doDerivs) const {

    int i, j, n, counterIJ;
    const double *beta0MX_coeff;
    const double *beta1MX_coeff;
    const double *CphiMX_coeff;
    double T = temperature();
    double Tr = m_TempPitzerRef;
    double tinv = 0.0, tln = 0.0, tlin = 0.0, tquad = 0.0;
    if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
      tlin = T - Tr;
    } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
      tlin = T - Tr;
      tquad = T * T - Tr * Tr;
      tln = log(T/ Tr);
      tinv = 1.0/T - 1.0/Tr;
    }

    for (i = 1; i < (m_kk - 1); i++) {
      for (j = (i+1); j < m_kk; j++) {

	    
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	    
	beta0MX_coeff = m_Beta0MX_ij_coeff.ptrColumn(counterIJ);
	beta1MX_coeff = m_Beta1MX_ij_coeff.ptrColumn(counterIJ);
	CphiMX_coeff = m_CphiMX_ij_coeff.ptrColumn(counterIJ);

	switch (m_formPitzerTemp) {
	case PITZER_TEMP_CONSTANT:
	  break;
	case PITZER_TEMP_LINEAR:
	  m_Beta0MX_ij[counterIJ] = beta0MX_coeff[0] 
	    + beta0MX_coeff[1]*tlin;
	  m_Beta0MX_ij_L[counterIJ] = beta0MX_coeff[1];
	  m_Beta0MX_ij_LL[counterIJ] = 0.0;
	  m_Beta1MX_ij[counterIJ]   = beta1MX_coeff[0]
	    + beta1MX_coeff[1]*tlin;
	  m_Beta1MX_ij_L[counterIJ] = beta1MX_coeff[1];
	  m_Beta1MX_ij_LL[counterIJ] = 0.0;
	  m_CphiMX_ij [counterIJ]   = CphiMX_coeff[0]
	    + CphiMX_coeff[1]*tlin;
	  m_CphiMX_ij_L[counterIJ]  = CphiMX_coeff[1];
	  m_CphiMX_ij_LL[counterIJ] = 0.0;
	  break;

	case PITZER_TEMP_COMPLEX1:
	  m_Beta0MX_ij[counterIJ] = beta0MX_coeff[0] 
	    + beta0MX_coeff[1]*tlin
	    + beta0MX_coeff[2]*tquad
	    + beta0MX_coeff[3]*tinv
	    + beta0MX_coeff[4]*tln;
		
	  m_Beta1MX_ij[counterIJ] = beta1MX_coeff[0] 
	    + beta1MX_coeff[1]*tlin
	    + beta1MX_coeff[2]*tquad
	    + beta1MX_coeff[3]*tinv
	    + beta1MX_coeff[4]*tln;

	  m_CphiMX_ij[counterIJ] = CphiMX_coeff[0] 
	    + CphiMX_coeff[1]*tlin
	    + CphiMX_coeff[2]*tquad
	    + CphiMX_coeff[3]*tinv
	    + CphiMX_coeff[4]*tln;

	  m_Beta0MX_ij_L[counterIJ] =  beta0MX_coeff[1]
	    + beta0MX_coeff[2]*2.0*T
	    - beta0MX_coeff[3]/(T*T)
	    + beta0MX_coeff[4]/T;

	  m_Beta1MX_ij_L[counterIJ] =  beta1MX_coeff[1]
	    + beta1MX_coeff[2]*2.0*T
	    - beta1MX_coeff[3]/(T*T)
	    + beta1MX_coeff[4]/T;


	  m_CphiMX_ij_L[counterIJ] =  CphiMX_coeff[1]
	    + CphiMX_coeff[2]*2.0*T
	    - CphiMX_coeff[3]/(T*T)
	    + CphiMX_coeff[4]/T;
	
	  doDerivs = 2;
	  if (doDerivs > 1) {
	    m_Beta0MX_ij_LL[counterIJ] = 
	      + beta0MX_coeff[2]*2.0
	      + 2.0*beta0MX_coeff[3]/(T*T*T)
	      - beta0MX_coeff[4]/(T*T);
		  
	    m_Beta1MX_ij_LL[counterIJ] =
	      + beta1MX_coeff[2]*2.0
	      + 2.0*beta1MX_coeff[3]/(T*T*T)
	      - beta1MX_coeff[4]/(T*T);
		  
	    m_CphiMX_ij_LL[counterIJ] = 
	      + CphiMX_coeff[2]*2.0
	      + 2.0*CphiMX_coeff[3]/(T*T*T)
	      - CphiMX_coeff[4]/(T*T);
	  }

#ifdef DEBUG_HKM
	  /*
	   * Turn terms off for debugging
	   */
	  //m_Beta0MX_ij_L[counterIJ] = 0;
	  //m_Beta0MX_ij_LL[counterIJ] = 0;
	  //m_Beta1MX_ij_L[counterIJ] = 0;
	  //m_Beta1MX_ij_LL[counterIJ] = 0;
	  //m_CphiMX_ij_L[counterIJ] = 0;
	  //m_CphiMX_ij_LL[counterIJ] = 0;
#endif
	  break;
	}
	    
	   

      }
    }

  }
  /**
   * Calculate the Pitzer portion of the activity coefficients.
   *
   * This is the main routine in the whole module. It calculates the
   * molality based activity coefficients for the solutes, and
   * the activity of water.
   */
  void HMWSoln::
  s_updatePitzerSublnMolalityActCoeff() const {

    /*
     * HKM -> Assumption is made that the solvent is
     *        species 0.
     */
    if (m_indexSolvent != 0) {
      printf("Wrong index solvent value!\n");
      exit(-1);
    }

#ifdef DEBUG_MODE
    int printE = 0;
    if (temperature() == 323.15) {
      printE = 0;
    }
#endif
    double wateract;
    string sni,  snj, snk; 

    /*
     * This is the molality of the species in solution. 
     */
    const double *molality = DATA_PTR(m_molalities);
    /*
     * These are the charges of the species accessed from Constituents.h
     */
    const double *charge = DATA_PTR(m_speciesCharge);

    /*
     * These are data inputs about the Pitzer correlation. They come
     * from the input file for the Pitzer model.
     */
    const double *beta0MX =  DATA_PTR(m_Beta0MX_ij);
    const double *beta1MX =  DATA_PTR(m_Beta1MX_ij);
    const double *beta2MX =  DATA_PTR(m_Beta2MX_ij);
    const double *CphiMX  =  DATA_PTR(m_CphiMX_ij);
    const double *thetaij =  DATA_PTR(m_Theta_ij);
    const double *alphaMX =  DATA_PTR(m_Alpha1MX_ij);

    const double *psi_ijk =  DATA_PTR(m_Psi_ijk);
    //n = k + j * m_kk + i * m_kk * m_kk;


    double *gamma = DATA_PTR(m_gamma);
    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation, 
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double *g        =  DATA_PTR(m_gfunc_IJ);
    double *hfunc    =  DATA_PTR(m_hfunc_IJ);
    double *BMX      =  DATA_PTR(m_BMX_IJ);
    double *BprimeMX =  DATA_PTR(m_BprimeMX_IJ);
    double *BphiMX   =  DATA_PTR(m_BphiMX_IJ);
    double *Phi      =  DATA_PTR(m_Phi_IJ);
    double *Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double *Phiphi   =  DATA_PTR(m_PhiPhi_IJ);
    double *CMX      =  DATA_PTR(m_CMX_IJ);


    double x, g12rooti, gprime12rooti;
    double Aphi, F, zsqF;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, osmotic_coef, lnwateract;

    int z1, z2;
    int n, i, j, k, m, counterIJ,  counterIJ2;

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf("\n Debugging information from hmw_act \n");
    }
#endif
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();

    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
      //      ionic strength
      Is += charge[n] * charge[n] * molality[n];
      //      total molar charge
      molarcharge +=  fabs(charge[n]) * molality[n];
      molalitysum += molality[n];
    }
    Is *= 0.5;
    if (Is > m_maxIionicStrength) {
      Is = m_maxIionicStrength;
    }
    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 1: \n");
      printf(" ionic strenth      = %14.7le \n total molar "
	     "charge = %14.7le \n", Is, molarcharge);
    }
#endif

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the 
     * ionic strength (Is)
     */
    calc_lambdas(Is);

    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive 
     *                unlike charges up to 4
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 2: \n");
    }
#endif
    for (z1 = 1; z1 <=4; z1++) {
      for (z2 =1; z2 <=4; z2++) {
	calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n", 
		 z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
	}
#endif
      }
    }

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 3: \n");
      printf(" Species          Species            g(x) "
	     " hfunc(x)   \n");
    }
#endif
	
    /*
     *
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
      for (j = (i+1); j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * Only loop over oppositely charge species
	 */
	if (charge[i]*charge[j] < 0) {
	  /*
	   * x is a reduced function variable
	   */
	  x = sqrtIs * alphaMX[counterIJ];
	  if (x > 1.0E-100) {
	    g[counterIJ]     =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
	    hfunc[counterIJ] = -2.0*
	      (1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
	  }
	  else {
	    g[counterIJ]     = 0.0;
	    hfunc[counterIJ] = 0.0;
	  }
	} 
	else {
	  g[counterIJ]     = 0.0;
	  hfunc[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(), 
		 g[counterIJ], hfunc[counterIJ]);
	}
#endif
      }
    }

    /*
     * --------- SUBSECTION TO CALCULATE BMX, BprimeMX, BphiMX ----------
     * --------- Agrees with Pitzer, Eq. (49), (51), (55)
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 4: \n");
      printf(" Species          Species            BMX    "
	     "BprimeMX    BphiMX   \n");
    }
#endif
    x = 12.0 * sqrtIs; 
    if (x > 1.0E-100) {
      g12rooti      =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
      gprime12rooti = -2.0*(1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
    } else {
      g12rooti = 0.0;
      gprime12rooti = 0.0;
    }

    for (i = 1; i < m_kk - 1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];

#ifdef DEBUG_MODE
	if (printE) {
	  if (counterIJ == 2) {
	    printf("%s %s\n", speciesName(i).c_str(),
		   speciesName(j).c_str());
	    printf("beta0MX[%d] = %g\n", counterIJ, beta0MX[counterIJ]);
	    printf("beta1MX[%d] = %g\n", counterIJ, beta1MX[counterIJ]);
	  }
	}
#endif
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {	
	  BMX[counterIJ]  = beta0MX[counterIJ]
	    + beta1MX[counterIJ] * g[counterIJ]
	    + beta2MX[counterIJ] * g12rooti;
#ifdef DEBUG_MODE
	  if (m_debugCalc) {
	    printf("%d %g: %g %g %g\n",
		   counterIJ,  BMX[counterIJ], beta0MX[counterIJ],
		   beta1MX[counterIJ], g[counterIJ]);
	  }
#endif
	  if (Is > 1.0E-150) {
	    BprimeMX[counterIJ] = (beta1MX[counterIJ] * hfunc[counterIJ]/Is +
				   beta2MX[counterIJ] * gprime12rooti/Is);
	  } else {
	    BprimeMX[counterIJ] = 0.0;
	  }
	  BphiMX[counterIJ]   = BMX[counterIJ] + Is*BprimeMX[counterIJ];
	} 
	else {
	  BMX[counterIJ]      = 0.0;
	  BprimeMX[counterIJ] = 0.0;
	  BphiMX[counterIJ]   = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f %11.7f %11.7f \n", 
		 sni.c_str(), snj.c_str(), 
		 BMX[counterIJ], BprimeMX[counterIJ], BphiMX[counterIJ] );
	}
#endif
      }
    }	

    /*
     * --------- SUBSECTION TO CALCULATE CMX ----------
     * --------- Agrees with Pitzer, Eq. (53).
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 5: \n");
      printf(" Species          Species            CMX \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {
	  CMX[counterIJ] = CphiMX[counterIJ]/ 
	    (2.0* sqrt(fabs(charge[i]*charge[j])));
	} 
	else {
	  CMX[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (printE) {
	  if (counterIJ == 2) {
	    printf("%s %s\n", speciesName(i).c_str(), 
		   speciesName(j).c_str());
	    printf("CphiMX[%d] = %g\n", counterIJ, CphiMX[counterIJ]);
	  }
	}
#endif
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
		 CMX[counterIJ]);
	}
#endif
      }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------- Agrees with Pitzer, Eq. 72, 73, 74
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 6: \n");
      printf(" Species          Species            Phi_ij "
	     " Phiprime_ij  Phi^phi_ij \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] > 0) {
	  z1 = (int) fabs(charge[i]);
	  z2 = (int) fabs(charge[j]);
	  Phi[counterIJ] = thetaij[counterIJ] + etheta[z1][z2];
	  Phiprime[counterIJ] = etheta_prime[z1][z2];
	  Phiphi[counterIJ] = Phi[counterIJ] + Is * Phiprime[counterIJ];
	} 
	else {
	  Phi[counterIJ]      = 0.0;
	  Phiprime[counterIJ] = 0.0;
	  Phiphi[counterIJ]   = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %10.6f %10.6f %10.6f \n", 
		 sni.c_str(), snj.c_str(),
		 Phi[counterIJ], Phiprime[counterIJ], Phiphi[counterIJ] );
	}
#endif
      }
    }

    /*
     * ------------- SUBSECTION FOR CALCULATION OF F ----------------------
     * ------------ Agrees with Pitzer Eqn. (65) --------------------------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 7: \n");
    }
#endif
    // A_Debye_Huckel = 0.5092; (units = sqrt(kg/gmol))
    // A_Debye_Huckel = 0.5107; <- This value is used to match GWB data
    //                             ( A * ln(10) = 1.17593)
    // Aphi = A_Debye_Huckel * 2.30258509 / 3.0;
    Aphi = m_A_Debye / 3.0;
    F = -Aphi * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
		  + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
#ifdef DEBUG_MODE
    if (printE) {
      printf("Aphi = %20.13g\n", Aphi);
    }
#endif
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" initial value of F = %10.6f \n", F );		
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0) {
	  F = F + molality[i]*molality[j] * BprimeMX[counterIJ];
	}
	/*
	 * Both species have a non-zero charge, and they
	 * have the same sign
	 */
	if (charge[i]*charge[j] > 0) {
	  F = F + molality[i]*molality[j] * Phiprime[counterIJ];
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) printf(" F = %10.6f \n", F );
#endif
      }
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 8: \n");
    }
#endif

    for (i = 1; i < m_kk; i++) {

      /*
       * -------- SUBSECTION FOR CALCULATING THE ACTCOEFF FOR CATIONS -----
       * -------- -> equations agree with my notes, Eqn. (118).
       *          -> Equations agree with Pitzer, eqn.(63)
       */
      if (charge[i] > 0 ) {
	// species i is the cation (positive) to calc the actcoeff
	zsqF = charge[i]*charge[i]*F;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  if (charge[j] < 0.0) {
	    // sum over all anions
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX[counterIJ]+molarcharge*CMX[counterIJ]);
	    if (j < m_kk-1) {
	      /*
	       * This term is the ternary interaction involving the 
	       * non-duplicate sum over double anions, j, k, with
	       * respect to the cation, i.
	       */
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all anions
		if (charge[k] < 0.0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk[n];
		}
	      }
	    }
	  }

	     
	  if (charge[j] > 0.0) {
	    // sum over all cations
	    if (j != i) sum2 = sum2 + molality[j]*(2.0*Phi[counterIJ]);
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] < 0.0) {
		// two inner sums over anions

		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk[n];
		/*
		 * Find the counterIJ for the j,k interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + (fabs(charge[i])*
			       molality[j]*molality[k]*CMX[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * Handle neutral j species
	   */
	  if (charge[j] == 0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij(j,i);
	  }
	}
	/*
	 * Add all of the contributions up to yield the log of the
	 * solute activity coefficients (molality scale)
	 */
	m_lnActCoeffMolal[i] = zsqF + sum1 + sum2 + sum3 + sum4 + sum5;
	gamma[i] = exp(m_lnActCoeffMolal[i]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f \n", 
		 sni.c_str(), m_lnActCoeffMolal[i], gamma[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqF, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }

      /*
       * -------- SUBSECTION FOR CALCULATING THE ACTCOEFF FOR ANIONS ------
       * -------- -> equations agree with my notes, Eqn. (119).
       *          -> Equations agree with Pitzer, eqn.(64)
       */
      if (charge[i] < 0 ) {
	//          species i is an anion (negative)
	zsqF = charge[i]*charge[i]*F;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  /*
	   * For Anions, do the cation interactions.
	   */
	  if (charge[j] > 0) {
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX[counterIJ]+molarcharge*CMX[counterIJ]);
	    if (j < m_kk-1) {
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all cations
		if (charge[k] > 0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk[n];
		}
	      }
	    }
	  }

	  /*
	   * For Anions, do the other anion interactions.
	   */
	  if (charge[j] < 0.0) {
	    //  sum over all anions
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] > 0.0) {
		// two inner sums over cations
		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk[n];
		/*
		 * Find the counterIJ for the symmetric binary interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + 
		  (fabs(charge[i])*
		   molality[j]*molality[k]*CMX[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * for Anions, do the neutral species interaction
	   */
	  if (charge[j] == 0.0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij(j,i);
	  }
	}
	m_lnActCoeffMolal[i] = zsqF + sum1 + sum2 + sum3 + sum4 + sum5;
	gamma[i] = exp(m_lnActCoeffMolal[i]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f\n", 
		 sni.c_str(), m_lnActCoeffMolal[i], gamma[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqF, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }
      /*
       * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
       * ------ -> equations agree with my notes,
       *        -> Equations agree with Pitzer,
       */
      if (charge[i] == 0.0 ) {
	sum1 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  sum1 = sum1 + molality[j]*2.0*m_Lambda_ij(i,j);
	}
	m_lnActCoeffMolal[i] = sum1;
	gamma[i] = exp(m_lnActCoeffMolal[i]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f \n", 
		 sni.c_str(), m_lnActCoeffMolal[i], gamma[i]);
	}
#endif
      }

    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 9: \n");
    }
#endif
    /*
     * -------- SUBSECTION FOR CALCULATING THE OSMOTIC COEFF ---------
     * -------- -> equations agree with my notes, Eqn. (117).
     *          -> Equations agree with Pitzer, eqn.(62)
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    /*
     * term1 is the DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer 
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -Aphi * pow(Is,1.5) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
      /*
       * Loop Over Cations
       */
      if (charge[j] > 0.0) {
	for (k = 1; k < m_kk; k++){
	  if (charge[k] < 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	
	    sum1 = sum1 + molality[j]*molality[k]*
	      (BphiMX[counterIJ] + molarcharge*CMX[counterIJ]);
	  }
	}

	for (k = j+1; k < m_kk; k++) {
	  if (j == (m_kk-1)) {
	    // we should never reach this step
	    printf("logic error 1 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] > 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between 2 cations.
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	    sum2 = sum2 + molality[j]*molality[k]*Phiphi[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] < 0.0) {
		// species m is an anion
		n = m + k * m_kk + j * m_kk * m_kk;
		sum2 = sum2 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Anions
       */
      if (charge[j] < 0) {
	for (k = j+1; k < m_kk; k++) {
	  if (j == m_kk-1) {
	    // we should never reach this step
	    printf("logic error 2 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] < 0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between two anions
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];

	    sum3 = sum3 + molality[j]*molality[k]*Phiphi[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] > 0.0) {
		n = m + k * m_kk + j * m_kk * m_kk;
		sum3 = sum3 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Neutral Species
       */
      if (charge[j] == 0) {
	for (k = 1; k < m_kk; k++) {
	  if (charge[k] < 0.0) {
	    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_ij(j,k);
	  }
	  if (charge[k] > 0.0) {
	    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_ij(j,k);
	  }
	  if (charge[k] == 0.0) {
	    if (k > j) {
	      sum6 = sum6 + molality[j]*molality[k]*m_Lambda_ij(j,k);
	    } else if (k == j) {
	      sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_ij(j,k);
	    }
	  }
	}
      }
    }
    sum_m_phi_minus_1 = 2.0 * 
      (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6);
    /*
     * Calculate the osmotic coefficient from 
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
      osmotic_coef = 1.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
      osmotic_coef = 1.0;
    }
#ifdef DEBUG_MODE
    if (printE) {
      printf("OsmCoef - 1 = %20.13g\n", osmotic_coef - 1.0);
    }
#endif
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
	     "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
	     term1, sum1, sum2, sum3, sum4, sum5);
      printf("     sum_m_phi_minus_1=%10.6f        osmotic_coef=%10.6f\n", 
	     sum_m_phi_minus_1, osmotic_coef);
    }

    if (m_debugCalc) {
      printf(" Step 10: \n");
    }
#endif
    lnwateract = -(m_weightSolvent/1000.0) * molalitysum * osmotic_coef;
    wateract = exp(lnwateract);

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    m_lnActCoeffMolal[0] = lnwateract - log(xmolSolvent);
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Weight of Solvent = %16.7g\n", m_weightSolvent);
      printf(" molalitySum = %16.7g\n", molalitysum);
      printf(" ln_a_water=%10.6f a_water=%10.6f\n\n", 
	     lnwateract, wateract);
    }
#endif
  }

  /**
   * s_update_dlnMolalityActCoeff_dT()         (private, const )
   *
   *   Using internally stored values, this function calculates
   *   the temperature derivative of the logarithm of the
   *   activity coefficient for all species in the mechanism.
   *
   *   We assume that the activity coefficients are current.
   *
   *   solvent activity coefficient is on the molality
   *   scale. It's derivative is too.
   */
  void HMWSoln::s_update_dlnMolalityActCoeff_dT() const {

    for (int k = 0; k < m_kk; k++) {
      m_dlnActCoeffMolaldT[k] = 0.0;
    }
    s_Pitzer_dlnMolalityActCoeff_dT();
  }

  /*************************************************************************************/

  /**
   * Calculate the Pitzer portion of the temperature
   * derivative of the log activity coefficients.
   * This is an internal routine.
   *
   * It may be assumed that the 
   * Pitzer activity coefficient routine is called immediately
   * preceding the calling of this routine. Therefore, some
   * quantities do not need to be recalculated in this routine.
   *
   */
  void HMWSoln::s_Pitzer_dlnMolalityActCoeff_dT() const {

    /*
     * HKM -> Assumption is made that the solvent is
     *        species 0.
     */
#ifdef DEBUG_MODE
    m_debugCalc = 0;
#endif
    if (m_indexSolvent != 0) {
      printf("Wrong index solvent value!\n");
      exit(-1);
    }

    double d_wateract_dT;
    string sni, snj, snk; 

    const double *molality  =  DATA_PTR(m_molalities);
    const double *charge    =  DATA_PTR(m_speciesCharge);
    const double *beta0MX_L =  DATA_PTR(m_Beta0MX_ij_L);
    const double *beta1MX_L =  DATA_PTR(m_Beta1MX_ij_L);
    const double *beta2MX_L =  DATA_PTR(m_Beta2MX_ij_L);
    const double *CphiMX_L  =  DATA_PTR(m_CphiMX_ij_L);
    const double *thetaij_L =  DATA_PTR(m_Theta_ij_L);
    const double *alphaMX   =  DATA_PTR(m_Alpha1MX_ij);
    const double *psi_ijk_L =  DATA_PTR(m_Psi_ijk_L);
    double *gamma     =  DATA_PTR(m_gamma);
    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation, 
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double *g        =  DATA_PTR(m_gfunc_IJ);
    double *hfunc    =  DATA_PTR(m_hfunc_IJ);
    double *BMX_L    =  DATA_PTR(m_BMX_IJ_L);
    double *BprimeMX_L= DATA_PTR(m_BprimeMX_IJ_L);
    double *BphiMX_L =  DATA_PTR(m_BphiMX_IJ_L);
    double *Phi_L    =  DATA_PTR(m_Phi_IJ_L);
    double *Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double *Phiphi_L =  DATA_PTR(m_PhiPhi_IJ_L);
    double *CMX_L    =  DATA_PTR(m_CMX_IJ_L);

    double x, g12rooti, gprime12rooti;
    double Aphi, dFdT, zsqdFdT;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, d_osmotic_coef_dT, d_lnwateract_dT;

    int z1, z2;
    int n, i, j, k, m, counterIJ,  counterIJ2;

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf("\n Debugging information from "
	     "s_Pitzer_dlnMolalityActCoeff_dT()\n");
    }
#endif
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();

    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
      //      ionic strength
      Is += charge[n] * charge[n] * molality[n];
      //      total molar charge
      molarcharge +=  fabs(charge[n]) * molality[n];
      molalitysum += molality[n];
    }
    Is *= 0.5;
    if (Is > m_maxIionicStrength) {
      Is = m_maxIionicStrength;
    }
    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 1: \n");
      printf(" ionic strenth      = %14.7le \n total molar "
	     "charge = %14.7le \n", Is, molarcharge);
    }
#endif

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the 
     * ionic strength (Is)
     */
    calc_lambdas(Is);

    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive 
     *                unlike charges up to 4
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 2: \n");
    }
#endif
    for (z1 = 1; z1 <=4; z1++) {
      for (z2 =1; z2 <=4; z2++) {
	calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n", 
		 z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
	}
#endif
      }
    }

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 3: \n");
      printf(" Species          Species            g(x) "
	     " hfunc(x)   \n");
    }
#endif
	
    /*
     *
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
      for (j = (i+1); j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * Only loop over oppositely charge species
	 */
	if (charge[i]*charge[j] < 0) {
	  /*
	   * x is a reduced function variable
	   */
	  x = sqrtIs * alphaMX[counterIJ];
	  if (x > 1.0E-100) {
	    g[counterIJ]     =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
	    hfunc[counterIJ] = -2.0*
	      (1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
	  }
	  else {
	    g[counterIJ]     = 0.0;
	    hfunc[counterIJ] = 0.0;
	  }
	} 
	else {
	  g[counterIJ]     = 0.0;
	  hfunc[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(), 
		 g[counterIJ], hfunc[counterIJ]);
	}
#endif
      }
    }

    /*
     * ------- SUBSECTION TO CALCULATE BMX_L, BprimeMX_L, BphiMX_L ----------
     * ------- These are now temperature derivatives of the
     *         previously calculated quantities.
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 4: \n");
      printf(" Species          Species            BMX    "
	     "BprimeMX    BphiMX   \n");
    }
#endif
    x = 12.0 * sqrtIs; 
    if (x > 1.0E-100) {
      g12rooti      =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
      gprime12rooti = -2.0*(1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
    } else {
      g12rooti = 0.0;
      gprime12rooti = 0.0;
    }

    for (i = 1; i < m_kk - 1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {	
	  BMX_L[counterIJ]  = beta0MX_L[counterIJ]
	    + beta1MX_L[counterIJ] * g[counterIJ]
	    + beta2MX_L[counterIJ] * g12rooti;
#ifdef DEBUG_MODE
	  if (m_debugCalc) {
	    printf("%d %g: %g %g %g\n",
		   counterIJ,  BMX_L[counterIJ], beta0MX_L[counterIJ],
		   beta1MX_L[counterIJ], g[counterIJ]);
	  }
#endif
	  if (Is > 1.0E-150) {
	    BprimeMX_L[counterIJ] = (beta1MX_L[counterIJ] * hfunc[counterIJ]/Is +
				     beta2MX_L[counterIJ] * gprime12rooti/Is);
	  } else {
	    BprimeMX_L[counterIJ] = 0.0;
	  }
	  BphiMX_L[counterIJ] = BMX_L[counterIJ] + Is*BprimeMX_L[counterIJ];
	} 
	else {
	  BMX_L[counterIJ]      = 0.0;
	  BprimeMX_L[counterIJ] = 0.0;
	  BphiMX_L[counterIJ]     = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f %11.7f %11.7f \n", 
		 sni.c_str(), snj.c_str(), 
		 BMX_L[counterIJ], BprimeMX_L[counterIJ], BphiMX_L[counterIJ]);
	}
#endif
      }
    }	

    /*
     * --------- SUBSECTION TO CALCULATE CMX_L ----------
     * ---------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 5: \n");
      printf(" Species          Species            CMX \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {
	  CMX_L[counterIJ] = CphiMX_L[counterIJ]/ 
	    (2.0* sqrt(fabs(charge[i]*charge[j])));
	} 
	else {
	  CMX_L[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
		 CMX_L[counterIJ]);
	}
#endif
      }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 6: \n");
      printf(" Species          Species            Phi_ij "
	     " Phiprime_ij  Phi^phi_ij \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] > 0) {
	  z1 = (int) fabs(charge[i]);
	  z2 = (int) fabs(charge[j]);
	  //Phi[counterIJ] = thetaij_L[counterIJ] + etheta[z1][z2];
	  Phi_L[counterIJ] = thetaij_L[counterIJ];
	  //Phiprime[counterIJ] = etheta_prime[z1][z2];
	  Phiprime[counterIJ] = 0.0;
	  Phiphi_L[counterIJ] = Phi_L[counterIJ] + Is * Phiprime[counterIJ];
	} 
	else {
	  Phi_L[counterIJ]      = 0.0;
	  Phiprime[counterIJ] = 0.0;
	  Phiphi_L[counterIJ]   = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %10.6f %10.6f %10.6f \n", 
		 sni.c_str(), snj.c_str(),
		 Phi_L[counterIJ], Phiprime[counterIJ], Phiphi_L[counterIJ] );
	}
#endif
      }
    }

    /*
     * ----------- SUBSECTION FOR CALCULATION OF dFdT ---------------------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 7: \n");
    }
#endif
    // A_Debye_Huckel = 0.5092; (units = sqrt(kg/gmol))
    // A_Debye_Huckel = 0.5107; <- This value is used to match GWB data
    //                             ( A * ln(10) = 1.17593)
    // Aphi = A_Debye_Huckel * 2.30258509 / 3.0;
    Aphi = m_A_Debye / 3.0;

    double dA_DebyedT = dA_DebyedT_TP();
    double dAphidT = dA_DebyedT /3.0;
#ifdef DEBUG_HKM
    //dAphidT = 0.0;
#endif
    //F = -Aphi * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
    //      + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    //dAphidT = Al / (4.0 * GasConstant * T * T);
    dFdT = -dAphidT * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
			+ (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" initial value of dFdT = %10.6f \n", dFdT );		
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0) {
	  dFdT = dFdT + molality[i]*molality[j] * BprimeMX_L[counterIJ];
	}
	/*
	 * Both species have a non-zero charge, and they
	 * have the same sign, e.g., both positive or both negative.
	 */
	if (charge[i]*charge[j] > 0) {
	  dFdT = dFdT + molality[i]*molality[j] * Phiprime[counterIJ];
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) printf(" dFdT = %10.6f \n", dFdT);
#endif
      }
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 8: \n");
    }
#endif

    for (i = 1; i < m_kk; i++) {

      /*
       * -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
       * --
       */
      if (charge[i] > 0 ) {
	// species i is the cation (positive) to calc the actcoeff
	zsqdFdT = charge[i]*charge[i]*dFdT;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  if (charge[j] < 0.0) {
	    // sum over all anions
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX_L[counterIJ] + molarcharge*CMX_L[counterIJ]);
	    if (j < m_kk-1) {
	      /*
	       * This term is the ternary interaction involving the 
	       * non-duplicate sum over double anions, j, k, with
	       * respect to the cation, i.
	       */
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all anions
		if (charge[k] < 0.0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk_L[n];
		}
	      }
	    }
	  }

	     
	  if (charge[j] > 0.0) {
	    // sum over all cations
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi_L[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] < 0.0) {
		// two inner sums over anions

		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk_L[n];
		/*
		 * Find the counterIJ for the j,k interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + (fabs(charge[i])*
			       molality[j]*molality[k]*CMX_L[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * Handle neutral j species
	   */
	  if (charge[j] == 0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij_L(j,i);
	  }
	}
	/*
	 * Add all of the contributions up to yield the log of the
	 * solute activity coefficients (molality scale)
	 */
	m_dlnActCoeffMolaldT[i] =
	  zsqdFdT + sum1 + sum2 + sum3 + sum4 + sum5;
	gamma[i] = exp(m_dlnActCoeffMolaldT[i]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f \n", 
		 sni.c_str(), m_dlnActCoeffMolaldT[i], gamma[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqdFdT, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }

      /*
       * ------ SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR ANIONS ------
       *
       */
      if (charge[i] < 0 ) {
	//          species i is an anion (negative)
	zsqdFdT = charge[i]*charge[i]*dFdT;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  /*
	   * For Anions, do the cation interactions.
	   */
	  if (charge[j] > 0) {
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX_L[counterIJ] + molarcharge*CMX_L[counterIJ]);
	    if (j < m_kk-1) {
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all cations
		if (charge[k] > 0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk_L[n];
		}
	      }
	    }
	  }

	  /*
	   * For Anions, do the other anion interactions.
	   */
	  if (charge[j] < 0.0) {
	    //  sum over all anions
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi_L[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] > 0.0) {
		// two inner sums over cations
		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk_L[n];
		/*
		 * Find the counterIJ for the symmetric binary interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + 
		  (fabs(charge[i])*
		   molality[j]*molality[k]*CMX_L[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * for Anions, do the neutral species interaction
	   */
	  if (charge[j] == 0.0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij_L(j,i);
	  }
	}
	m_dlnActCoeffMolaldT[i] = 
	  zsqdFdT + sum1 + sum2 + sum3 + sum4 + sum5;
	gamma[i] = exp(m_dlnActCoeffMolaldT[i]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f\n", 
		 sni.c_str(), m_dlnActCoeffMolaldT[i], gamma[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqdFdT, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }
      /*
       * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
       * ------ -> equations agree with my notes,
       *        -> Equations agree with Pitzer,
       */
      if (charge[i] == 0.0 ) {
	sum1 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  sum1 = sum1 + molality[j]*2.0*m_Lambda_ij_L(i,j);
	}
	m_dlnActCoeffMolaldT[i] = sum1;
	gamma[i] = exp(m_dlnActCoeffMolaldT[i]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f gamma[i]=%10.6f \n", 
		 sni.c_str(), m_dlnActCoeffMolaldT[i], gamma[i]);
	}
#endif
      }

    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 9: \n");
    }
#endif
    /*
     * ------ SUBSECTION FOR CALCULATING THE d OSMOTIC COEFF dT ---------
     *
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    /*
     * term1 is the temperature derivative of the
     * DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer 
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -dAphidT * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
      /*
       * Loop Over Cations
       */
      if (charge[j] > 0.0) {
	for (k = 1; k < m_kk; k++){
	  if (charge[k] < 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	
	    sum1 = sum1 + molality[j]*molality[k]*
	      (BphiMX_L[counterIJ] + molarcharge*CMX_L[counterIJ]);
	  }
	}

	for (k = j+1; k < m_kk; k++) {
	  if (j == (m_kk-1)) {
	    // we should never reach this step
	    printf("logic error 1 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] > 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between 2 cations.
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	    sum2 = sum2 + molality[j]*molality[k]*Phiphi_L[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] < 0.0) {
		// species m is an anion
		n = m + k * m_kk + j * m_kk * m_kk;
		sum2 = sum2 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk_L[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Anions
       */
      if (charge[j] < 0) {
	for (k = j+1; k < m_kk; k++) {
	  if (j == m_kk-1) {
	    // we should never reach this step
	    printf("logic error 2 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] < 0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between two anions
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];

	    sum3 = sum3 + molality[j]*molality[k]*Phiphi_L[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] > 0.0) {
		n = m + k * m_kk + j * m_kk * m_kk;
		sum3 = sum3 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk_L[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Neutral Species
       */
      if (charge[j] == 0) {
	for (k = 1; k < m_kk; k++) {
	  if (charge[k] < 0.0) {
	    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_ij_L(j,k);
	  }
	  if (charge[k] > 0.0) {
	    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_ij_L(j,k);
	  }
	  if (charge[k] == 0.0) {
	    if (k > j) {
	      sum6 = sum6 + molality[j]*molality[k]*m_Lambda_ij_L(j,k);
	    } else if (k == j) {
	      sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_ij_L(j,k);
	    }
	  }
	}
      }
    }
    sum_m_phi_minus_1 = 2.0 * 
      (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6);
    /*
     * Calculate the osmotic coefficient from 
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
      d_osmotic_coef_dT = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
      d_osmotic_coef_dT = 0.0;
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
	     "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
	     term1, sum1, sum2, sum3, sum4, sum5);
      printf("     sum_m_phi_minus_1=%10.6f        d_osmotic_coef_dT =%10.6f\n", 
	     sum_m_phi_minus_1, d_osmotic_coef_dT);
    }

    if (m_debugCalc) {
      printf(" Step 10: \n");
    }
#endif
    d_lnwateract_dT = -(m_weightSolvent/1000.0) * molalitysum * d_osmotic_coef_dT;
    d_wateract_dT = exp(d_lnwateract_dT);

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    //double xmolSolvent = moleFraction(m_indexSolvent);
    m_dlnActCoeffMolaldT[0] = d_lnwateract_dT;
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" d_ln_a_water_dT = %10.6f d_a_water_dT=%10.6f\n\n", 
	     d_lnwateract_dT, d_wateract_dT); 
    }
#endif
  }

  /*************************************************************************************/


  /**
   * s_update_d2lnMolalityActCoeff_dT2()         (private, const )
   *
   *   Using internally stored values, this function calculates
   *   the temperature 2nd derivative of the logarithm of the
   *   activity coefficient for all species in the mechanism.
   *   This is an internal routine
   *
   *   We assume that the activity coefficients and first temperature
   *   derivatives of the activity coefficients  are current.
   *
   * It may be assumed that the 
   * Pitzer activity coefficient and first deriv routine are called immediately
   * preceding the calling of this routine. Therefore, some
   * quantities do not need to be recalculated in this routine.
   *
   *   solvent activity coefficient is on the molality
   *   scale. It's derivatives are too.
   */
  void HMWSoln::s_update_d2lnMolalityActCoeff_dT2() const {

    /*
     * HKM -> Assumption is made that the solvent is
     *        species 0.
     */
#ifdef DEBUG_MODE
    m_debugCalc = 0;
#endif
    if (m_indexSolvent != 0) {
      printf("Wrong index solvent value!\n");
      exit(-1);
    }

    double d2_wateract_dT2;
    string sni, snj, snk; 

    const double *molality  =  DATA_PTR(m_molalities);
    const double *charge    =  DATA_PTR(m_speciesCharge);
    const double *beta0MX_LL=  DATA_PTR(m_Beta0MX_ij_LL);
    const double *beta1MX_LL=  DATA_PTR(m_Beta1MX_ij_LL);
    const double *beta2MX_LL=  DATA_PTR(m_Beta2MX_ij_LL);
    const double *CphiMX_LL =  DATA_PTR(m_CphiMX_ij_LL);
    const double *thetaij_LL=  DATA_PTR(m_Theta_ij_LL);
    const double *alphaMX   =  DATA_PTR(m_Alpha1MX_ij);
    const double *psi_ijk_LL=  DATA_PTR(m_Psi_ijk_LL);

    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation, 
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double *g        =  DATA_PTR(m_gfunc_IJ);
    double *hfunc    =  DATA_PTR(m_hfunc_IJ);
    double *BMX_LL   =  DATA_PTR(m_BMX_IJ_LL);
    double *BprimeMX_LL=DATA_PTR(m_BprimeMX_IJ_LL);
    double *BphiMX_LL=  DATA_PTR(m_BphiMX_IJ_LL);
    double *Phi_LL   =  DATA_PTR(m_Phi_IJ_LL);
    double *Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double *Phiphi_LL=  DATA_PTR(m_PhiPhi_IJ_LL);
    double *CMX_LL   =  DATA_PTR(m_CMX_IJ_LL);


    double x, g12rooti, gprime12rooti;
    double d2FdT2, zsqd2FdT2;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, d2_osmotic_coef_dT2, d2_lnwateract_dT2;

    int z1, z2;
    int n, i, j, k, m, counterIJ,  counterIJ2;

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf("\n Debugging information from "
	     "s_Pitzer_d2lnMolalityActCoeff_dT2()\n");
    }
#endif
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();


    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
      //      ionic strength
      Is += charge[n] * charge[n] * molality[n];
      //      total molar charge
      molarcharge +=  fabs(charge[n]) * molality[n];
      molalitysum += molality[n];
    }
    Is *= 0.5;
    if (Is > m_maxIionicStrength) {
      Is = m_maxIionicStrength;
    }
    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 1: \n");
      printf(" ionic strenth      = %14.7le \n total molar "
	     "charge = %14.7le \n", Is, molarcharge);
    }
#endif

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the 
     * ionic strength (Is)
     */
    calc_lambdas(Is);

    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive 
     *                unlike charges up to 4
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 2: \n");
    }
#endif
    for (z1 = 1; z1 <=4; z1++) {
      for (z2 =1; z2 <=4; z2++) {
	calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n", 
		 z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
	}
#endif
      }
    }

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 3: \n");
      printf(" Species          Species            g(x) "
	     " hfunc(x)   \n");
    }
#endif
	
    /*
     *
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
      for (j = (i+1); j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * Only loop over oppositely charge species
	 */
	if (charge[i]*charge[j] < 0) {
	  /*
	   * x is a reduced function variable
	   */
	  x = sqrtIs * alphaMX[counterIJ];
	  if (x > 1.0E-100) {
	    g[counterIJ]     =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
	    hfunc[counterIJ] = -2.0*
	      (1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
	  }
	  else {
	    g[counterIJ]     = 0.0;
	    hfunc[counterIJ] = 0.0;
	  }
	} 
	else {
	  g[counterIJ]     = 0.0;
	  hfunc[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(), 
		 g[counterIJ], hfunc[counterIJ]);
	}
#endif
      }
    }
    /*
     * ------- SUBSECTION TO CALCULATE BMX_L, BprimeMX_LL, BphiMX_L ----------
     * ------- These are now temperature derivatives of the
     *         previously calculated quantities.
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 4: \n");
      printf(" Species          Species            BMX    "
	     "BprimeMX    BphiMX   \n");
    }
#endif
    x = 12.0 * sqrtIs; 
    if (x > 1.0E-100) {
      g12rooti      =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
      gprime12rooti = -2.0*(1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
    } else {
      g12rooti = 0.0;
      gprime12rooti = 0.0;
    }

    for (i = 1; i < m_kk - 1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {
	  BMX_LL[counterIJ]  = beta0MX_LL[counterIJ]
	    + beta1MX_LL[counterIJ] * g[counterIJ]
	    + beta2MX_LL[counterIJ] * g12rooti;
#ifdef DEBUG_MODE
	  if (m_debugCalc) {
	    printf("%d %g: %g %g %g\n",
		   counterIJ,  BMX_LL[counterIJ], beta0MX_LL[counterIJ],
		   beta1MX_LL[counterIJ], g[counterIJ]);
	  }
#endif
	  if (Is > 1.0E-150) {
	    BprimeMX_LL[counterIJ] = (beta1MX_LL[counterIJ] * hfunc[counterIJ]/Is +
				      beta2MX_LL[counterIJ] * gprime12rooti/Is);
	  } else {
	    BprimeMX_LL[counterIJ] = 0.0;
	  }
	  BphiMX_LL[counterIJ] = BMX_LL[counterIJ] + Is*BprimeMX_LL[counterIJ];
	} 
	else {
	  BMX_LL[counterIJ]      = 0.0;
	  BprimeMX_LL[counterIJ] = 0.0;
	  BphiMX_LL[counterIJ]     = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f %11.7f %11.7f \n", 
		 sni.c_str(), snj.c_str(), 
		 BMX_LL[counterIJ], BprimeMX_LL[counterIJ], BphiMX_LL[counterIJ]);
	}
#endif
      }
    }	

    /*
     * --------- SUBSECTION TO CALCULATE CMX_LL ----------
     * ---------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 5: \n");
      printf(" Species          Species            CMX \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {
	  CMX_LL[counterIJ] = CphiMX_LL[counterIJ]/ 
	    (2.0* sqrt(fabs(charge[i]*charge[j])));
	} else {
	  CMX_LL[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
		 CMX_LL[counterIJ]);
	}
#endif
      }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 6: \n");
      printf(" Species          Species            Phi_ij "
	     " Phiprime_ij  Phi^phi_ij \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] > 0) {
	  z1 = (int) fabs(charge[i]);
	  z2 = (int) fabs(charge[j]);
	  //Phi[counterIJ] = thetaij[counterIJ] + etheta[z1][z2];
	  //Phi_L[counterIJ] = thetaij_L[counterIJ];
	  Phi_LL[counterIJ] = thetaij_LL[counterIJ];
	  //Phiprime[counterIJ] = etheta_prime[z1][z2];
	  Phiprime[counterIJ] = 0.0;
	  //Phiphi[counterIJ] = Phi[counterIJ] + Is * Phiprime[counterIJ];
	  //Phiphi_L[counterIJ] = Phi_L[counterIJ] + Is * Phiprime[counterIJ];
	  Phiphi_LL[counterIJ] = Phi_LL[counterIJ];
	} 
	else {
	  Phi_LL[counterIJ]      = 0.0;
	  Phiprime[counterIJ] = 0.0;
	  Phiphi_LL[counterIJ]   = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  //printf(" %-16s %-16s %10.6f %10.6f %10.6f \n", 
	  //	     sni.c_str(), snj.c_str(),
	  //     Phi_L[counterIJ], Phiprime[counterIJ], Phiphi_L[counterIJ] );
	  printf(" %-16s %-16s %10.6f %10.6f %10.6f \n", 
		 sni.c_str(), snj.c_str(),
		 Phi_LL[counterIJ], Phiprime[counterIJ], Phiphi_LL[counterIJ] );
	}
#endif
      }
    }

    /*
     * ----------- SUBSECTION FOR CALCULATION OF d2FdT2 ---------------------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 7: \n");
    }
#endif
    // A_Debye_Huckel = 0.5092; (units = sqrt(kg/gmol))
    // A_Debye_Huckel = 0.5107; <- This value is used to match GWB data
    //                             ( A * ln(10) = 1.17593)
    // Aphi = A_Debye_Huckel * 2.30258509 / 3.0;
    // Aphi = m_A_Debye / 3.0;

    //double dA_DebyedT = dA_DebyedT_TP();
    //double dAphidT = dA_DebyedT /3.0;
    double d2AphidT2 = d2A_DebyedT2_TP() / 3.0;
#ifdef DEBUG_HKM
    //d2AphidT2 = 0.0;
#endif
    //F = -Aphi * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
    //      + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    //dAphidT = Al / (4.0 * GasConstant * T * T);
    //dFdT = -dAphidT * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
    //		    + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    d2FdT2 = -d2AphidT2 * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
			    + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
#ifdef DEBUG_MODE
    if (m_debugCalc) {	
      printf(" initial value of d2FdT2 = %10.6f \n", d2FdT2 );		
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0) {
	  d2FdT2 = d2FdT2 + molality[i]*molality[j] * BprimeMX_LL[counterIJ];
	}
	/*
	 * Both species have a non-zero charge, and they
	 * have the same sign, e.g., both positive or both negative.
	 */
	if (charge[i]*charge[j] > 0) {
	  d2FdT2 = d2FdT2 + molality[i]*molality[j] * Phiprime[counterIJ];
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) printf(" d2FdT2 = %10.6f \n", d2FdT2);
#endif
      }
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 8: \n");
    }
#endif

    for (i = 1; i < m_kk; i++) {

      /*
       * -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
       * --
       */
      if (charge[i] > 0 ) {
	// species i is the cation (positive) to calc the actcoeff
	zsqd2FdT2 = charge[i]*charge[i]*d2FdT2;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  if (charge[j] < 0.0) {
	    // sum over all anions
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX_LL[counterIJ] + molarcharge*CMX_LL[counterIJ]);
	    if (j < m_kk-1) {
	      /*
	       * This term is the ternary interaction involving the 
	       * non-duplicate sum over double anions, j, k, with
	       * respect to the cation, i.
	       */
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all anions
		if (charge[k] < 0.0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk_LL[n];
		}
	      }
	    }
	  }

	     
	  if (charge[j] > 0.0) {
	    // sum over all cations
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi_LL[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] < 0.0) {
		// two inner sums over anions

		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk_LL[n];
		/*
		 * Find the counterIJ for the j,k interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + (fabs(charge[i])*
			       molality[j]*molality[k]*CMX_LL[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * Handle neutral j species
	   */
	  if (charge[j] == 0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij_LL(j,i);
	  }
	}
	/*
	 * Add all of the contributions up to yield the log of the
	 * solute activity coefficients (molality scale)
	 */
	m_d2lnActCoeffMolaldT2[i] =
	  zsqd2FdT2 + sum1 + sum2 + sum3 + sum4 + sum5;
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s d2lngammadT2[i]=%10.6f \n", 
		 sni.c_str(), m_d2lnActCoeffMolaldT2[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqd2FdT2, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }


      /*
       * ------ SUBSECTION FOR CALCULATING THE d2ACTCOEFFdT2 FOR ANIONS ------
       *
       */
      if (charge[i] < 0 ) {
	//          species i is an anion (negative)
	zsqd2FdT2 = charge[i]*charge[i]*d2FdT2;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  /*
	   * For Anions, do the cation interactions.
	   */
	  if (charge[j] > 0) {
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX_LL[counterIJ] + molarcharge*CMX_LL[counterIJ]);
	    if (j < m_kk-1) {
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all cations
		if (charge[k] > 0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk_LL[n];
		}
	      }
	    }
	  }

	  /*
	   * For Anions, do the other anion interactions.
	   */
	  if (charge[j] < 0.0) {
	    //  sum over all anions
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi_LL[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] > 0.0) {
		// two inner sums over cations
		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk_LL[n];
		/*
		 * Find the counterIJ for the symmetric binary interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + 
		  (fabs(charge[i])*
		   molality[j]*molality[k]*CMX_LL[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * for Anions, do the neutral species interaction
	   */
	  if (charge[j] == 0.0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij_LL(j,i);
	  }
	}
	m_d2lnActCoeffMolaldT2[i] = 
	  zsqd2FdT2 + sum1 + sum2 + sum3 + sum4 + sum5;
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s d2lngammadT2[i]=%10.6f\n", 
		 sni.c_str(), m_d2lnActCoeffMolaldT2[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqd2FdT2, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }
      /*
       * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
       * ------ -> equations agree with my notes,
       *        -> Equations agree with Pitzer,
       */
      if (charge[i] == 0.0 ) {
	sum1 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  sum1 = sum1 + molality[j]*2.0*m_Lambda_ij_LL(i,j);
	}
	m_d2lnActCoeffMolaldT2[i] = sum1;
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s d2lngammadT2[i]=%10.6f \n", 
		 sni.c_str(), m_d2lnActCoeffMolaldT2[i]);
	}
#endif
      }

    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 9: \n");
    }
#endif

    /*
     * ------ SUBSECTION FOR CALCULATING THE d2 OSMOTIC COEFF dT2 ---------
     *
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    /*
     * term1 is the temperature derivative of the
     * DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer 
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -d2AphidT2 * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
      /*
       * Loop Over Cations
       */
      if (charge[j] > 0.0) {
	for (k = 1; k < m_kk; k++){
	  if (charge[k] < 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	
	    sum1 = sum1 + molality[j]*molality[k]*
	      (BphiMX_LL[counterIJ] + molarcharge*CMX_LL[counterIJ]);
	  }
	}

	for (k = j+1; k < m_kk; k++) {
	  if (j == (m_kk-1)) {
	    // we should never reach this step
	    printf("logic error 1 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] > 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between 2 cations.
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	    sum2 = sum2 + molality[j]*molality[k]*Phiphi_LL[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] < 0.0) {
		// species m is an anion
		n = m + k * m_kk + j * m_kk * m_kk;
		sum2 = sum2 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk_LL[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Anions
       */
      if (charge[j] < 0) {
	for (k = j+1; k < m_kk; k++) {
	  if (j == m_kk-1) {
	    // we should never reach this step
	    printf("logic error 2 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] < 0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between two anions
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];

	    sum3 = sum3 + molality[j]*molality[k]*Phiphi_LL[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] > 0.0) {
		n = m + k * m_kk + j * m_kk * m_kk;
		sum3 = sum3 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk_LL[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Neutral Species
       */
      if (charge[j] == 0) {
	for (k = 1; k < m_kk; k++) {
	  if (charge[k] < 0.0) {
	    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_ij_LL(j,k);
	  }
	  if (charge[k] > 0.0) {
	    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_ij_LL(j,k);
	  }
	  if (charge[k] == 0.0) {
	    if (k > j) {
	      sum6 = sum6 + molality[j]*molality[k]*m_Lambda_ij_LL(j,k);
	    } else if (k == j) {
	      sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_ij_LL(j,k);
	    }
	  }
	}
      }
    }
    sum_m_phi_minus_1 = 2.0 * 
      (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6);
    /*
     * Calculate the osmotic coefficient from 
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
      d2_osmotic_coef_dT2 = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
      d2_osmotic_coef_dT2 = 0.0;
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
	     "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
	     term1, sum1, sum2, sum3, sum4, sum5);
      printf("     sum_m_phi_minus_1=%10.6f        d2_osmotic_coef_dT2=%10.6f\n", 
	     sum_m_phi_minus_1, d2_osmotic_coef_dT2);
    }

    if (m_debugCalc) {
      printf(" Step 10: \n");
    }
#endif
    d2_lnwateract_dT2 = -(m_weightSolvent/1000.0) * molalitysum * d2_osmotic_coef_dT2;

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    m_d2lnActCoeffMolaldT2[0] = d2_lnwateract_dT2;

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      d2_wateract_dT2 = exp(d2_lnwateract_dT2);
      printf(" d2_ln_a_water_dT2 = %10.6f d2_a_water_dT2=%10.6f\n\n", 
	     d2_lnwateract_dT2, d2_wateract_dT2); 
    }
#endif
  }

  /********************************************************************************************/

  /**
   * s_Pitzer_dlnMolalityActCoeff_dP()         (private, const )
   *
   *   Using internally stored values, this function calculates
   *   the pressure derivative of the logarithm of the
   *   activity coefficient for all species in the mechanism.
   *
   *   We assume that the activity coefficients are current.
   *
   *   solvent activity coefficient is on the molality
   *   scale. It's derivative is too.
   */
  void HMWSoln::s_Pitzer_dlnMolalityActCoeff_dP() const {
 
    for (int k = 0; k < m_kk; k++) {
      m_dlnActCoeffMolaldP[k] = 0.0;
    }
    s_update_dlnMolalityActCoeff_dP();
  }

  /**
   * s_update_dlnMolalityActCoeff_dP()         (private, const )
   *
   *   Using internally stored values, this function calculates
   *   the pressure derivative of the logarithm of the
   *   activity coefficient for all species in the mechanism.
   *   This is an internal routine
   *
   *   We assume that the activity coefficients are current.
   *
   * It may be assumed that the 
   * Pitzer activity coefficient and first deriv routine are called immediately
   * preceding the calling of this routine. Therefore, some
   * quantities do not need to be recalculated in this routine.
   *
   *   solvent activity coefficient is on the molality
   *   scale. It's derivatives are too.
   */
  void HMWSoln::s_update_dlnMolalityActCoeff_dP() const {
   

    /*
     * HKM -> Assumption is made that the solvent is
     *        species 0.
     */
#ifdef DEBUG_MODE
    m_debugCalc = 0;
#endif
    if (m_indexSolvent != 0) {
      printf("Wrong index solvent value!\n");
      exit(-1);
    }

    double d_wateract_dP;
    string sni, snj, snk; 

    const double *molality  =  DATA_PTR(m_molalities);
    const double *charge    =  DATA_PTR(m_speciesCharge);
    const double *beta0MX_P =  DATA_PTR(m_Beta0MX_ij_P);
    const double *beta1MX_P =  DATA_PTR(m_Beta1MX_ij_P);
    const double *beta2MX_P =  DATA_PTR(m_Beta2MX_ij_P);
    const double *CphiMX_P  =  DATA_PTR(m_CphiMX_ij_P);
    const double *thetaij_P =  DATA_PTR(m_Theta_ij_P);
    const double *alphaMX   =  DATA_PTR(m_Alpha1MX_ij);
    const double *psi_ijk_P =  DATA_PTR(m_Psi_ijk_P);

    /*
     * Local variables defined by Coltrin
     */
    double etheta[5][5], etheta_prime[5][5], sqrtIs;
    /*
     * Molality based ionic strength of the solution
     */
    double Is = 0.0;
    /*
     * Molarcharge of the solution: In Pitzer's notation, 
     * this is his variable called "Z".
     */
    double molarcharge = 0.0;
    /*
     * molalitysum is the sum of the molalities over all solutes,
     * even those with zero charge.
     */
    double molalitysum = 0.0;

    double *g        =  DATA_PTR(m_gfunc_IJ);
    double *hfunc    =  DATA_PTR(m_hfunc_IJ);
    double *BMX_P    =  DATA_PTR(m_BMX_IJ_P);
    double *BprimeMX_P= DATA_PTR(m_BprimeMX_IJ_P);
    double *BphiMX_P =  DATA_PTR(m_BphiMX_IJ_P);
    double *Phi_P    =  DATA_PTR(m_Phi_IJ_P);
    double *Phiprime =  DATA_PTR(m_Phiprime_IJ);
    double *Phiphi_P =  DATA_PTR(m_PhiPhi_IJ_P);
    double *CMX_P    =  DATA_PTR(m_CMX_IJ_P);

    double x, g12rooti, gprime12rooti;
    double Aphi, dFdP, zsqdFdP;
    double sum1, sum2, sum3, sum4, sum5, term1;
    double sum_m_phi_minus_1, d_osmotic_coef_dP, d_lnwateract_dP;

    int z1, z2;
    int n, i, j, k, m, counterIJ,  counterIJ2;

    double currTemp = temperature();
    double currPres = pressure();

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf("\n Debugging information from "
	     "s_Pitzer_dlnMolalityActCoeff_dP()\n");
    }
#endif
    /*
     * Make sure the counter variables are setup
     */
    counterIJ_setup();

    /*
     * ---------- Calculate common sums over solutes ---------------------
     */
    for (n = 1; n < m_kk; n++) {
      //      ionic strength
      Is += charge[n] * charge[n] * molality[n];
      //      total molar charge
      molarcharge +=  fabs(charge[n]) * molality[n];
      molalitysum += molality[n];
    }
    Is *= 0.5;
    if (Is > m_maxIionicStrength) {
      Is = m_maxIionicStrength;
    }
    /*
     * Store the ionic molality in the object for reference.
     */
    m_IionicMolality = Is;
    sqrtIs = sqrt(Is);
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 1: \n");
      printf(" ionic strenth      = %14.7le \n total molar "
	     "charge = %14.7le \n", Is, molarcharge);
    }
#endif

    /*
     * The following call to calc_lambdas() calculates all 16 elements
     * of the elambda and elambda1 arrays, given the value of the 
     * ionic strength (Is)
     */
    calc_lambdas(Is);


    /*
     * ----- Step 2:  Find the coefficients E-theta and -------------------
     *                E-thetaprime for all combinations of positive 
     *                unlike charges up to 4
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 2: \n");
    }
#endif
    for (z1 = 1; z1 <=4; z1++) {
      for (z2 =1; z2 <=4; z2++) {
	calc_thetas(z1, z2, &etheta[z1][z2], &etheta_prime[z1][z2]);
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  printf(" z1=%3d z2=%3d E-theta(I) = %f, E-thetaprime(I) = %f\n", 
		 z1, z2, etheta[z1][z2], etheta_prime[z1][z2]);
	}
#endif
      }
    }

#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 3: \n");
      printf(" Species          Species            g(x) "
	     " hfunc(x)   \n");
    }
#endif
	
	
    /*
     *
     *  calculate g(x) and hfunc(x) for each cation-anion pair MX
     *   In the original literature, hfunc, was called gprime. However,
     *   it's not the derivative of g(x), so I renamed it.
     */
    for (i = 1; i < (m_kk - 1); i++) {
      for (j = (i+1); j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * Only loop over oppositely charge species
	 */
	if (charge[i]*charge[j] < 0) {
	  /*
	   * x is a reduced function variable
	   */
	  x = sqrtIs * alphaMX[counterIJ];
	  if (x > 1.0E-100) {
	    g[counterIJ]     =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
	    hfunc[counterIJ] = -2.0*
	      (1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
	  }
	  else {
	    g[counterIJ]     = 0.0;
	    hfunc[counterIJ] = 0.0;
	  }
	} 
	else {
	  g[counterIJ]     = 0.0;
	  hfunc[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %9.5f %9.5f \n", sni.c_str(), snj.c_str(), 
		 g[counterIJ], hfunc[counterIJ]);
	}
#endif
      }
    }


    /*
     * ------- SUBSECTION TO CALCULATE BMX_L, BprimeMX_L, BphiMX_L ----------
     * ------- These are now temperature derivatives of the
     *         previously calculated quantities.
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 4: \n");
      printf(" Species          Species            BMX    "
	     "BprimeMX    BphiMX   \n");
    }
#endif
    x = 12.0 * sqrtIs; 
    if (x > 1.0E-100) {
      g12rooti      =  2.0*(1.0-(1.0 + x) * exp(-x)) / (x*x);
      gprime12rooti = -2.0*(1.0-(1.0 + x + 0.5*x*x) * exp(-x)) / (x*x);
    } else {
      g12rooti = 0.0;
      gprime12rooti = 0.0;
    }

    for (i = 1; i < m_kk - 1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {	
	  BMX_P[counterIJ]  = beta0MX_P[counterIJ]
	    + beta1MX_P[counterIJ] * g[counterIJ]
	    + beta2MX_P[counterIJ] * g12rooti;
#ifdef DEBUG_MODE
	  if (m_debugCalc) {
	    printf("%d %g: %g %g %g\n",
		   counterIJ,  BMX_P[counterIJ], beta0MX_P[counterIJ],
		   beta1MX_P[counterIJ], g[counterIJ]);
	  }
#endif
	  if (Is > 1.0E-150) {
	    BprimeMX_P[counterIJ] = (beta1MX_P[counterIJ] * hfunc[counterIJ]/Is +
				     beta2MX_P[counterIJ] * gprime12rooti/Is);
	  } else {
	    BprimeMX_P[counterIJ] = 0.0;
	  }
	  BphiMX_P[counterIJ] = BMX_P[counterIJ] + Is*BprimeMX_P[counterIJ];
	} 
	else {
	  BMX_P[counterIJ]      = 0.0;
	  BprimeMX_P[counterIJ] = 0.0;
	  BphiMX_P[counterIJ]     = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f %11.7f %11.7f \n", 
		 sni.c_str(), snj.c_str(), 
		 BMX_P[counterIJ], BprimeMX_P[counterIJ], BphiMX_P[counterIJ]);
	}
#endif
      }
    }	


    /*
     * --------- SUBSECTION TO CALCULATE CMX_L ----------
     * ---------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 5: \n");
      printf(" Species          Species            CMX \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0.0) {
	  CMX_P[counterIJ] = CphiMX_P[counterIJ]/ 
	    (2.0* sqrt(fabs(charge[i]*charge[j])));
	} 
	else {
	  CMX_P[counterIJ] = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %11.7f \n", sni.c_str(), snj.c_str(),
		 CMX_P[counterIJ]);
	}
#endif
      }
    }

    /*
     * ------- SUBSECTION TO CALCULATE Phi, PhiPrime, and PhiPhi ----------
     * --------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 6: \n");
      printf(" Species          Species            Phi_ij "
	     " Phiprime_ij  Phi^phi_ij \n");
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] > 0) {
	  z1 = (int) fabs(charge[i]);
	  z2 = (int) fabs(charge[j]);
	  //Phi[counterIJ] = thetaij_L[counterIJ] + etheta[z1][z2];
	  Phi_P[counterIJ] = thetaij_P[counterIJ];
	  //Phiprime[counterIJ] = etheta_prime[z1][z2];
	  Phiprime[counterIJ] = 0.0;
	  Phiphi_P[counterIJ] = Phi_P[counterIJ] + Is * Phiprime[counterIJ];
	} 
	else {
	  Phi_P[counterIJ]      = 0.0;
	  Phiprime[counterIJ] = 0.0;
	  Phiphi_P[counterIJ]   = 0.0;
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  snj = speciesName(j);
	  printf(" %-16s %-16s %10.6f %10.6f %10.6f \n", 
		 sni.c_str(), snj.c_str(),
		 Phi_P[counterIJ], Phiprime[counterIJ], Phiphi_P[counterIJ] );
	}
#endif
      }
    }

    /*
     * ----------- SUBSECTION FOR CALCULATION OF dFdT ---------------------
     */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 7: \n");
    }
#endif
    // A_Debye_Huckel = 0.5092; (units = sqrt(kg/gmol))
    // A_Debye_Huckel = 0.5107; <- This value is used to match GWB data
    //                             ( A * ln(10) = 1.17593)
    // Aphi = A_Debye_Huckel * 2.30258509 / 3.0;
    Aphi = m_A_Debye / 3.0;

    double dA_DebyedP = dA_DebyedP_TP(currTemp, currPres);
    double dAphidP = dA_DebyedP /3.0;
#ifdef DEBUG_MODE
    //dAphidT = 0.0;
#endif
    //F = -Aphi * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
    //      + (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
    //dAphidT = Al / (4.0 * GasConstant * T * T);
    dFdP = -dAphidP * ( sqrt(Is) / (1.0 + 1.2*sqrt(Is)) 
			+ (2.0/1.2) * log(1.0+1.2*(sqrtIs)));
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" initial value of dFdP = %10.6f \n", dFdP );		
    }
#endif
    for (i = 1; i < m_kk-1; i++) {
      for (j = i+1; j < m_kk; j++) {
	/*
	 * Find the counterIJ for the symmetric binary interaction
	 */
	n = m_kk*i + j;
	counterIJ = m_CounterIJ[n];
	/*
	 * 	both species have a non-zero charge, and one is positive
	 *  and the other is negative
	 */
	if (charge[i]*charge[j] < 0) {
	  dFdP = dFdP + molality[i]*molality[j] * BprimeMX_P[counterIJ];
	}
	/*
	 * Both species have a non-zero charge, and they
	 * have the same sign, e.g., both positive or both negative.
	 */
	if (charge[i]*charge[j] > 0) {
	  dFdP = dFdP + molality[i]*molality[j] * Phiprime[counterIJ];
	}
#ifdef DEBUG_MODE
	if (m_debugCalc) printf(" dFdP = %10.6f \n", dFdP);
#endif
      }
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 8: \n");
    }
#endif


    for (i = 1; i < m_kk; i++) {

      /*
       * -------- SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR CATIONS -----
       * --
       */
      if (charge[i] > 0 ) {
	// species i is the cation (positive) to calc the actcoeff
	zsqdFdP = charge[i]*charge[i]*dFdP;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  if (charge[j] < 0.0) {
	    // sum over all anions
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX_P[counterIJ] + molarcharge*CMX_P[counterIJ]);
	    if (j < m_kk-1) {
	      /*
	       * This term is the ternary interaction involving the 
	       * non-duplicate sum over double anions, j, k, with
	       * respect to the cation, i.
	       */
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all anions
		if (charge[k] < 0.0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk_P[n];
		}
	      }
	    }
	  }


	     
	  if (charge[j] > 0.0) {
	    // sum over all cations
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi_P[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] < 0.0) {
		// two inner sums over anions

		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk_P[n];
		/*
		 * Find the counterIJ for the j,k interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + (fabs(charge[i])*
			       molality[j]*molality[k]*CMX_P[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * Handle neutral j species
	   */
	  if (charge[j] == 0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij_L(j,i);
	  }
	}

	/*
	 * Add all of the contributions up to yield the log of the
	 * solute activity coefficients (molality scale)
	 */
	m_dlnActCoeffMolaldP[i] =
	  zsqdFdP + sum1 + sum2 + sum3 + sum4 + sum5;

#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lngamma[i]=%10.6f \n", 
		 sni.c_str(), m_dlnActCoeffMolaldP[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqdFdP, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }

      /*
       * ------ SUBSECTION FOR CALCULATING THE dACTCOEFFdT FOR ANIONS ------
       *
       */
      if (charge[i] < 0 ) {
	//          species i is an anion (negative)
	zsqdFdP = charge[i]*charge[i]*dFdP;
	sum1 = 0.0;
	sum2 = 0.0;
	sum3 = 0.0;
	sum4 = 0.0;
	sum5 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  /*
	   * Find the counterIJ for the symmetric binary interaction
	   */
	  n = m_kk*i + j;
	  counterIJ = m_CounterIJ[n];

	  /*
	   * For Anions, do the cation interactions.
	   */
	  if (charge[j] > 0) {
	    sum1 = sum1 + molality[j]*
	      (2.0*BMX_P[counterIJ] + molarcharge*CMX_P[counterIJ]);
	    if (j < m_kk-1) {
	      for (k = j+1; k < m_kk; k++) {
		// an inner sum over all cations
		if (charge[k] > 0) {
		  n = k + j * m_kk + i * m_kk * m_kk;
		  sum3 = sum3 + molality[j]*molality[k]*psi_ijk_P[n];
		}
	      }
	    }
	  }

	  /*
	   * For Anions, do the other anion interactions.
	   */
	  if (charge[j] < 0.0) {
	    //  sum over all anions
	    if (j != i) {
	      sum2 = sum2 + molality[j]*(2.0*Phi_P[counterIJ]);
	    }
	    for (k = 1; k < m_kk; k++) {
	      if (charge[k] > 0.0) {
		// two inner sums over cations
		n = k + j * m_kk + i * m_kk * m_kk;
		sum2 = sum2 + molality[j]*molality[k]*psi_ijk_P[n];
		/*
		 * Find the counterIJ for the symmetric binary interaction
		 */
		n = m_kk*j + k;
		counterIJ2 = m_CounterIJ[n];
		sum4 = sum4 + 
		  (fabs(charge[i])*
		   molality[j]*molality[k]*CMX_P[counterIJ2]);
	      }
	    }
	  }

	  /*
	   * for Anions, do the neutral species interaction
	   */
	  if (charge[j] == 0.0) {
	    sum5 = sum5 + molality[j]*2.0*m_Lambda_ij_L(j,i);
	  }
	}
	m_dlnActCoeffMolaldP[i] = 
	  zsqdFdP + sum1 + sum2 + sum3 + sum4 + sum5;
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s lndactcoeffmolaldP[i]=%10.6f \n", 
		 sni.c_str(), m_dlnActCoeffMolaldP[i]);
	  printf("                   %12g %12g %12g %12g %12g %12g\n",
		 zsqdFdP, sum1, sum2, sum3, sum4, sum5);
	}
#endif
      }


      /*
       * ------ SUBSECTION FOR CALCULATING NEUTRAL SOLUTE ACT COEFF -------
       * ------ -> equations agree with my notes,
       *        -> Equations agree with Pitzer,
       */
      if (charge[i] == 0.0 ) {
	sum1 = 0.0;
	for (j = 1; j < m_kk; j++) {
	  sum1 = sum1 + molality[j]*2.0*m_Lambda_ij_L(i,j);
	}
	m_dlnActCoeffMolaldP[i] = sum1;
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  sni = speciesName(i);
	  printf(" %-16s dlnActCoeffMolaldP[i]=%10.6f \n", 
		 sni.c_str(), m_dlnActCoeffMolaldP[i]);
	}
#endif
      }

    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Step 9: \n");
    }
#endif

    /*
     * ------ SUBSECTION FOR CALCULATING THE d OSMOTIC COEFF dT ---------
     *
     */
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;
    sum4 = 0.0;
    sum5 = 0.0;
    double sum6 = 0.0;
    /*
     * term1 is the temperature derivative of the
     * DH term in the osmotic coefficient expression
     * b = 1.2 sqrt(kg/gmol) <- arbitrarily set in all Pitzer 
     *                          implementations.
     * Is = Ionic strength on the molality scale (units of (gmol/kg))
     * Aphi = A_Debye / 3   (units of sqrt(kg/gmol))
     */
    term1 = -dAphidP * Is * sqrt(Is) / (1.0 + 1.2 * sqrt(Is));

    for (j = 1; j < m_kk; j++) {
      /*
       * Loop Over Cations
       */
      if (charge[j] > 0.0) {
	for (k = 1; k < m_kk; k++){
	  if (charge[k] < 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	
	    sum1 = sum1 + molality[j]*molality[k]*
	      (BphiMX_P[counterIJ] + molarcharge*CMX_P[counterIJ]);
	  }
	}

	for (k = j+1; k < m_kk; k++) {
	  if (j == (m_kk-1)) {
	    // we should never reach this step
	    printf("logic error 1 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] > 0.0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between 2 cations.
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];
	    sum2 = sum2 + molality[j]*molality[k]*Phiphi_P[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] < 0.0) {
		// species m is an anion
		n = m + k * m_kk + j * m_kk * m_kk;
		sum2 = sum2 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk_P[n];
	      }
	    }
	  }
	}
      }
	 
	  
      /*
       * Loop Over Anions
       */
      if (charge[j] < 0) {
	for (k = j+1; k < m_kk; k++) {
	  if (j == m_kk-1) {
	    // we should never reach this step
	    printf("logic error 2 in Step 9 of hmw_act");
	    exit(1);
	  }
	  if (charge[k] < 0) {
	    /*
	     * Find the counterIJ for the symmetric j,k binary interaction
	     * between two anions
	     */
	    n = m_kk*j + k;
	    counterIJ = m_CounterIJ[n];

	    sum3 = sum3 + molality[j]*molality[k]*Phiphi_P[counterIJ];
	    for (m = 1; m < m_kk; m++) {
	      if (charge[m] > 0.0) {
		n = m + k * m_kk + j * m_kk * m_kk;
		sum3 = sum3 + 
		  molality[j]*molality[k]*molality[m]*psi_ijk_P[n];
	      }
	    }
	  }
	}
      }
	  
      /*
       * Loop Over Neutral Species
       */
      if (charge[j] == 0) {
	for (k = 1; k < m_kk; k++) {
	  if (charge[k] < 0.0) {
	    sum4 = sum4 + molality[j]*molality[k]*m_Lambda_ij_P(j,k);
	  }
	  if (charge[k] > 0.0) {
	    sum5 = sum5 + molality[j]*molality[k]*m_Lambda_ij_P(j,k);
	  }
	  if (charge[k] == 0.0) {
	    if (k > j) {
	      sum6 = sum6 + molality[j]*molality[k]*m_Lambda_ij_P(j,k);
	    } else if (k == j) {
	      sum6 = sum6 + 0.5 * molality[j]*molality[k]*m_Lambda_ij_P(j,k);
	    }
	  }
	}
      }
    }
    sum_m_phi_minus_1 = 2.0 * 
      (term1 + sum1 + sum2 + sum3 + sum4 + sum5 + sum6);


    /*
     * Calculate the osmotic coefficient from 
     *       osmotic_coeff = 1 + dGex/d(M0noRT) / sum(molality_i)
     */
    if (molalitysum > 1.0E-150) {
      d_osmotic_coef_dP = 0.0 + (sum_m_phi_minus_1 / molalitysum);
    } else {
      d_osmotic_coef_dP = 0.0;
    }
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" term1=%10.6f sum1=%10.6f sum2=%10.6f "
	     "sum3=%10.6f sum4=%10.6f sum5=%10.6f\n",
	     term1, sum1, sum2, sum3, sum4, sum5);
      printf("     sum_m_phi_minus_1=%10.6f        d_osmotic_coef_dP =%10.6f\n", 
	     sum_m_phi_minus_1, d_osmotic_coef_dP);
    }

    if (m_debugCalc) {
      printf(" Step 10: \n");
    }
#endif
    d_lnwateract_dP = -(m_weightSolvent/1000.0) * molalitysum * d_osmotic_coef_dP;
    d_wateract_dP = exp(d_lnwateract_dP);

    /*
     * In Cantera, we define the activity coefficient of the solvent as
     *
     *     act_0 = actcoeff_0 * Xmol_0
     *
     * We have just computed act_0. However, this routine returns
     *  ln(actcoeff[]). Therefore, we must calculate ln(actcoeff_0).
     */
    //double xmolSolvent = moleFraction(m_indexSolvent);
    m_dlnActCoeffMolaldP[0] = d_lnwateract_dP;
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" d_ln_a_water_dP = %10.6f d_a_water_dP=%10.6f\n\n", 
	     d_lnwateract_dP, d_wateract_dP); 
    }
#endif



  }

  /***********************************************************************************************/

  /*
   * Calculate the lambda interactions. 
   *
   * Calculate E-lambda terms for charge combinations of like sign,
   *   using method of Pitzer (1975).
   *
   *  This code snipet is included from Bethke, Appendix 2.
   */
  void HMWSoln::calc_lambdas(double is) const {
    double aphi, dj, jfunc, jprime, t, x, zprod;
    int i, ij, j;
    /*
     * Coefficients c1-c4 are used to approximate 
     * the integral function "J";
     * aphi is the Debye-Huckel constant at 25 C
     */

    double c1 = 4.581, c2 = 0.7237, c3 = 0.0120, c4 = 0.528;

    aphi = 0.392;   /* Value at 25 C */
#ifdef DEBUG_MODE
    if (m_debugCalc) {
      printf(" Is = %g\n", is);
    }
#endif
    if (is < 1.0E-150) {
      for (i = 0; i < 17; i++) {
	elambda[i] = 0.0;
	elambda1[i] = 0.0;
      }
      return;
    }
    /*
     * Calculate E-lambda terms for charge combinations of like sign,
     * using method of Pitzer (1975). Charges up to 4 are calculated.
     */

    for (i=1; i<=4; i++) {
      for (j=i; j<=4; j++) {
	ij = i*j;
	/*
	 * calculate the product of the charges
	 */
	zprod = (double)ij;
	/*
	 * calculate Xmn (A1) from Harvie, Weare (1980).
	 */
	x = 6.0* zprod * aphi * sqrt(is);                      /* eqn 23 */
	    
	jfunc = x / (4.0 + c1*pow(x,-c2)*exp(-c3*pow(x,c4)));  /* eqn 47 */

	t = c3 * c4 * pow(x,c4);
	dj = c1* pow(x,(-c2-1.0)) * (c2+t) * exp(-c3*pow(x,c4));
	jprime = (jfunc/x)*(1.0 + jfunc*dj);

	elambda[ij] = zprod*jfunc / (4.0*is);                  /* eqn 14 */
	elambda1[ij] = (3.0*zprod*zprod*aphi*jprime/(4.0*sqrt(is))
			- elambda[ij])/is;
#ifdef DEBUG_MODE
	if (m_debugCalc) {
	  printf(" ij = %d, elambda = %g, elambda1 = %g\n",
		 ij, elambda[ij], elambda1[ij]);
	}
#endif
      }
    }
  }

  /*
   * Calculate the etheta interaction. 
   * This interaction accounts for the mixing effects of like-signed
   * ions with different charges. There is fairly extensive literature
   * on this effect. See the notes.
   * This interaction will be nonzero for species with the same charge.
   *
   *  This code snipet is included from Bethke, Appendix 2.
   */
  void HMWSoln::calc_thetas(int z1, int z2,
			    double *etheta, double *etheta_prime) const {
    int i, j;
    double f1, f2;
	
    /*
     * Calculate E-theta(i) and E-theta'(I) using method of
     * Pitzer (1987) 
     */
    i = abs(z1);
    j = abs(z2);

#ifdef DEBUG_MODE
    if (i > 4 || j > 4) {
      printf("we shouldn't be here\n");
      exit(-1);
    }
#endif

    if ((i == 0) || (j == 0)) {
      printf("ERROR calc_thetas called with one species being neutral\n");
      exit(-1);
    }

    /*
     *  Check to see if the charges are of opposite sign. If they are of 
     *  opposite sign then their etheta interaction is zero.
     */
    if (z1*z2 < 0) {
      *etheta = 0.0;
      *etheta_prime = 0.0;
    }
    /*
     * Actually calculate the interaction.
     */
    else {
      f1 = (double)i / (2.0  * j);
      f2 = (double)j / (2.0  * i);
      *etheta = elambda[i*j] - f1*elambda[j*j] - f2*elambda[i*i];
      *etheta_prime = elambda1[i*j] - f1*elambda1[j*j] - f2*elambda1[i*i];
    }
  }
    
  /**
   * This routine prints out the input pitzer coefficients for the 
   * current mechanism
   */
  void HMWSoln::printCoeffs() const {
    int i, j, k;
    string sni, snj;
    calcMolalities();
    const double *charge = DATA_PTR(m_speciesCharge);
    double *molality = DATA_PTR(m_molalities);
    double *moleF = DATA_PTR(m_tmpV);
    /*
     * Update the coefficients wrt Temperature
     * Calculate the derivatives as well
     */
    s_updatePitzerCoeffWRTemp(2);
    getMoleFractions(moleF);

    printf("Index  Name                  MoleF      Molality      Charge\n");
    for (k = 0; k < m_kk; k++) {
      sni = speciesName(k);
      printf("%2d     %-16s %14.7le %14.7le %5.1f \n",
	     k, sni.c_str(), moleF[k], molality[k], charge[k]);
    }

    printf("\n Species          Species            beta0MX  "
	   "beta1MX   beta2MX   CphiMX    alphaMX thetaij    \n");
    for (i = 1; i < m_kk - 1; i++) {
      sni = speciesName(i);
      for (j = i+1; j < m_kk; j++) {
	snj = speciesName(j);
	int n  = i * m_kk + j;
	int ct = m_CounterIJ[n];
	printf(" %-16s %-16s %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n",
	       sni.c_str(), snj.c_str(),
	       m_Beta0MX_ij[ct], m_Beta1MX_ij[ct],
	       m_Beta2MX_ij[ct], m_CphiMX_ij[ct],
	       m_Alpha1MX_ij[ct], m_Theta_ij[ct] );


      }
    }

    printf("\n Species          Species          Species       "
	   "psi   \n");
    for (i = 1; i < m_kk; i++) {
      sni = speciesName(i);
      for (j = 1; j < m_kk; j++) {
	snj = speciesName(j);
	for (k = 1; k < m_kk; k++) {
	  string snk = speciesName(k);
	  int n = k + j * m_kk + i * m_kk * m_kk;
	  if (m_Psi_ijk[n] != 0.0) {
	    printf(" %-16s %-16s %-16s %9.5f \n",
		   sni.c_str(), snj.c_str(), 
		   snk.c_str(), m_Psi_ijk[n]);
	  }
	}
      }
    }
  }

 
  /*****************************************************************************/
}
/*****************************************************************************/
