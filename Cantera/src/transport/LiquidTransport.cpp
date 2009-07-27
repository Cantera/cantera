/**
 *  @file LiquidTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
/* 
 * $Revision: 1.10 $
 * $Date: 2009/03/24 20:44:30 $
 */

#include "ThermoPhase.h"
#include "LiquidTransport.h"

#include "utilities.h"
#include "LiquidTransportParams.h"
#include "TransportFactory.h"

#include "ctlapack.h"

#include <iostream>
using namespace std;

/** 
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#define MIN_X 1.e-14


namespace Cantera {

  //////////////////// class LiquidTransport methods //////////////


  LiquidTransport::LiquidTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim),
    m_nsp(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_logt(0.0),
    m_sqrt_t(-1.0),
    m_t14(-1.0),
    m_t32(-1.0),
    m_sqrt_kbt(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_visc_conc_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false),
    m_mode(-1000),
    m_debug(false)
  {
  }


  LiquidTransport::LiquidTransport(const LiquidTransport &right) :
    Transport(),
    m_nsp(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_logt(0.0),
    m_sqrt_t(-1.0),
    m_t14(-1.0),
    m_t32(-1.0),
    m_sqrt_kbt(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_visc_conc_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false),
    m_mode(-1000),
    m_debug(false)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = right;
  }

  LiquidTransport& LiquidTransport::operator=(const LiquidTransport& right) {
    if (&right != this) {
      return *this; 
    }
    Transport::operator=(right);
    m_nsp                                 = right.m_nsp;
    m_tmin                                = right.m_tmin;
    m_tmax                                = right.m_tmax;
    m_mw                                  = right.m_mw;
    m_poly                                = right.m_poly;
    viscCoeffsVector_                     = right.viscCoeffsVector_;
    m_condcoeffs                          = right.m_condcoeffs;
    m_diffcoeffs                          = right.m_diffcoeffs;
    m_Grad_X                              = right.m_Grad_X;
    m_Grad_T                              = right.m_Grad_T;
    m_Grad_V                              = right.m_Grad_V;
    m_ck_Grad_mu                          = right.m_ck_Grad_mu;
    m_bdiff                               = right.m_bdiff;
    viscSpecies_                          = right.viscSpecies_;
    m_sqvisc                              = right.m_sqvisc;
    m_cond                                = right.m_cond;
    m_polytempvec                         = right.m_polytempvec;
    m_iStateMF = -1;
    m_molefracs                           = right.m_molefracs;
    m_concentrations                      = right.m_concentrations;
    m_chargeSpecies                       = right.m_chargeSpecies;
    m_DiffCoeff_StefMax                   = right.m_DiffCoeff_StefMax;
    viscosityModel_                       = right.viscosityModel_;
    m_phi                                 = right.m_phi;
    m_wratjk                              = right.m_wratjk;
    m_wratkj1                             = right.m_wratkj1;
    m_B                                   = right.m_B;
    m_A                                   = right.m_A;
    m_eps                                 = right.m_eps;
    m_alpha                               = right.m_alpha;
    m_temp                                = right.m_temp;
    m_logt                                = right.m_logt;
    m_sqrt_t                              = right.m_sqrt_t;
    m_t14                                 = right.m_t14;
    m_t32                                 = right.m_t32;
    m_sqrt_kbt                            = right.m_sqrt_kbt;
    m_press                               = right.m_press;
    m_flux                                = right.m_flux;
    m_lambda                              = right.m_lambda;
    m_viscmix                             = right.m_viscmix;
    m_spwork                              = right.m_spwork;
    m_visc_mix_ok    = false;
    m_visc_temp_ok   = false;
    m_visc_conc_ok   = false;
    m_diff_mix_ok    = false;
    m_diff_temp_ok   = false;
    m_cond_temp_ok   = false;
    m_cond_mix_ok    = false;
    m_mode                                = right.m_mode;
    m_diam                                = right.m_diam;
    m_debug                               = right.m_debug;
    m_nDim                                = right.m_nDim;

    return *this; 
  }


  Transport *LiquidTransport::duplMyselfAsTransport() const {
    LiquidTransport* tr = new LiquidTransport(*this);
    return (dynamic_cast<Transport *>(tr));
  }

  // Initialize the object
  /*
   *  This is where we dimension everything.
   */
  bool LiquidTransport::initLiquid(LiquidTransportParams& tr) {

    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(), 
	 m_thermo->molecularWeights().end(), m_mw.begin());

    // copy polynomials and parameters into local storage
    m_poly       = tr.poly;
    viscCoeffsVector_ = tr.viscCoeffsVector_;
    m_condcoeffs = tr.condcoeffs;
    m_diffcoeffs = tr.diffcoeffs;

    m_mode       = tr.mode;
    m_diam       = tr.diam;
    m_eps        = tr.eps;
    m_alpha      = tr.alpha;

    m_phi.resize(m_nsp, m_nsp, 0.0);


    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    int j, k;
    for (j = 0; j < m_nsp; j++) 
      for (k = j; k < m_nsp; k++) {
	m_wratjk(j,k) = sqrt(m_mw[j]/m_mw[k]);
	m_wratjk(k,j) = sqrt(m_wratjk(j,k));
	m_wratkj1(j,k) = sqrt(1.0 + m_mw[k]/m_mw[j]);
      }
    
    m_polytempvec.resize(5);
    viscSpecies_.resize(m_nsp);
    m_sqvisc.resize(m_nsp);
    m_cond.resize(m_nsp);
    m_bdiff.resize(m_nsp, m_nsp);

    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);
    m_ck_Grad_mu.resize(m_nDim * m_nsp, 0.0);


    // set all flags to false
    m_visc_mix_ok   = false;
    m_visc_temp_ok  = false;
    m_visc_conc_ok  = false;

    m_cond_temp_ok = false;
    m_cond_mix_ok  = false;
    m_diff_temp_ok   = false;
    m_diff_mix_ok  = false;

    return true;
  }



  /******************  viscosity ******************************/

  /*
   * The viscosity is computed using the Wilke mixture rule.
   * \f[
   * \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
   * \f]
   * Here \f$ \mu_k \f$ is the viscosity of pure species \e k,
   * and 
   * \f[
   * \Phi_{k,j} = \frac{\left[1 
   * + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
   * {\sqrt{8}\sqrt{1 + M_k/M_j}}
   * \f] 
   * @see updateViscosity_T();
   */ 
  doublereal LiquidTransport::viscosity() {
        
    update_temp();
    update_conc();

    if (m_visc_mix_ok) return m_viscmix;
  
    // update viscSpecies_[] and m_phi[] if necessary
    if (!m_visc_temp_ok) {
      updateViscosity_temp();
    }

    if (!m_visc_conc_ok) {
      updateViscosities_conc();
    }

    if (viscosityModel_ == LVISC_CONSTANT) {
      return m_viscmix;
    } else if (viscosityModel_ == LVISC_MIXTUREAVG) {
      m_viscmix = dot_product(viscSpecies_, m_molefracs);
    } else if (viscosityModel_ == LVISC_WILKES) {
      multiply(m_phi, DATA_PTR(m_molefracs), DATA_PTR(m_spwork));
      m_viscmix = 0.0;
      for (int k = 0; k < m_nsp; k++) {
	m_viscmix += m_molefracs[k] * viscSpecies_[k]/m_spwork[k]; 
      }
    }
    
    return m_viscmix;
  }

  void LiquidTransport::getSpeciesViscosities(doublereal* visc) { 
    update_temp();
    if (!m_visc_temp_ok) {
      updateViscosity_temp();
    }
    copy(viscSpecies_.begin(), viscSpecies_.end(), visc); 
  }


  /******************* binary diffusion coefficients **************/


  void LiquidTransport::getBinaryDiffCoeffs(int ld, doublereal* d) {
    int i,j;

    update_temp();

    // if necessary, evaluate the binary diffusion coefficents
    // from the polynomial fits
    if (!m_diff_temp_ok) updateDiff_temp();
    doublereal pres = m_thermo->pressure();

    doublereal rp = 1.0/pres;
    for (i = 0; i < m_nsp; i++) 
      for (j = 0; j < m_nsp; j++) {
	d[ld*j + i] = rp * m_bdiff(i,j);
      }
  }


  void LiquidTransport::getMobilities(doublereal* const mobil) {
    // this needs to be checked out. 
    int k;
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (k = 0; k < m_nsp; k++) {
      mobil[k] = c1 * m_spwork[k] * m_thermo->charge(k);
    }
  } 
  

  
  void LiquidTransport::set_Grad_V(const doublereal* const grad_V) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_V[a] = grad_V[a];
    }
  }

  void LiquidTransport::set_Grad_T(const doublereal* const grad_T) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_T[a] = grad_T[a];
    }
  }

 void LiquidTransport::set_Grad_X(const doublereal* const grad_X) {
   int itop = m_nDim * m_nsp;
   for (int i = 0; i < itop; i++) {
     m_Grad_X[i] = grad_X[i];
   }
   update_Grad_lnAC();
 }


  /****************** thermal conductivity **********************/

  /*
   * The thermal conductivity is computed from the following mixture rule:
   * \[
   * \lambda = 0.5 \left( \sum_k X_k \lambda_k 
   * + \frac{1}{\sum_k X_k/\lambda_k}\right)
   * \]
   */
  doublereal LiquidTransport::thermalConductivity() {
   
    update_temp();
    update_conc();

    if (!m_cond_temp_ok) {
      updateCond_temp();
    } 
    if (!m_cond_mix_ok) {
      doublereal sum1 = 0.0, sum2 = 0.0;
      for (int k = 0; k < m_nsp; k++) {
	sum1 += m_molefracs[k] * m_cond[k];
	sum2 += m_molefracs[k] / m_cond[k];
      }
      m_lambda = 0.5*(sum1 + 1.0/sum2);
      m_cond_mix_ok = true;
    }

    return m_lambda;
  }


  /****************** thermal diffusion coefficients ************/

  /**
   * Thermal diffusion is not considered in this mixture-averaged
   * model. To include thermal diffusion, use transport manager
   * MultiTransport instead. This methods fills out array dt with
   * zeros.
   */
  void LiquidTransport::getThermalDiffCoeffs(doublereal* const dt) {
    for (int k = 0; k < m_nsp; k++) {
      dt[k] = 0.0;
    }
  }

  /**
   * @param ndim The number of spatial dimensions (1, 2, or 3).
   * @param grad_T The temperature gradient (ignored in this model).
   * @param ldx  Leading dimension of the grad_X array.
   * The diffusive mass flux of species \e k is computed from
   *
   * \f[
   *      \vec{j}_k = -n M_k D_k \nabla X_k.
   * \f]
   */
  void LiquidTransport::getSpeciesFluxes(int ndim, 
					 const doublereal* grad_T, 
					 int ldx, const doublereal* grad_X, 
					 int ldf, doublereal* fluxes) {
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesFluxesExt(ldf, fluxes);
  }

  /**
   * @param ndim The number of spatial dimensions (1, 2, or 3).
   * @param grad_T The temperature gradient (ignored in this model).
   * @param ldx  Leading dimension of the grad_X array.
   * The diffusive mass flux of species \e k is computed from
   *
   * \f[
   *      \vec{j}_k = -n M_k D_k \nabla X_k.
   * \f]
   */
  void LiquidTransport::getSpeciesFluxesExt(int ldf, doublereal* fluxes) {
    int n, k;

    update_temp();
    update_conc();


    getMixDiffCoeffs(DATA_PTR(m_spwork));


    const array_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();
    // Unroll wrt ndim
    vector_fp sum(m_nDim,0.0);
    for (n = 0; n < m_nDim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * m_Grad_X[n*m_nsp + k];
	sum[n] += fluxes[n*ldf + k];
      }
    }
    // add correction flux to enforce sum to zero
    for (n = 0; n < m_nDim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*ldf + k] -= y[k]*sum[n];
      }
    }
  }

  void LiquidTransport::getSpeciesDiffusiveMassFluxes(doublereal* const fluxes) {
    int n, k;

    update_temp();
    update_conc();


    getMixDiffCoeffs(DATA_PTR(m_spwork));


    const array_fp& mw = m_thermo->molecularWeights();
    const doublereal* const y  = m_thermo->massFractions();
    const doublereal rhon = m_thermo->molarDensity();
    // Unroll wrt ndim
    vector_fp sum(m_nDim,0.0);
    for (n = 0; n < m_nDim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*m_nsp + k] = -rhon * mw[k] * m_spwork[k] * m_Grad_X[n*m_nsp + k];
	sum[n] += fluxes[n*m_nsp + k];
      }
    }
    // add correction flux to enforce sum to zero
    for (n = 0; n < m_nDim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*m_nsp + k] -= y[k]*sum[n];
      }
    }
  }

  /**
   * Mixture-averaged diffusion coefficients [m^2/s]. 
   *
   * For the single species case or the pure fluid case
   * the routine returns the self-diffusion coefficient.
   * This is need to avoid a Nan result in the formula
   * below.
   */
  void LiquidTransport::getMixDiffCoeffs(doublereal* const d) {

    update_temp();
    update_conc();

    // update the binary diffusion coefficients if necessary
    if (!m_diff_temp_ok) {
      updateDiff_temp();
    }
 
    int k, j;
    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal sumxw_tran = 0.0;
    doublereal sum2;
 
    if (m_nsp == 1) {
      d[0] = m_bdiff(0,0);
    } else {
      for (k = 0; k < m_nsp; k++) {
	sumxw_tran += m_molefracs_tran[k] * m_mw[k];
      }
      for (k = 0; k < m_nsp; k++) {
	sum2 = 0.0;
	for (j = 0; j < m_nsp; j++) {
	  if (j != k) {
	    sum2 += m_molefracs_tran[j] / m_bdiff(j,k);
	  }
	}
	// Because we use m_molefracs_tran, sum2 must be positive definate
	// if (sum2 <= 0.0) {
	//  d[k] = m_bdiff(k,k);
	//  } else {
	d[k] = (sumxw_tran - m_molefracs_tran[k] * m_mw[k])/(mmw * sum2);
	//  }
      }
    }
  }

                 
  // Handles the effects of changes in the Temperature, internally
  // within the object.
  /*
   *  This is called whenever a transport property is
   *  requested.  
   *  The first task is to check whether the temperature has changed
   *  since the last call to update_temp().
   *  If it hasn't then an immediate return is carried out.
   *
   *     @internal
   */ 
  void LiquidTransport::update_temp()
  {
    // First make a decision about whether we need to recalculate
    doublereal t = m_thermo->temperature();
    if (t == m_temp) return;

    // Next do a reality check on temperature value
    if (t < 0.0) {
      throw CanteraError("LiquidTransport::update_temp()",
			 "negative temperature "+fp2str(t));
    }

    // Compute various direct functions of temperature
    m_temp = t;
    m_logt = log(m_temp);
    m_kbt = Boltzmann * m_temp;
    m_sqrt_t = sqrt(m_temp);
    m_t14 = sqrt(m_sqrt_t);
    m_t32 = m_temp * m_sqrt_t;
    m_sqrt_kbt = sqrt(Boltzmann*m_temp);

    // compute powers of log(T)
    // -> may move this
    m_polytempvec[0] = 1.0;
    m_polytempvec[1] = m_logt;
    m_polytempvec[2] = m_logt*m_logt;
    m_polytempvec[3] = m_logt*m_logt*m_logt;
    m_polytempvec[4] = m_logt*m_logt*m_logt*m_logt;

    // temperature has changed so temp flags are flipped
    m_visc_temp_ok  = false;
    m_diff_temp_ok  = false;

    // temperature has changed, so polynomial temperature 
    // interpolations will need to be reevaluated.
    // This means that many concentration 
    m_visc_conc_ok  = false;
    m_cond_temp_ok  = false;

    // Mixture stuff needs to be evaluated 
    m_visc_mix_ok = false;
    m_diff_mix_ok = false;
    //  m_cond_mix_ok = false; (don't need it because a lower lvl flag is set    

  }                 


  // Handles the effects of changes in the mixture concentration
  /*
   *   This is called for every interface call to check whether
   *   the concentrations have changed. Concentrations change
   *   whenever the pressure or the mole fraction has changed.
   *   If it has changed, the recalculations should be done.
   *
   *   Note this should be a lightweight function since it's
   *   part of all of the interfaces.
   *
   *   @internal
   */ 
  void LiquidTransport::update_conc() {
    // If the pressure has changed then the concentrations 
    // have changed.
    doublereal pres = m_thermo->pressure();
    bool qReturn = true;
    if (pres != m_press) {
      qReturn = false;
      m_press = pres;
    } 
    int iStateNew = m_thermo->stateMFNumber();
    if (iStateNew != m_iStateMF) {
      qReturn = false;
      m_thermo->getMoleFractions(DATA_PTR(m_molefracs));
      m_thermo->getConcentrations(DATA_PTR(m_concentrations));
      concTot_ = 0.0;
      concTot_tran_ = 0.0;
      for (int k = 0; k < m_nsp; k++) {
	m_molefracs[k] = fmaxx(0.0, m_molefracs[k]);
	m_molefracs_tran[k] = fmaxx(MIN_X, m_molefracs[k]);
	concTot_tran_ += m_molefracs_tran[k];
	concTot_ += m_concentrations[k];
      }
      dens_ = m_thermo->density();
      meanMolecularWeight_ =  m_thermo->meanMolecularWeight();
      concTot_tran_ *= concTot_;
    }
    if (qReturn) {
      return;
    }

    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.
    m_visc_conc_ok = false;
  
    // Mixture stuff needs to be evaluated
    m_visc_mix_ok = false;
    m_diff_mix_ok = false;
    m_cond_mix_ok = false;
  }


  // We formulate the directional derivative
  /*
   *     We only calculate the change in ac due to composition.
   *  The pressure and the temperature are taken care of in
   *  other parts of the expression.
   *
   */
  void LiquidTransport::update_Grad_lnAC() {
    int k;
    

    for (int a = 0; a < m_nDim; a++) {
      // We form the directional derivative
      double * ma_Grad_X = &m_Grad_X[a*m_nsp];
      double sum = 0.0;
      for (k = 0; k < m_nsp; k++) {
        sum += ma_Grad_X[k] * ma_Grad_X[k];
      }
      if (sum == 0.0) {
	for (k = 0; k < m_nsp; k++) {
	  m_Grad_lnAC[m_nsp * a + k] = 0.0;
	}
	continue;
      }
      double mag = 1.0E-7 / sum;
    
	for (k = 0; k < m_nsp; k++) {
	  Xdelta_[k] = m_molefracs[k] + mag * ma_Grad_X[k];
	  if (Xdelta_[k] > 1.0) {
	    Xdelta_[k] = 1.0;
	  }
	  if (Xdelta_[k] < 0.0) {
	    Xdelta_[k] = 0.0;
	  }
	}
      m_thermo->setMoleFractions(DATA_PTR(Xdelta_));
      m_thermo->getActivityCoefficients(DATA_PTR(lnActCoeffMolarDelta_));
      for (k = 0; k < m_nsp; k++) {
	lnActCoeffMolarDelta_[k] = log(lnActCoeffMolarDelta_[k]);
      }

      for (k = 0; k < m_nsp; k++) {
	m_Grad_lnAC[m_nsp * a + k] =
	  sum * (lnActCoeffMolarDelta_[k] - log(actCoeffMolar_[k])) / mag;
      }
    }
    m_thermo->setMoleFractions(DATA_PTR(m_molefracs));

  }

  /*************************************************************************
   *
   *    methods to update temperature-dependent properties
   *
   *************************************************************************/

  /**
   * Update the temperature-dependent parts of the mixture-averaged 
   * thermal conductivity. 
   */
  void LiquidTransport::updateCond_temp() {

    int k;
    if (m_mode == CK_Mode) {
      for (k = 0; k < m_nsp; k++) {
	m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
      }
    } else {
      for (k = 0; k < m_nsp; k++) {
	m_cond[k] = m_sqrt_t * dot5(m_polytempvec, m_condcoeffs[k]);
      }
    }
    m_cond_temp_ok = true;
    m_cond_mix_ok = false;
  }


  /**
   * Update the binary diffusion coefficients. These are evaluated
   * from the polynomial fits at unit pressure (1 Pa).
   */
  void LiquidTransport::updateDiff_temp() {

    // evaluate binary diffusion coefficients at unit pressure
    int i,j;
    int ic = 0;
    if (m_mode == CK_Mode) {
      for (i = 0; i < m_nsp; i++) {
	for (j = i; j < m_nsp; j++) {
	  m_bdiff(i,j) = exp(dot4(m_polytempvec, m_diffcoeffs[ic]));
	  m_bdiff(j,i) = m_bdiff(i,j);
	  ic++;
	}
      }
    }       
    else {
      for (i = 0; i < m_nsp; i++) {
	for (j = i; j < m_nsp; j++) {
	  m_bdiff(i,j) = m_temp * m_sqrt_t*dot5(m_polytempvec, 
						m_diffcoeffs[ic]);
	  m_bdiff(j,i) = m_bdiff(i,j);
	  ic++;
	}
      }
    }

    m_diff_temp_ok = true;
    m_diff_mix_ok = false;
  }


  /**
   * Update the pure-species viscosities.
   */
  void LiquidTransport::updateViscosities_conc() {
    m_visc_conc_ok = true;
  }


  /**
   * Update the temperature-dependent viscosity terms.
   * Updates the array of pure species viscosities, and the 
   * weighting functions in the viscosity mixture rule.
   * The flag m_visc_ok is set to true.
   */
  void LiquidTransport::updateViscosity_temp() {
    int k;
    doublereal vratiokj, wratiojk, factor1;

    if (m_mode == CK_Mode) {
      for (k = 0; k < m_nsp; k++) {
	viscSpecies_[k] = exp(dot4(m_polytempvec, viscCoeffsVector_[k]));
	m_sqvisc[k] = sqrt(viscSpecies_[k]);
      }
    }
    else {
      for (k = 0; k < m_nsp; k++) {
	// the polynomial fit is done for sqrt(visc/sqrt(T))
	m_sqvisc[k] = m_t14*dot5(m_polytempvec, viscCoeffsVector_[k]);
	viscSpecies_[k] = (m_sqvisc[k]*m_sqvisc[k]);
      }
    }

    // see Eq. (9-5.15) of Reid, Prausnitz, and Poling
    int j;
    for (j = 0; j < m_nsp; j++) {
      for (k = j; k < m_nsp; k++) {
	vratiokj = viscSpecies_[k]/viscSpecies_[j];
	wratiojk = m_mw[j]/m_mw[k];

	// Note that m_wratjk(k,j) holds the square root of
	// m_wratjk(j,k)!
	factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
	m_phi(k,j) = factor1*factor1 /
	  (SqrtEight * m_wratkj1(j,k)); 
	m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
      }
    }

    m_visc_temp_ok = true;
    m_visc_mix_ok = false;
  }


  /*
   *
   *    Solve for the diffusional velocities in the Stefan-Maxwell equations
   *
   */
  void LiquidTransport::stefan_maxwell_solve() {
    int i, j, a;
    doublereal tmp;
    int VIM = m_nDim;
    m_B.resize(m_nsp, VIM);
    //! grab a local copy of the molecular weights
    const vector_fp& M =  m_thermo->molecularWeights();
    
 
    /*
     * Update the concentrations in the mixture.
     */
    update_conc();

    double T = m_thermo->temperature();

 
    m_thermo->getStandardVolumes(DATA_PTR(volume_specPM_));
    m_thermo->getActivityCoefficients(DATA_PTR(actCoeffMolar_));

    /* 
     *  Calculate the electrochemical potential gradient. This is the
     *  driving force for relative diffusional transport.
     *
     *  Here we calculate
     *
     *          c_i * (grad (mu_i) + S_i grad T - M_i / dens * grad P
     *
     *   This is  Eqn. 13-1 p. 318 Newman. The original equation is from
     *   Hershfeld, Curtis, and Bird.
     *
     *   S_i is the partial molar entropy of species i. This term will cancel
     *   out a lot of the grad T terms in grad (mu_i), therefore simplifying
     *   the expression.
     *
     *  Ok I think there may be many ways to do this. One way is to do it via basis
     *  functions, at the nodes, as a function of the variables in the problem.
     *
     *  For calculation of molality based thermo systems, we current get
     *  the molar based values. This may change.
     *
     *  Note, we have broken the symmetry of the matrix here, due to 
     *  consideratins involving species concentrations going to zero.
     *
     */
    for (i = 0; i < m_nsp; i++) {
      double xi_denom = m_molefracs_tran[i];
      for (a = 0; a < VIM; a++) {
	m_ck_Grad_mu[a*m_nsp + i] =
	  m_chargeSpecies[i] * concTot_ * Faraday * m_Grad_V[a]
	  + concTot_ * (volume_specPM_[i] - M[i]/dens_) * m_Grad_P[a]
	  + concTot_ * GasConstant * T * m_Grad_lnAC[a*m_nsp+i] / actCoeffMolar_[i]
	  + concTot_ * GasConstant * T * m_Grad_X[a*m_nsp+i] / xi_denom;
      }
    }

    if (m_thermo->activityConvention() == cAC_CONVENTION_MOLALITY) {
      int iSolvent = 0;
      double mwSolvent = m_thermo->molecularWeight(iSolvent);
      double mnaught = mwSolvent/ 1000.;
      double lnmnaught = log(mnaught);
      for (i = 1; i < m_nsp; i++) {
	for (a = 0; a < VIM; a++) {
	  m_ck_Grad_mu[a*m_nsp + i] -=
	    m_concentrations[i] * GasConstant * m_Grad_T[a] * lnmnaught;
	}
      }
    }

    /*
     * Just for Note, m_A(i,j) refers to the ith row and jth column.
     * They are still fortran ordered, so that i varies fastest.
     */
    switch (VIM) {
    case 1:  /* 1-D approximation */
      m_B(0,0) = 0.0;
      for (j = 0; j < m_nsp; j++) {
	m_A(0,j) = M[j] * m_concentrations[j];
      }
      for (i = 1; i < m_nsp; i++){
	m_B(i,0) = m_ck_Grad_mu[i] / (GasConstant * T);
	m_A(i,i) = 0.0;
	for (j = 0; j < m_nsp; j++){
	  if (j != i) {
	    tmp = m_concentrations[j] / m_DiffCoeff_StefMax(i,j);
	    m_A(i,i) +=   tmp;
	    m_A(i,j)  = - tmp;
	  }
	}
      }

      //! invert and solve the system  Ax = b. Answer is in m_B
      solve(m_A, m_B);
  	
      break;
    case 2:  /* 2-D approximation */
      m_B(0,0) = 0.0;
      m_B(0,1) = 0.0;
      for (j = 0; j < m_nsp; j++) {
	m_A(0,j) = M[j] * m_concentrations[j];
      }
      for (i = 1; i < m_nsp; i++){
	m_B(i,0) =  m_ck_Grad_mu[i]         / (GasConstant * T);
	m_B(i,1) =  m_ck_Grad_mu[m_nsp + i] / (GasConstant * T);
	m_A(i,i) = 0.0;
	for (j = 0; j < m_nsp; j++) {
	  if (j != i) {
	    tmp =  m_concentrations[j] / m_DiffCoeff_StefMax(i,j);
	    m_A(i,i) +=   tmp;
	    m_A(i,j)  = - tmp;
	  }
	}
      }

      //! invert and solve the system  Ax = b. Answer is in m_B
      solve(m_A, m_B);
	 
 	
      break;

    case 3:  /* 3-D approximation */
      m_B(0,0) = 0.0;
      m_B(0,1) = 0.0;
      m_B(0,2) = 0.0;
      for (j = 0; j < m_nsp; j++) {
	m_A(0,j) = M[j] * m_concentrations[j];
      }
      for (i = 1; i < m_nsp; i++){
	m_B(i,0) = m_ck_Grad_mu[i]           / (GasConstant * T);
	m_B(i,1) = m_ck_Grad_mu[m_nsp + i]   / (GasConstant * T);
	m_B(i,2) = m_ck_Grad_mu[2*m_nsp + i] / (GasConstant * T);
	m_A(i,i) = 0.0;
	for (j = 0; j < m_nsp; j++) {
	  if (j != i) {
	    tmp =  m_concentrations[j] / m_DiffCoeff_StefMax(i,j);
	    m_A(i,i) +=   tmp;
	    m_A(i,j)  = - tmp;
	  }
	}
      }

      //! invert and solve the system  Ax = b. Answer is in m_B
      solve(m_A, m_B);

      break;
    default:
      printf("uninmplemetnd\n");
      throw CanteraError("routine", "not done");
      break;
    }

    for (a = 0; a < VIM; a++) {
      for (j = 0; j < m_nsp; j++) {
	m_flux(j,a) =  M[j] * m_concentrations[j] * m_B(j,a);
      }
    }
  }
}
