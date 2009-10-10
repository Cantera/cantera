/**
 *  @file SimpleTransport.cpp
 *  Simple mostly constant transport properties
 */
/* 
 * $Revision: 1.10 $
 * $Date: 2009/03/24 20:44:30 $
 */

#include "ThermoPhase.h"
#include "SimpleTransport.h"

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
  //================================================================================================
  SimpleTransport::SimpleTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim),
    m_nsp(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false)
  {
  }
  //================================================================================================
  SimpleTransport::SimpleTransport(const SimpleTransport &right) :
    Transport(),
    m_nsp(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = right;
  }
  //================================================================================================
  SimpleTransport& SimpleTransport::operator=(const SimpleTransport& right) {
    if (&right != this) {
      return *this; 
    }
    Transport::operator=(right);
    m_nsp                                 = right.m_nsp;
    m_tmin                                = right.m_tmin;
    m_tmax                                = right.m_tmax;
    m_mw                                  = right.m_mw;

    m_Grad_X                              = right.m_Grad_X;
    m_Grad_T                              = right.m_Grad_T;
    m_Grad_V                              = right.m_Grad_V;

    m_viscSpecies                         = right.m_viscSpecies;
    m_condSpecies                         = right.m_condSpecies;
    m_iStateMF = -1;
    m_molefracs                           = right.m_molefracs;
    m_concentrations                      = right.m_concentrations;
    m_chargeSpecies                       = right.m_chargeSpecies;
 

    m_temp                                = right.m_temp;
    m_press                               = right.m_press;
    m_lambda                              = right.m_lambda;
    m_viscmix                             = right.m_viscmix;
    m_spwork                              = right.m_spwork;
    m_visc_mix_ok    = false;
    m_visc_temp_ok   = false;
    m_diff_mix_ok    = false;
    m_diff_temp_ok   = false;
    m_cond_temp_ok   = false;
    m_cond_mix_ok    = false;
    m_nDim                                = right.m_nDim;

    return *this; 
  }

  //================================================================================================
  Transport *SimpleTransport::duplMyselfAsTransport() const {
    SimpleTransport* tr = new SimpleTransport(*this);
    return (dynamic_cast<Transport *>(tr));
  }
  //================================================================================================
  // Initialize the object
  /*
   *  This is where we dimension everything.
   */
  bool SimpleTransport::initLiquid(LiquidTransportParams& tr) {

    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(), 
	 m_thermo->molecularWeights().end(), m_mw.begin());

    //save logarithm of pre-exponential for easier computation

    //m_diffcoeffs = tr.diffcoeffs;


    m_viscSpecies.resize(m_nsp);
    m_condSpecies.resize(m_nsp);


    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);



    // set all flags to false
    m_visc_mix_ok   = false;
    m_visc_temp_ok  = false;

    m_cond_temp_ok = false;
    m_cond_mix_ok  = false;

    m_diff_temp_ok   = false;
    m_diff_mix_ok  = false;

    return true;
  }

  //================================================================================================
  // Returns the mixture viscosity of the solution
  /*
   * The viscosity is computed using the general mixture rules
   * specified in the variable compositionDepType_.
   * 
   * Solvent-only:
   *    \f[
   *         \mu = \mu_0
   *    \f]
   * Mixture-average:
   *    \f[
   *         \mu = \sum_k {\mu_k X_k}
   *    \f]
   *  
   * Here \f$ \mu_k \f$ is the viscosity of pure species \e k.
   *
   * @see updateViscosity_T();
   */ 
  doublereal SimpleTransport::viscosity() {
        
    update_T();
    update_C();

    if (m_visc_mix_ok) return m_viscmix;
  
    // update m_viscSpecies[] if necessary
    if (!m_visc_temp_ok) {
      updateViscosity_T();
    }

    if (compositionDepType_ == 0) {
      m_viscmix = m_viscSpecies[0];
    } else if (compositionDepType_ == 1) {
      m_viscmix = 0.0;
      for (int k = 0; k < m_nsp; k++) {
	m_viscmix += m_viscSpecies[k] * m_molefracs[k];
      }
    }
    m_visc_mix_ok = true;
    return m_viscmix;
  }
  //================================================================================================
  void SimpleTransport::getSpeciesViscosities(doublereal* visc) { 
    update_T();
    if (!m_visc_temp_ok) {
      updateViscosity_T();
    }
    copy(m_viscSpecies.begin(), m_viscSpecies.end(), visc); 
  }
  //================================================================================================
  void SimpleTransport::getBinaryDiffCoeffs(int ld, doublereal* d) {
    int i, j;
    double bdiff;
    update_T();

    // if necessary, evaluate the species diffusion coefficents
    // from the polynomial fits
    if (!m_diff_temp_ok) updateDiff_T();
 
    for (i = 0; i < m_nsp; i++) {
      for (j = 0; j < m_nsp; j++) {
        bdiff = 0.5 * (m_diffSpecies[i] + m_diffSpecies[j]);
	d[i*m_nsp+j] = bdiff;
      }
    }
  }
  //================================================================================================
  void SimpleTransport::getMobilities(doublereal* const mobil) {
    // this needs to be checked out. 
    int k;
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (k = 0; k < m_nsp; k++) {
      mobil[k] = c1 * m_spwork[k] * m_thermo->charge(k);
    }
  } 
  //================================================================================================
  void SimpleTransport::set_Grad_V(const doublereal* const grad_V) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_V[a] = grad_V[a];
    }
  }
  //================================================================================================
  void SimpleTransport::set_Grad_T(const doublereal* const grad_T) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_T[a] = grad_T[a];
    }
  }
  //================================================================================================
  void SimpleTransport::set_Grad_X(const doublereal* const grad_X) {
    int itop = m_nDim * m_nsp;
    for (int i = 0; i < itop; i++) {
      m_Grad_X[i] = grad_X[i];
    }
  }

  //================================================================================================
  // Returns the mixture thermal conductivity of the solution
  /*
   * The thermal is computed using the general mixture rules
   * specified in the variable compositionDepType_.
   * 
   * Solvent-only:
   *    \f[
   *         \lambda = \lambda_0
   *    \f]
   * Mixture-average:
   *    \f[
   *         \lambda = \sum_k {\lambda_k X_k}
   *    \f]
   *  
   * Here \f$ \lambda_k \f$ is the thermal conductivity of pure species \e k.
   *
   * @see updateCond_T();
   */ 
  doublereal SimpleTransport::thermalConductivity() {
    update_T();
    update_C();
    if (!m_cond_temp_ok) {
      updateCond_T();
    } 
    if (!m_cond_mix_ok) {
      if (compositionDepType_ == 0) {
	m_lambda = m_condSpecies[0];
      } else if (compositionDepType_ == 1) {
	m_lambda = 0.0;
	for (int k = 0; k < m_nsp; k++) {
	  m_lambda += m_condSpecies[k] * m_molefracs[k];
	}
      }
      m_cond_mix_ok = true;
    }
    return m_lambda;
  }
 //================================================================================================

  /*
   * Thermal diffusion is not considered in this mixture-averaged
   * model. To include thermal diffusion, use transport manager
   * MultiTransport instead. This methods fills out array dt with
   * zeros.
   */
  void SimpleTransport::getThermalDiffCoeffs(doublereal* const dt) {
    for (int k = 0; k < m_nsp; k++) {
      dt[k] = 0.0;
    }
  }
//================================================================================================
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
  void SimpleTransport::getSpeciesFluxes(int ndim, 
					 const doublereal* grad_T, 
					 int ldx, const doublereal* grad_X, 
					 int ldf, doublereal* fluxes) {
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesFluxesExt(ldf, fluxes);
  }
  //================================================================================================
  //  Return the species diffusive mass fluxes wrt to
  //  the mass averaged velocity,
  /*
   *
   *  units = kg/m2/s
   *
   * Internally, gradients in the in mole fraction, temperature
   * and electrostatic potential contribute to the diffusive flux
   *  
   *
   * The diffusive mass flux of species \e k is computed from the following 
   * formula
   *
   *    \f[
   *         j_k = - \rho M_k D_k \nabla X_k - Y_k V_c
   *    \f]
   *
   *    where V_c is the correction velocity
   *
   *    \f[
   *         V_c =  - \sum_j {\rho M_j D_j \nabla X_j}
   *    \f]
   *
   *  @param ldf     stride of the fluxes array. Must be equal to
   *                 or greater than the number of species.
   *  @param fluxes  Vector of calculated fluxes
   */
  void SimpleTransport::getSpeciesFluxesExt(int ldf, doublereal* fluxes) {
    int n, k;
    AssertThrow(ldf >= m_nsp ,"SimpleTransport::getSpeciesFluxesExt: Stride must be greater than m_nsp");
    update_T();
    update_C();

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
  //================================================================================================
  // Mixture-averaged diffusion coefficients [m^2/s]. 
  /*
   *  Returns the simple diffusion coefficients input into the model. Nothing fancy here.
   */
  void SimpleTransport::getMixDiffCoeffs(doublereal* const d) {
    update_T();
    update_C();
    // update the binary diffusion coefficients if necessary
    if (!m_diff_temp_ok) {
      updateDiff_T();
    }  
    for (int k = 0; k < m_nsp; k++) {
      d[k] = m_diffSpecies[k];
    }
  }
  //================================================================================================
           
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
  bool SimpleTransport::update_C() {
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
      for (int k = 0; k < m_nsp; k++) {
	m_molefracs[k] = fmaxx(0.0, m_molefracs[k]);
	concTot_ += m_concentrations[k];
      }
      dens_ = m_thermo->density();
      meanMolecularWeight_ =  m_thermo->meanMolecularWeight();
    }
    if (qReturn) {
      return false;
    }

  
    // Mixture stuff needs to be evaluated
    m_visc_mix_ok = false;
    m_diff_mix_ok = false;
    m_cond_mix_ok = false;

    return true;
  }

  //================================================================================================
  /**
   * Update the temperature-dependent parts of the mixture-averaged 
   * thermal conductivity. 
   */
  void SimpleTransport::updateCond_T() {
    int k;
    if (tempDepType_ == 0) {
      for (k = 0; k < m_nsp; k++) {
	Coeff_T_ &coeff = m_coeffLambda_Ns[k];
	m_condSpecies[k] = coeff[0];
      }
    } else if (tempDepType_ == 1) {
      for (k = 0; k < m_nsp; k++) {
	Coeff_T_ &coeff = m_coeffLambda_Ns[k];
	m_condSpecies[k] = coeff[0] * pow(m_temp,coeff[1]) * exp(-coeff[2]/m_temp);
      }
    }
    m_cond_temp_ok = true;
    m_cond_mix_ok = false;
  }
  //================================================================================================
  /**
   * Update the species diffusion coefficients.
   */
  void SimpleTransport::updateDiff_T() {
    int k;
    if (tempDepType_ == 0) {
      for (k = 0; k < m_nsp; k++) {
	Coeff_T_ &coeff = m_coeffDiff_Ns[k];
	m_diffSpecies[k] = coeff[0];
      }
    } else if (tempDepType_ == 1) {
      for (k = 0; k < m_nsp; k++) {
	Coeff_T_ &coeff = m_coeffDiff_Ns[k];
	m_viscSpecies[k] = coeff[0] * pow(m_temp,coeff[1]) * exp(-coeff[2]/m_temp);
      }
    }
    m_diff_temp_ok = true;
    m_diff_mix_ok = false;
  }
  //================================================================================================

  /**
   * Update the pure-species viscosities.
   */
  void SimpleTransport::updateViscosities_C() {

  }
  //================================================================================================
  /**
   * Update the temperature-dependent viscosity terms.
   * Updates the array of pure species viscosities, and the 
   * weighting functions in the viscosity mixture rule.
   * The flag m_visc_ok is set to true.
   */
  void SimpleTransport::updateViscosity_T() {
    int k;
    if (tempDepType_ == 0) {
      for (k = 0; k < m_nsp; k++) {
	Coeff_T_ &coeff = m_coeffVisc_Ns[k];
	m_viscSpecies[k] = coeff[0];
      }
    } else if (tempDepType_ == 1) {
      for (k = 0; k < m_nsp; k++) {
	Coeff_T_ &coeff = m_coeffVisc_Ns[k];
	m_viscSpecies[k] = coeff[0] * pow(m_temp,coeff[1]) * exp(-coeff[2]/m_temp);
      }
    }
    m_visc_temp_ok = true;
    m_visc_mix_ok = false;
  }
  //=================================================================================================
  bool SimpleTransport::update_T()
  {
    doublereal t = m_thermo->temperature();
    if (t == m_temp) return false;
    if (t < 0.0) {
      throw CanteraError("SimpleTransport::update_T",
                         "negative temperature "+fp2str(t));
    }

    // Compute various functions of temperature
    m_temp = t;
  
    // temperature has changed, so polynomial temperature
    // interpolations will need to be reevaluated.
    // Set all of these flags to false
    m_visc_mix_ok = false;
    m_visc_temp_ok  = false;

    m_cond_temp_ok = true;
    m_cond_mix_ok = false;

    m_diff_mix_ok = false;
    m_diff_temp_ok = false;

    return true;
  }

  /**
   * Throw an exception if this method is invoked. 
   * This probably indicates something is not yet implemented.
   */
  doublereal SimpleTransport::err(std::string msg) const {
    throw CanteraError("SimpleTransport Class",
		       "\n\n\n**** Method "+ msg +" not implemented in model "
		       + int2str(model()) + " ****\n"
		       "(Did you forget to specify a transport model?)\n\n\n");
      
    return 0.0;
  }
  //================================================================================================

}
//================================================================================================
