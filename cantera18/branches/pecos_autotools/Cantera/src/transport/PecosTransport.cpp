/**
 *  @file PecosTransport.cpp
 *  Mixture-averaged transport properties.
 */

/* $Author$
 * $Revision$
 * $Date$
 */

#include "ThermoPhase.h"
#include "PecosTransport.h"

#include "utilities.h"
#include "TransportParams.h"
#include "TransportFactory.h"

#include <iostream>
using namespace std;

/** 
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#define MIN_X 1.e-20

namespace Cantera {

  //////////////////// class PecosTransport methods //////////////


  PecosTransport::PecosTransport() :
    m_nsp(0),
    m_tmin(-1.0),
    m_tmax(100000.),
    m_temp(-1.0),
    m_logt(0.0)
  {


  }

  bool PecosTransport::initGas( GasTransportParams& tr ) {

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
    m_visccoeffs = tr.visccoeffs;
    m_condcoeffs = tr.condcoeffs;
    m_diffcoeffs = tr.diffcoeffs;

    m_zrot       = tr.zrot;
    m_crot       = tr.crot;
    m_epsilon    = tr.epsilon;
    m_mode       = tr.mode_;
    m_diam       = tr.diam;
    m_eps        = tr.eps;
    m_alpha      = tr.alpha;
    m_dipoleDiag.resize(m_nsp);
    for (int i = 0; i < m_nsp; i++) {
      m_dipoleDiag[i] = tr.dipole(i,i);
    }

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
    m_visc.resize(m_nsp);
    m_sqvisc.resize(m_nsp);
    m_cond.resize(m_nsp);
    m_bdiff.resize(m_nsp, m_nsp);

    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);

    // set flags all false
    m_viscmix_ok = false;
    m_viscwt_ok = false;
    m_spvisc_ok = false;
    m_spcond_ok = false;
    m_condmix_ok = false;
    m_spcond_ok = false;
    m_diffmix_ok = false;
    m_abc_ok = false;

    return true;
  }


  /*********************************************************
   *
   *                Public methods
   *
   *********************************************************/


  /******************  viscosity ******************************/

  /**
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
  doublereal PecosTransport::viscosity() {
        
    update_T();
    update_C();

    if (m_viscmix_ok) return m_viscmix;

    doublereal vismix = 0.0;
    int k;
    // update m_visc and m_phi if necessary
    if (!m_viscwt_ok) updateViscosity_T();

    multiply(m_phi, DATA_PTR(m_molefracs), DATA_PTR(m_spwork));

    for (k = 0; k < m_nsp; k++) {
      vismix += m_molefracs[k] * m_visc[k]/m_spwork[k]; //denom;
    }
    m_viscmix = vismix;
    return vismix;
  }


  /******************* binary diffusion coefficients **************/


  void PecosTransport::getBinaryDiffCoeffs(const int ld, doublereal* const d) {
    int i,j;

    update_T();

    // if necessary, evaluate the binary diffusion coefficents
    // from the polynomial fits
    if (!m_bindiff_ok) updateDiff_T();

    doublereal rp = 1.0/pressure_ig();
    for (i = 0; i < m_nsp; i++) 
      for (j = 0; j < m_nsp; j++) {
	d[ld*j + i] = rp * m_bdiff(i,j);
      }
  }


  void PecosTransport::getMobilities(doublereal* const mobil) {
    int k;
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (k = 0; k < m_nsp; k++) {
      mobil[k] = c1 * m_spwork[k] * m_thermo->charge(k);
    }
  } 
        

  /****************** thermal conductivity **********************/

  /**
   * The thermal conductivity is computed from the following mixture rule:
   * \[
   * \lambda = 0.5 \left( \sum_k X_k \lambda_k 
   * + \frac{1}{\sum_k X_k/\lambda_k}\right)
   * \]
   */
  doublereal PecosTransport::thermalConductivity() {
    int k;

    update_T();
    update_C();

    if (!m_spcond_ok)  updateCond_T(); 
    if (!m_condmix_ok) {
      doublereal sum1 = 0.0, sum2 = 0.0;
      for (k = 0; k < m_nsp; k++) {
	sum1 += m_molefracs[k] * m_cond[k];
	sum2 += m_molefracs[k] / m_cond[k];
      }
      m_lambda = 0.5*(sum1 + 1.0/sum2);
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
  void PecosTransport::getThermalDiffCoeffs(doublereal* const dt) {
    int k;
    for (k = 0; k < m_nsp; k++) {
      dt[k] = 0.0;
    }
  }

  /**
   * @param ndim The number of spatial dimensions (1, 2, or 3).
   * @param grad_T The temperature gradient (ignored in this model).
   * @param ldx  Leading dimension of the grad_X array.
   * The diffusive mass flux of species \e k is computed from
   * \f[
   * \vec{j}_k = -n M_k D_k \nabla X_k.
   * \f]
   */
  void PecosTransport::getSpeciesFluxes(int ndim, 
				      const doublereal* grad_T, int ldx, const doublereal* grad_X, 
				      int ldf, doublereal* fluxes) {
    int n, k;

    update_T();
    update_C();

    getMixDiffCoeffs(DATA_PTR(m_spwork));

    const array_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();

    vector_fp sum(ndim,0.0);
    for (n = 0; n < ndim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
	sum[n] += fluxes[n*ldf + k];
      }
    }
    // add correction flux to enforce sum to zero
    for (n = 0; n < ndim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*ldf + k] -= y[k]*sum[n];
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
  void PecosTransport::getMixDiffCoeffs(doublereal* const d) {

    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) updateDiff_T();

    int k, j;
    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal sumxw = 0.0, sum2;
    doublereal p = pressure_ig();
    if (m_nsp == 1) {
      d[0] = m_bdiff(0,0) / p;
    } else {
      for (k = 0; k < m_nsp; k++) sumxw += m_molefracs[k] * m_mw[k];
      for (k = 0; k < m_nsp; k++) {
	sum2 = 0.0;
	for (j = 0; j < m_nsp; j++) {
	  if (j != k) {
	    sum2 += m_molefracs[j] / m_bdiff(j,k);
	  }
	}
	if (sum2 <= 0.0) {
	  d[k] = m_bdiff(k,k) / p;
	} else {
	  d[k] = (sumxw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
	}
      }
    }
  }

                 
  /**
   *  @internal This is called whenever a transport property is
   *  requested from ThermoSubstance if the temperature has changed
   *  since the last call to update_T.
   */ 
  void PecosTransport::update_T() 
  {
    doublereal t = m_thermo->temperature();
    if (t == m_temp) return;
    if (t < 0.0) {
      throw CanteraError("PecosTransport::update_T",
			 "negative temperature "+fp2str(t));
    }
    m_temp = t;
    m_logt = log(m_temp);
    m_kbt = Boltzmann * m_temp;
    m_sqrt_t = sqrt(m_temp);
    m_t14 = sqrt(m_sqrt_t);
    m_t32 = m_temp * m_sqrt_t;
    m_sqrt_kbt = sqrt(Boltzmann*m_temp);

    // compute powers of log(T)
    m_polytempvec[0] = 1.0;
    m_polytempvec[1] = m_logt;
    m_polytempvec[2] = m_logt*m_logt;
    m_polytempvec[3] = m_logt*m_logt*m_logt;
    m_polytempvec[4] = m_logt*m_logt*m_logt*m_logt;

    // temperature has changed, so polynomial fits will need to be
    // redone.
    m_viscmix_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_spcond_ok = false;
    m_diffmix_ok = false;
    m_bindiff_ok = false;
    m_abc_ok  = false;
    m_condmix_ok = false;                 
  }                 

  /**
   *  @internal This is called the first time any transport property
   *  is requested from Mixture after the concentrations
   *  have changed.
   */ 
  void PecosTransport::update_C()  
  {
    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_viscmix_ok = false;
    m_diffmix_ok = false;
    m_condmix_ok = false;

    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // add an offset to avoid a pure species condition
    int k;
    for (k = 0; k < m_nsp; k++) {
      m_molefracs[k] = fmaxx(MIN_X, m_molefracs[k]);
    }
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
  void PecosTransport::updateCond_T() {

    int k;
    if (m_mode == CK_Mode) {
      for (k = 0; k < m_nsp; k++) {
	m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
      }
    }
    else {
      for (k = 0; k < m_nsp; k++) {
	m_cond[k] = m_sqrt_t*dot5(m_polytempvec, m_condcoeffs[k]);
      }
    }
    m_spcond_ok = true;
    m_condmix_ok = false;
  }


  /**
   * Update the binary diffusion coefficients. These are evaluated
   * from the polynomial fits at unit pressure (1 Pa).
   */
  void PecosTransport::updateDiff_T() {

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

    m_bindiff_ok = true;
    m_diffmix_ok = false;
  }


  /**
   * Update the pure-species viscosities.
   */
  void PecosTransport::updateSpeciesViscosities() {

    int k;
    if (m_mode == CK_Mode) {
      for (k = 0; k < m_nsp; k++) {
	m_visc[k] = exp(dot4(m_polytempvec, m_visccoeffs[k]));
	m_sqvisc[k] = sqrt(m_visc[k]);
      }
    }
    else {
      for (k = 0; k < m_nsp; k++) {
	// the polynomial fit is done for sqrt(visc/sqrt(T))
	m_sqvisc[k] = m_t14*dot5(m_polytempvec, m_visccoeffs[k]);
	m_visc[k] = (m_sqvisc[k]*m_sqvisc[k]);
      }
    }
    m_spvisc_ok = true;
  }


  /**
   * Update the temperature-dependent viscosity terms.
   * Updates the array of pure species viscosities, and the 
   * weighting functions in the viscosity mixture rule.
   * The flag m_visc_ok is set to true.
   */
  void PecosTransport::updateViscosity_T() {
    doublereal vratiokj, wratiojk, factor1;

    if (!m_spvisc_ok) updateSpeciesViscosities();

    // see Eq. (9-5.15) of Reid, Prausnitz, and Poling
    int j, k;
    for (j = 0; j < m_nsp; j++) {
      for (k = j; k < m_nsp; k++) {
	vratiokj = m_visc[k]/m_visc[j];
	wratiojk = m_mw[j]/m_mw[k];

	// Note that m_wratjk(k,j) holds the square root of
	// m_wratjk(j,k)!
	factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
	m_phi(k,j) = factor1*factor1 /
	  (SqrtEight * m_wratkj1(j,k)); 
	m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
      }
    }
    m_viscwt_ok = true;
  }

  /**
   * This function returns a Transport data object for a given species.
   *
   */
  struct GasTransportData PecosTransport::
    getGasTransportData(int kSpecies) 
  {
    struct GasTransportData td;
    td.speciesName = m_thermo->speciesName(kSpecies);

    td.geometry = 2;
    if (m_crot[kSpecies] == 0.0) {
      td.geometry = 0;
    } else if (m_crot[kSpecies] == 1.0) {
      td.geometry = 1;
    }
    td.wellDepth = m_eps[kSpecies] / Boltzmann;
    td.dipoleMoment = m_dipoleDiag[kSpecies] * 1.0E25 / SqrtTen;
    td.diameter = m_diam(kSpecies, kSpecies) * 1.0E10;
    td.polarizability = m_alpha[kSpecies] * 1.0E30;
    td.rotRelaxNumber = m_zrot[kSpecies];

    return td;
  }

}

