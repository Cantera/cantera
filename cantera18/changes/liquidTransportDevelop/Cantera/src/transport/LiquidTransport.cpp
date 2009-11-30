/**
 *  @file LiquidTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
/* 
 * $Revision$
 * $Date$
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
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_visc_conc_ok(false),
    m_radi_mix_ok(false),
    m_radi_temp_ok(false),
    m_radi_conc_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false),
    m_mode(-1000),
    m_debug(false),
    m_nDim(1)
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
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_visc_mix_ok(false),
    m_visc_temp_ok(false),
    m_visc_conc_ok(false),
    m_radi_mix_ok(false),
    m_radi_temp_ok(false),
    m_radi_conc_ok(false),
    m_diff_mix_ok(false),
    m_diff_temp_ok(false),
    m_cond_temp_ok(false),
    m_cond_mix_ok(false),
    m_mode(-1000),
    m_debug(false),
    m_nDim(1)
  {
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy construtor.
     */
    *this = right;
  }

  LiquidTransport& LiquidTransport::operator=(const LiquidTransport& right) {
    if (&right == this) {
      return *this; 
    }
    Transport::operator=(right);
    m_nsp                                 = right.m_nsp;
    m_tmin                                = right.m_tmin;
    m_tmax                                = right.m_tmax;
    m_mw                                  = right.m_mw;
    m_viscTempDep_Ns                      = right.m_viscTempDep_Ns;
    m_lambdaTempDep_Ns                    = right.m_lambdaTempDep_Ns;
    m_diffTempDep_Ns                      = right.m_diffTempDep_Ns;
    m_radiusTempDep_Ns                    = right.m_radiusTempDep_Ns;
    m_hydrodynamic_radius                 = right.m_hydrodynamic_radius;
    m_Grad_X                              = right.m_Grad_X;
    m_Grad_T                              = right.m_Grad_T;
    m_Grad_V                              = right.m_Grad_V;
    m_Grad_mu                             = right.m_Grad_mu;
    m_bdiff                               = right.m_bdiff;
    m_viscSpecies                         = right.m_viscSpecies;
    m_hydrodynamic_radius                 = right.m_hydrodynamic_radius;
    m_lambdaSpecies                       = right.m_lambdaSpecies;
    m_viscMixModel                        = right.m_viscMixModel;
    m_lambdaMixModel                      = right.m_lambdaMixModel;
    m_diffMixModel                        = right.m_diffMixModel;
    m_iStateMF = -1;
    m_molefracs                           = right.m_molefracs;
    m_molefracs_tran                      = right.m_molefracs_tran;
    m_concentrations                      = right.m_concentrations;
    m_actCoeff                            = right.m_actCoeff;
    m_Grad_lnAC                           = right.m_Grad_lnAC;
    m_chargeSpecies                       = right.m_chargeSpecies;
    m_volume_spec                         = right.m_volume_spec;
    m_B                                   = right.m_B;
    m_A                                   = right.m_A;
    m_temp                                = right.m_temp;
    m_logt                                = right.m_logt;
    m_press                               = right.m_press;
    m_flux                                = right.m_flux;
    m_Vdiff                               = right.m_Vdiff;
    m_lambda                              = right.m_lambda;
    m_viscmix                             = right.m_viscmix;
    m_spwork                              = right.m_spwork;
    m_visc_mix_ok    = false;
    m_visc_temp_ok   = false;
    m_visc_conc_ok   = false;
    m_radi_mix_ok    = false;
    m_radi_temp_ok   = false;
    m_radi_conc_ok   = false;
    m_diff_mix_ok    = false;
    m_diff_temp_ok   = false;
    m_cond_temp_ok   = false;
    m_cond_mix_ok    = false;
    m_mode                                = right.m_mode;
    m_debug                               = right.m_debug;
    m_nDim                                = right.m_nDim;

    return *this; 
  }


  Transport *LiquidTransport::duplMyselfAsTransport() const {
    LiquidTransport* tr = new LiquidTransport(*this);
    return (dynamic_cast<Transport *>(tr));
  }

   LiquidTransport::~LiquidTransport() {

     //These are constructed in TransportFactory::newLTP
     for ( int k = 0; k < m_nsp; k++) {
       if ( m_viscTempDep_Ns[k]   ) delete m_viscTempDep_Ns[k];
       if ( m_lambdaTempDep_Ns[k] ) delete m_lambdaTempDep_Ns[k];
       if ( m_radiusTempDep_Ns[k] ) delete m_radiusTempDep_Ns[k];
       if ( m_diffTempDep_Ns[k]   ) delete m_diffTempDep_Ns[k];
     }
     //These are constructed in TransportFactory::newLTI
     if ( m_viscMixModel   ) delete m_viscMixModel;
     if ( m_lambdaMixModel ) delete m_lambdaMixModel;
     if ( m_diffMixModel   ) delete m_diffMixModel;
     if ( m_radiusMixModel ) delete m_radiusMixModel;
     
   }

  // Initialize the object
  /*
   *  This is where we dimension everything.
   */
  bool LiquidTransport::initLiquid(LiquidTransportParams& tr) {

    int k;
    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(), 
	 m_thermo->molecularWeights().end(), m_mw.begin());

    /*
     *  Get the input Viscosities
     */
    m_viscSpecies.resize(m_nsp);
    m_viscTempDep_Ns.resize(m_nsp);
    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      m_viscTempDep_Ns[k] =  ltd.viscosity;
    }

    /*
     *  Get the input Thermal Conductivities
     */
    m_lambdaSpecies.resize(m_nsp);
    m_lambdaTempDep_Ns.resize(m_nsp);
    //for each species, assign thermal conductivity model 
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      m_lambdaTempDep_Ns[k] =  ltd.thermalCond;
    }

    /*
     *  Get the input Hydrodynamic Radii
     */
    m_hydrodynamic_radius.resize(m_nsp);
    m_radiusTempDep_Ns.resize(m_nsp);
    //for each species, assign model for hydrodynamic radius
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      m_radiusTempDep_Ns[k] =  ltd.hydroRadius;
    }

    /*
     *  Get the input Species Diffusivities
     *  Note that species diffusivities are not what is needed.
     *  Rather the Stefan Boltzmann interaction parameters are 
     *  needed for the current model.  This section may, therefore,
     *  be extraneous.
     */
    m_diffTempDep_Ns.resize(m_nsp);
    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      if ( ltd.speciesDiffusivity >= 0 ) {
	cout << "Warning: diffusion coefficient data for " 
	     << m_thermo->speciesName(k)
	     <<  endl 
	     << "in the input file is not used for LiquidTransport model."
	     <<  endl 
	     << "LiquidTransport model uses Stefan-Maxwell interaction "
	     <<  endl 
	     << "parameters defined in the <transport> input block." 
	     << endl;
      }
    }

    /*
     * Here we get interaction parameters from LiquidTransportParams 
     * that were filled in  TransportFactory::getLiquidInteractionsTransportData
     * Interaction models are provided here for viscosity, thermal conductivity,
     * species diffusivity and hydrodynamics radius (perhaps not needed in the 
     * present class).
     */
    m_viscMixModel = tr.viscosity;
    m_lambdaMixModel = tr.thermalCond;
    m_radiusMixModel = tr.hydroRadius;
    m_diffMixModel = tr.speciesDiffusivity;
    m_bdiff.resize(m_nsp,m_nsp);
    //Don't really need to update this here.  
    //It is updated in updateDiff_T()
    m_bdiff = m_diffMixModel->getMatrixTransProp(); 

    m_mode       = tr.mode_;

    m_molefracs.resize(m_nsp);
    m_molefracs_tran.resize(m_nsp);
    m_concentrations.resize(m_nsp);
    m_actCoeff.resize(m_nsp);
    m_chargeSpecies.resize(m_nsp);
    for ( int i = 0; i < m_nsp; i++ )
      m_chargeSpecies[i] = m_thermo->charge( i );
    m_volume_spec.resize(m_nsp);
    m_Grad_lnAC.resize(m_nsp); 
    m_spwork.resize(m_nsp);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);
    m_Grad_mu.resize(m_nDim * m_nsp, 0.0);

    m_flux.resize(m_nsp, m_nDim);
    m_Vdiff.resize(m_nsp, m_nDim);


    // set all flags to false
    m_visc_mix_ok   = false;
    m_visc_temp_ok  = false;
    m_visc_conc_ok  = false;
    m_radi_temp_ok  = false;
    m_radi_conc_ok  = false;
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
        
    update_T();
    update_C();

    if (m_visc_mix_ok) return m_viscmix;
  
    // update m_viscSpecies[] if necessary
    if (!m_visc_temp_ok) {
      updateViscosity_T();
    }

    if (!m_visc_conc_ok) {
      updateViscosities_C();
    }

    ////// LiquidTranInteraction method
    m_viscmix = m_viscMixModel->getMixTransProp( m_viscTempDep_Ns );

    return m_viscmix;
  }

  void LiquidTransport::getSpeciesViscosities(doublereal* visc) { 
    update_T();
    if (!m_visc_temp_ok) {
      updateViscosity_T();
    }
    copy(m_viscSpecies.begin(), m_viscSpecies.end(), visc); 
  }

  //===============================================================
  // Returns the hydrodynamic radius for all species
  /*
   *  The pure species viscosities are to be given in an Arrhenius
   * form in accordance with activated-jump-process dominated transport.
   */
  void LiquidTransport::getSpeciesHydrodynamicRadius(doublereal* const radius) {
    update_T();
    if (!m_radi_temp_ok) {
      updateHydrodynamicRadius_T();
    }
    copy(m_hydrodynamic_radius.begin(), m_hydrodynamic_radius.end(), radius); 

  }

 //================================================================

  /******************* binary diffusion coefficients **************/


  void LiquidTransport::getBinaryDiffCoeffs(int ld, doublereal* d) {
    int i,j;

    if ( ld != m_nsp ) 
      throw CanteraError("LiquidTransport::getBinaryDiffCoeffs",
			 "First argument does not correspond to number of species in model.\nDiff Coeff matrix may be misdimensioned");
    update_T();

    // if necessary, evaluate the binary diffusion coefficents
    // from the polynomial fits
    if (!m_diff_temp_ok) updateDiff_T();

    for (i = 0; i < m_nsp; i++) 
      for (j = 0; j < m_nsp; j++) {
	d[ld*j + i] = m_bdiff(i,j);

      }
  }


  //================================================================================================
  //  Get the electrical Mobilities (m^2/V/s).
  /*
   *   This function returns the mobilities. In some formulations
   *   this is equal to the normal mobility multiplied by faraday's constant.
   *
   *   Frequently, but not always, the mobility is calculated from the
   *   diffusion coefficient using the Einstein relation
   *
   *     \f[
   *          \mu^e_k = \frac{F D_k}{R T}
   *     \f]
   *
   * @param mobil_e  Returns the mobilities of
   *               the species in array \c mobil_e. The array must be
   *               dimensioned at least as large as the number of species.
   */
  void LiquidTransport::getMobilities(doublereal* const mobil) {
    int k;
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (k = 0; k < m_nsp; k++) {
      mobil[k] = c1 * m_spwork[k];
    }
  } 

  //================================================================================================
  //! Get the fluid mobilities (s kmol/kg).
  /*!
   *   This function returns the fluid mobilities. Usually, you have
   *   to multiply Faraday's constant into the resulting expression
   *   to general a species flux expression.
   *
   *   Frequently, but not always, the mobility is calculated from the
   *   diffusion coefficient using the Einstein relation
   *
   *     \f[ 
   *          \mu^f_k = \frac{D_k}{R T}
   *     \f]
   *
   *
   * @param mobil_f  Returns the mobilities of
   *               the species in array \c mobil. The array must be
   *               dimensioned at least as large as the number of species.
   */
  void  LiquidTransport::getFluidMobilities(doublereal* const mobil_f) {
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = 1.0 / (GasConstant * m_temp);
    for (int k = 0; k < m_nsp; k++) {
      mobil_f[k] = c1 * m_spwork[k];
    }
  } 
  //==============================================================
  void LiquidTransport::set_Grad_T(const doublereal* const grad_T) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_T[a] = grad_T[a];
    }
  }
  //==============================================================
  void LiquidTransport::set_Grad_V(const doublereal* const grad_V) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_V[a] = grad_V[a];
    }
  }
  //==============================================================
  void LiquidTransport::set_Grad_X(const doublereal* const grad_X) {
    int itop = m_nDim * m_nsp;
    for (int i = 0; i < itop; i++) {
      m_Grad_X[i] = grad_X[i];
    }
  }
  //==============================================================
  /****************** thermal conductivity **********************/  
  /*
   * The thermal conductivity is computed from the following mixture rule:
   *   \[
   *    \lambda = \left( \sum_k Y_k \lambda_k \right) 
   *   \]
   */
  doublereal LiquidTransport::thermalConductivity() {
   
    update_T();
    update_C();

    if (!m_cond_temp_ok) {
      updateCond_T();
    } 
    if (!m_cond_mix_ok) {

      // mass-fraction weighted thermal conductivity
      {
	doublereal sum1 = 0.0, sum2 = 0.0;
	for (int k = 0; k < m_nsp; k++) {
	  sum1 += m_molefracs[k] * m_mw[k] * m_lambdaSpecies[k];
	  sum2 += m_molefracs[k] * m_mw[k] ;
	}
	m_lambda = sum1 / sum2 ;
      }

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
  void LiquidTransport::getSpeciesFluxesES(int ndim, 
					   const doublereal* grad_T, 
					   int ldx, 
					   const doublereal* grad_X, 
					   int ldf, 
					   const doublereal* grad_V, 
					   doublereal* fluxes) {
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_V);
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

    update_T();
    update_C();

    update_Grad_lnAC();

    stefan_maxwell_solve();

    for (n = 0; n < m_nDim; n++) {
      for (k = 0; k < m_nsp; k++) {
	fluxes[n*ldf + k] = m_flux(k,n);
      }
    }

    /*
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
    for (n = 0; n < m_nDim; n++) 
      for (k = 0; k < m_nsp; k++) 
	fluxes[n*ldf + k] -= sum[n];
    */

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

    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_diff_temp_ok) {
      updateDiff_T();
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
   *  since the last call to update_T().
   *  If it hasn't then an immediate return is carried out.
   *
   *     @internal
   */ 
  bool LiquidTransport::update_T()
  {
    // First make a decision about whether we need to recalculate
    doublereal t = m_thermo->temperature();
    if (t == m_temp) return false;

    // Next do a reality check on temperature value
    if (t < 0.0) {
      throw CanteraError("LiquidTransport::update_T()",
			 "negative temperature "+fp2str(t));
    }

    // Compute various direct functions of temperature
    m_temp = t;
    m_logt = log(m_temp);
    m_kbt = Boltzmann * m_temp;

    // temperature has changed so temp flags are flipped
    m_visc_temp_ok  = false;
    m_radi_temp_ok  = false;
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
    return true;
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
  bool LiquidTransport::update_C() {
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
      return false;
    }

    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.
    m_visc_conc_ok = false;
  
    // Mixture stuff needs to be evaluated
    m_visc_mix_ok = false;
    m_diff_mix_ok = false;
    m_cond_mix_ok = false;

    return true;
  }


  // We formulate the directional derivative
  /*
   *     We only calculate the change in ac due to composition.
   *  The pressure and the temperature are taken care of in
   *  other parts of the expression.
   *
   */
  /*
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
  */

  //! Evaluate the gradient of the activity coefficients 
  //! as they alter the diffusion coefficient.  
  /**
   * The required quantity is the derivitive of the logarithm of the 
   * activity coefficient with respect to the derivative of the 
   * logarithm of the mole fraction (or whatever concentration 
   * variable we are using to express chemical potential.
   *
   * Returns the vector over species i:
   * \[
   *    1 + \partial \left[ \ln ( \gamma_i ) \right] 
   *       / \partial \left[ \ln ( \X_i  ) \right] 
   * \]
   */
  void LiquidTransport::update_Grad_lnAC() {

    int k;
    
    vector_fp grad_lnAC(m_nsp);
    m_thermo->getdlnActCoeffdlnC( DATA_PTR(grad_lnAC) ); 

    for (k = 0; k < m_nsp; k++) {
      m_Grad_lnAC[k] = grad_lnAC[k];
      //      std::cout << k << " m_Grad_lnAC = " << m_Grad_lnAC[k] << std::endl;
    }

    return;
  }

  /*************************************************************************
   *
   *    methods to update species temperature-dependent properties 
   *
   *************************************************************************/

  /**
   * Update the temperature-dependent parts of the species
   * thermal conductivity. 
   */
  void LiquidTransport::updateCond_T() {

    int k;

    for (k = 0; k < m_nsp; k++) {
      m_lambdaSpecies[k] = m_lambdaTempDep_Ns[k]->getSpeciesTransProp() ;
    }
    m_cond_temp_ok = true;
    m_cond_mix_ok = false;
  }


  //! Update the StefanMaxwell interaction parameters.  
  /**
   * These are evaluated using the Stokes-Einstein 
   * relation from the viscosity and hydrodynamic radius.
   */
  void LiquidTransport::updateDiff_T() {

    m_bdiff = m_diffMixModel->getMatrixTransProp();
    m_diff_temp_ok = true;
    m_diff_mix_ok = false;
  }


  //! Update the pure-species viscosities functional dependence on concentration.
  void LiquidTransport::updateViscosities_C() {
    m_visc_conc_ok = true;
  }


  /**
   * Update the temperature-dependent viscosity terms.
   * Updates the array of pure species viscosities, and the 
   * weighting functions in the viscosity mixture rule.
   * The flag m_visc_ok is set to true.
   *
   * Note that for viscosity, a positive activation energy 
   * corresponds to the typical case of a positive argument
   * to the exponential so that the Arrhenius expression is
   *
   * \f[
   *      \mu = A T^n \exp( + E / R T )
   * \f]
   */
  void LiquidTransport::updateViscosity_T() {
    int k;

    for (k = 0; k < m_nsp; k++) {
      m_viscSpecies[k] = m_viscTempDep_Ns[k]->getSpeciesTransProp() ;
    }
    m_visc_temp_ok = true;
    m_visc_mix_ok = false;
  }


  //! Update the pure-species viscosities functional dependence on concentration.
  void LiquidTransport::updateHydrodynamicRadius_C() {
    m_radi_conc_ok = true;
  }


  /**
   * Update the temperature-dependent hydrodynamic radius terms.
   * Updates the array of pure species viscosities, and the 
   * weighting functions in the viscosity mixture rule.
   * The flag m_visc_ok is set to true.
   */
  void LiquidTransport::updateHydrodynamicRadius_T() {
    int k;

    for (k = 0; k < m_nsp; k++) {
      m_hydrodynamic_radius[k] = m_radiusTempDep_Ns[k]->getSpeciesTransProp() ;
    } 
    m_radi_temp_ok = true;
    m_radi_mix_ok = false;
  }


  /*
   *
   *    Solve for the diffusional velocities in the Stefan-Maxwell equations
   *
   */
  void LiquidTransport::stefan_maxwell_solve() {
    int i, j, a;
    doublereal tmp;
    m_B.resize(m_nsp, m_nDim);
    m_A.resize(m_nsp, m_nsp);

    //! grab a local copy of the molecular weights
    const vector_fp& M =  m_thermo->molecularWeights();
     
    /*
     * Update the concentrations and diffusion coefficients in the mixture.
     */
    update_C();
    if ( !m_diff_temp_ok ) updateDiff_T();

    double T = m_thermo->temperature();

    update_Grad_lnAC() ;
 
    //m_thermo->getStandardVolumes(DATA_PTR(m_volume_spec));
    m_thermo->getActivityCoefficients(DATA_PTR(m_actCoeff));

    /* 
     *  Calculate the electrochemical potential gradient. This is the
     *  driving force for relative diffusional transport.
     *
     *  Here we calculate
     *
     *          X_i * (grad (mu_i) + S_i grad T - M_i / dens * grad P
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
      for (a = 0; a < m_nDim; a++) {
	m_Grad_mu[a*m_nsp + i] =
	  m_chargeSpecies[i] *  Faraday * m_Grad_V[a]
	  //+  (m_volume_spec[i] - M[i]/dens_) * m_Grad_P[a]
	  +  GasConstant * T * m_Grad_X[a*m_nsp+i] 
	     * ( 1.0 * m_Grad_lnAC[i] ) / xi_denom;
      }
    }

    if (m_thermo->activityConvention() == cAC_CONVENTION_MOLALITY) {
      int iSolvent = 0;
      double mwSolvent = m_thermo->molecularWeight(iSolvent);
      double mnaught = mwSolvent/ 1000.;
      double lnmnaught = log(mnaught);
      for (i = 1; i < m_nsp; i++) {
	for (a = 0; a < m_nDim; a++) {
	  m_Grad_mu[a*m_nsp + i] -=
	    m_molefracs[i] * GasConstant * m_Grad_T[a] * lnmnaught;
	}
      }
    }

    /*
     * Just for Note, m_A(i,j) refers to the ith row and jth column.
     * They are still fortran ordered, so that i varies fastest.
     */
    switch (m_nDim) {
    case 1:  /* 1-D approximation */
      m_B(0,0) = 0.0;
      for (j = 0; j < m_nsp; j++) {
	m_A(0,j) = m_molefracs_tran[j];
      }
      for (i = 1; i < m_nsp; i++){
	m_B(i,0) = m_Grad_mu[i] / (GasConstant * T);
	m_A(i,i) = 0.0;
	for (j = 0; j < m_nsp; j++){
	  if (j != i) {
	    if ( !( m_bdiff(i,j) > 0.0 ) )
		 throw CanteraError("LiquidTransport::stefan_maxwell_solve",
			     "m_bdiff has zero entry in non-diagonal.");
	    tmp = m_molefracs_tran[j] / m_bdiff(i,j);
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
	m_A(0,j) = m_molefracs_tran[j];
      }
      for (i = 1; i < m_nsp; i++){
	m_B(i,0) =  m_Grad_mu[i]         / (GasConstant * T);
	m_B(i,1) =  m_Grad_mu[m_nsp + i] / (GasConstant * T);
	m_A(i,i) = 0.0;
	for (j = 0; j < m_nsp; j++) {
	  if (j != i) {
	    if ( !( m_bdiff(i,j) > 0.0 ) )
		 throw CanteraError("LiquidTransport::stefan_maxwell_solve",
			     "m_bdiff has zero entry in non-diagonal.");
	    tmp =  m_molefracs_tran[j] / m_bdiff(i,j);
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
	m_A(0,j) = m_molefracs_tran[j];
      }
      for (i = 1; i < m_nsp; i++){
	m_B(i,0) = m_Grad_mu[i]           / (GasConstant * T);
	m_B(i,1) = m_Grad_mu[m_nsp + i]   / (GasConstant * T);
	m_B(i,2) = m_Grad_mu[2*m_nsp + i] / (GasConstant * T);
	m_A(i,i) = 0.0;
	for (j = 0; j < m_nsp; j++) {
	  if (j != i) {
	    if ( !( m_bdiff(i,j) > 0.0 ) )
		 throw CanteraError("LiquidTransport::stefan_maxwell_solve",
			     "m_bdiff has zero entry in non-diagonal.");
	    tmp =  m_molefracs_tran[j] / m_bdiff(i,j);
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

    for (a = 0; a < m_nDim; a++) {
      for (j = 0; j < m_nsp; j++) {
	m_Vdiff(j,a) = m_B(j,a);
	m_flux(j,a) = concTot_ * M[j] * m_molefracs_tran[j] * m_B(j,a);
      }
    }
  }


  /**
   * Throw an exception if this method is invoked. 
   * This probably indicates something is not yet implemented.
   */
  doublereal LiquidTransport::err(std::string msg) const {
    throw CanteraError("Liquid Transport Class",
		       "\n\n\n**** Method "+ msg +" not implemented in model "
		       + int2str(model()) + " ****\n"
		       "(Did you forget to specify a transport model?)\n\n\n");
      
    return 0.0;
  }


}
