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
    m_compositionDepType(-1),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_logt(0.0),
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
    m_compositionDepType(-1),
    m_iStateMF(-1),
    m_temp(-1.0),
    m_logt(0.0),
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
    m_viscTempDepType_Ns                  = right.m_viscTempDepType_Ns;
    m_lambdaTempDepType_Ns                = right.m_lambdaTempDepType_Ns;
    m_diffTempDepType_Ns                  = right.m_diffTempDepType_Ns;
    m_radiusTempDepType_Ns                = right.m_radiusTempDepType_Ns;
    m_coeffVisc_Ns                        = right.m_coeffVisc_Ns;
    m_coeffLambda_Ns                      = right.m_coeffLambda_Ns;
    m_coeffDiff_Ns                        = right.m_coeffDiff_Ns;
    m_coeffRadius_Ns                      = right.m_coeffRadius_Ns;
    m_visc_Eij                            = right.m_visc_Eij; 
    m_visc_Sij                            = right.m_visc_Sij; 
    m_hydrodynamic_radius                 = right.m_hydrodynamic_radius;
    m_Grad_X                              = right.m_Grad_X;
    m_Grad_T                              = right.m_Grad_T;
    m_Grad_V                              = right.m_Grad_V;
    m_ck_Grad_mu                          = right.m_ck_Grad_mu;
    m_bdiff                               = right.m_bdiff;
    m_viscSpecies                         = right.m_viscSpecies;
    m_logViscSpecies                      = right.m_logViscSpecies;
    m_lambdaSpecies                       = right.m_lambdaSpecies;
    m_iStateMF = -1;
    m_molefracs                           = right.m_molefracs;
    m_concentrations                      = right.m_concentrations;
    m_chargeSpecies                       = right.m_chargeSpecies;
    m_DiffCoeff_StefMax                   = right.m_DiffCoeff_StefMax;
    viscosityModel_                       = right.viscosityModel_;
    m_B                                   = right.m_B;
    m_A                                   = right.m_A;
    m_temp                                = right.m_temp;
    m_logt                                = right.m_logt;
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

    int k;
    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    /*
     * Read the transport block in the phase XML Node
     * It's not an error if this block doesn't exist. Just use the defaults
     */
    XML_Node &phaseNode = m_thermo->xml();
    if (phaseNode.hasChild("transport")) {
      XML_Node& transportNode = phaseNode.child("transport");
      if ( transportNode.hasChild("viscosity")) {
	XML_Node& viscosityNode = transportNode.child("viscosity");
	string viscosityModel = viscosityNode.attrib("model");
	if (viscosityModel == "") {
	  throw CanteraError("LiquidTransport::initLiquid",
			     "transport::visosity XML node doesn't have a model string");
	}
      }



      string transportModel = transportNode.attrib("model");
      if (transportModel == "LiquidTransport") {
        /*
         * <compositionDependence model="Solvent_Only"/>
	 *      or
	 * <compositionDependence model="Mixture_Averaged"/>
	 */
	std::string modelName = "";
	if (getOptionalModel(transportNode, "compositionDependence",
			     modelName)) {
	  modelName = lowercase(modelName);
          if (modelName == "solvent_only") {
	    m_compositionDepType = 0;
	  } else if (modelName == "mixture_averaged") {
	    m_compositionDepType = 1;
	  } else {
	    throw CanteraError("LiquidTransport::initLiquid", "Unknown compositionDependence Model: " + modelName);
	  }
	}

      


      }
    }

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(), 
	 m_thermo->molecularWeights().end(), m_mw.begin());

    /*
     *  Get the input Viscosities
     */
    m_viscSpecies.resize(m_nsp);
    m_coeffVisc_Ns.clear(); 
    m_coeffVisc_Ns.resize(m_nsp);
    m_viscTempDepType_Ns.resize(m_nsp);

    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      //specify temperature dependence
      m_viscTempDepType_Ns[k] =  ltd.model_viscosity;
      //vector kentry corresponds to the k-th entry of m_coeffVisc_Ns
      vector_fp &kentry = m_coeffVisc_Ns[k]; 

      if ( m_viscTempDepType_Ns[k] == LTR_MODEL_CONSTANT
	   || m_viscTempDepType_Ns[k] == LTR_MODEL_POLY ) {
	kentry = ltd.viscCoeffs;

      } else if ( m_viscTempDepType_Ns[k] == LTR_MODEL_ARRHENIUS ) {
	kentry = ltd.viscCoeffs;
	//for Arrhenius form, also carry the logarithm of the pre-exponential
	kentry[3] = log( kentry[0] );

      } else if ( m_viscTempDepType_Ns[k] == LTR_MODEL_NOTSET ) {
	//we might be OK with viscosity not being set so
	// this error is repeated in updateViscosity_T()
	// and can be deleted from here if appropriate
	throw CanteraError("LiquidTransport::initLiquid",
			   "Viscosity Model is not set for species " 
			   + m_thermo->speciesName(k) 
			   + " in the input file");
      } else {
	throw CanteraError("LiquidTransport::initLiquid",
			   "Viscosity Model for species "
			   + m_thermo->speciesName(k)
			   + " is not handled by this object");
      }
    }

    /*
     *  Get the input Thermal Conductivities
     */
    m_lambdaSpecies.resize(m_nsp);
    m_coeffLambda_Ns.clear(); 
    m_coeffLambda_Ns.resize(m_nsp);
    m_lambdaTempDepType_Ns.resize(m_nsp);

    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      //specify temperature dependence
      m_lambdaTempDepType_Ns[k] =  ltd.model_thermalCond;
      //vector kentry corresponds to the k-th entry of m_coeffLambda_Ns
      vector_fp &kentry = m_coeffLambda_Ns[k]; 

      if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_CONSTANT
	   || m_lambdaTempDepType_Ns[k] == LTR_MODEL_POLY ) {
	kentry = ltd.thermalCondCoeffs;

      } else if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_ARRHENIUS ) {
	kentry = ltd.thermalCondCoeffs;
	//for Arrhenius form, also carry the logarithm of the pre-exponential
	kentry[3] = log( kentry[0] );

      } else if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_NOTSET ) {
	throw CanteraError("LiquidTransport::initLiquid",
			   "Thermal conductivity model is not set for species " 
			   + m_thermo->speciesName(k)
			   + " in the input file");
      } else {
	throw CanteraError("LiquidTransport::initLiquid",
			   "Thermal conductivity model for species "
			   + m_thermo->speciesName(k)
			   + " is not handled by this object");
      }
    }

    /*
     *  Get the input Hydrodynamic Radii
     */
    m_hydrodynamic_radius.resize(m_nsp);
    m_coeffRadius_Ns.clear(); 
    m_coeffRadius_Ns.resize(m_nsp);
    m_radiusTempDepType_Ns.resize(m_nsp);

    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      //specify temperature dependence
      m_radiusTempDepType_Ns[k] =  ltd.model_hydroradius;
      //vector kentry corresponds to the k-th entry of m_coeffRadius_Ns
      vector_fp &kentry = m_coeffRadius_Ns[k]; 

      if ( m_radiusTempDepType_Ns[k] == LTR_MODEL_CONSTANT
	   || m_radiusTempDepType_Ns[k] == LTR_MODEL_POLY ) {
	kentry = ltd.hydroRadiusCoeffs;

      } else if ( m_radiusTempDepType_Ns[k] == LTR_MODEL_ARRHENIUS ) {
	kentry = ltd.hydroRadiusCoeffs;
	//for Arrhenius form, also carry the logarithm of the pre-exponential
	kentry[3] = log( kentry[0] );

      } else if ( m_radiusTempDepType_Ns[k] == LTR_MODEL_NOTSET ) {
	throw CanteraError("LiquidTransport::initLiquid",
			   "Hydrodynamic radius model is not set for species " 
			   + m_thermo->speciesName(k)
			   + " in the input file");
      } else {
	throw CanteraError("LiquidTransport::initLiquid",
			   "Hydrodynamic radius model for species "
			   + m_thermo->speciesName(k)
			   + " is not handled by this object");
      }
    }

    /*
     *  Get the input Species Diffusivities
     *  Note that species diffusivities are not what is needed.
     *  Rather the Stefan Boltzmann interaction parameters are 
     *  needed for the current model.  This section may, therefore,
     *  be extraneous.
     */
    //    m_viscSpecies.resize(m_nsp);
    m_coeffDiff_Ns.clear(); 
    m_coeffDiff_Ns.resize(m_nsp);
    m_diffTempDepType_Ns.resize(m_nsp);

    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      //specify temperature dependence
      if ( ltd.model_speciesDiffusivity >= 0 
	   || ltd.speciesDiffusivityCoeffs.size() > 0 ) {
	cout << "Warning: diffusion coefficient data for " 
	     << m_thermo->speciesName(k)
	     <<  endl 
	     << "in the input file is not used for LiquidTransport model."
	     <<  endl 
	     << "LiquidTransport model uses hydrodynamicRadius, viscosity "
	     <<  endl 
	     << "and the Stokes-Einstein equation." 
	     << endl;
      }
    }

    m_mode       = tr.mode_;

    m_viscSpecies.resize(m_nsp);
    m_logViscSpecies.resize(m_nsp);
    m_lambdaSpecies.resize(m_nsp);
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

    /* We still need to implement interaction parameters */
    /* This constant viscosity model has no input */

    if (viscosityModel_ == LVISC_CONSTANT) {

      err("constant viscosity not implemented for LiquidTransport.");
      //return m_viscmix;

    } else if (viscosityModel_ == LVISC_AVG_ENERGIES) {

      m_viscmix = exp( dot_product(m_logViscSpecies, m_molefracs) );

    } else if (viscosityModel_ == LVISC_INTERACTION) {

      // log_visc_mix = sum_i (X_i log_visc_i) + sum_i sum_j X_i X_j G_ij
      double interaction = dot_product(m_logViscSpecies, m_molefracs);
      for ( int i = 0; i < m_nsp; i++ ) 
	for ( int j = 0; j < i; j++ ) 
	  interaction += m_molefracs[i] * m_molefracs[j] 
	    * ( m_visc_Sij(i,j) + m_visc_Eij(i,j) / m_temp );
      m_viscmix = exp( interaction );

    } else {
      err("Unknown viscosity model in LiquidTransport::viscosity().");
    }
    
    return m_viscmix;
  }

  void LiquidTransport::getSpeciesViscosities(doublereal* visc) { 
    update_T();
    if (!m_visc_temp_ok) {
      updateViscosity_T();
    }
    copy(m_viscSpecies.begin(), m_viscSpecies.end(), visc); 
  }
  //====================================================================================================================
  // Returns the hydrodynamic radius for all species
  /*
   *  The pure species viscosities are to be given in an Arrhenius
   * form in accordance with activated-jump-process dominated transport.
   */
  void LiquidTransport::getSpeciesHydrodynamicRadius(doublereal* const radius) {
  }
 //====================================================================================================================

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
  //================================================================================================
  void LiquidTransport::set_Grad_V(const doublereal* const grad_V) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_V[a] = grad_V[a];
    }
  }
  //================================================================================================
  void LiquidTransport::set_Grad_T(const doublereal* const grad_T) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_T[a] = grad_T[a];
    }
  }
  //================================================================================================
  void LiquidTransport::set_Grad_X(const doublereal* const grad_X) {
    int itop = m_nDim * m_nsp;
    for (int i = 0; i < itop; i++) {
      m_Grad_X[i] = grad_X[i];
    }
    update_Grad_lnAC();
  }
  //================================================================================================
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
  void LiquidTransport::getSpeciesFluxesExt(int ldf, doublereal* fluxes) {
    int n, k;

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
      vector_fp &coeffk = m_coeffLambda_Ns[k];

      if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_CONSTANT ) {
	m_lambdaSpecies[k] = coeffk[0] ;

      } else if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_ARRHENIUS ) {
	//m_coeffLambda_Ns[k][0] holds A
	//m_coeffLambda_Ns[k][1] holds n
	//m_coeffLambda_Ns[k][2] holds Tact
	//m_coeffLambda_Ns[k][3] holds log(A)
	m_lambdaSpecies[k] = coeffk[0] * exp( coeffk[1] * m_logt 
					      - coeffk[2] / m_temp );
      
      } else if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_POLY ) {
	m_lambdaSpecies[k] = coeffk[0]
	  + coeffk[1] * m_temp
	  + coeffk[2] * m_temp * m_temp
	  + coeffk[3] * m_temp * m_temp * m_temp
	  + coeffk[4] * m_temp * m_temp * m_temp * m_temp;

      } else if ( m_lambdaTempDepType_Ns[k] == LTR_MODEL_NOTSET ) {
	throw CanteraError("LiquidTransport::updateCond_T",
			   "Conductivity Model is not set for species " 
			   + m_thermo->speciesName(k) 
			   + " in the input file");
      } else {
	throw CanteraError("LiquidTransport::updateCond_T",
			   "Conductivity Model for species "
			   + m_thermo->speciesName(k) 
			   + " is not handled by this object");
      }
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

    double *viscSpec = new double(m_nsp);
    double *radiusSpec = new double(m_nsp);
    getSpeciesViscosities( viscSpec );
    getSpeciesHydrodynamicRadius( radiusSpec );

    int i,j;
    for (i = 0; i < m_nsp; i++) 
      for (j = 0; j < m_nsp; j++) {
	m_DiffCoeff_StefMax(i,j) = m_bdiff(i,j) = GasConstant * m_temp 
	  / ( 6.0 * Pi * radiusSpec[i] * viscSpec[j] ) ;
	cout << " D_ij = " << m_bdiff(i,j) << " for " 
	     << m_thermo->speciesName(i) << ", " 
	     << m_thermo->speciesName(j) << endl; 
      }
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
   */
  void LiquidTransport::updateViscosity_T() {
    int k;

    for (k = 0; k < m_nsp; k++) {
      vector_fp &coeffk = m_coeffVisc_Ns[k];

      if ( m_viscTempDepType_Ns[k] == LTR_MODEL_CONSTANT ) {
	m_logViscSpecies[k] = log( coeffk[0] );
	m_viscSpecies[k] = coeffk[0] ;

      } else if ( m_viscTempDepType_Ns[k] == LTR_MODEL_ARRHENIUS ) {
	//m_coeffVisc_Ns[k][0] holds A
	//m_coeffVisc_Ns[k][1] holds n
	//m_coeffVisc_Ns[k][2] holds Tact
	//m_coeffVisc_Ns[k][3] holds log(A)
	m_logViscSpecies[k] = coeffk[3] + coeffk[1] * m_logt 
	  - coeffk[2] / m_temp ;
	m_viscSpecies[k] = exp( m_logViscSpecies[k] );
      
      } else if ( m_viscTempDepType_Ns[k] == LTR_MODEL_POLY ) {
	m_viscSpecies[k] = coeffk[0]
	  + coeffk[1] * m_temp
	  + coeffk[2] * m_temp * m_temp
	  + coeffk[3] * m_temp * m_temp * m_temp
	  + coeffk[4] * m_temp * m_temp * m_temp * m_temp;
	m_logViscSpecies[k] = log( m_viscSpecies[k] );

      } else if ( m_viscTempDepType_Ns[k] == LTR_MODEL_NOTSET ) {
	throw CanteraError("LiquidTransport::updateViscosity_T",
			   "Viscosity Model is not set for species " 
			   + m_thermo->speciesName(k) 
			   + " in the input file");
      } else {
	throw CanteraError("LiquidTransport::updateViscosity_T",
			   "Viscosity Model for species "
			   + m_thermo->speciesName(k) 
			   + " is not handled by this object");
      }
      m_visc_temp_ok = true;
      m_visc_mix_ok = false;
    }
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
     * Update the concentrations and diffusion coefficients in the mixture.
     */
    update_C();
    if ( !m_diff_temp_ok ) updateDiff_T();

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
