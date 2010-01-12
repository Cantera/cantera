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
    m_massfracs                           = right.m_massfracs;
    m_massfracs_tran                      = right.m_massfracs_tran;
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
     //if ( m_radiusMixModel ) delete m_radiusMixModel;
     
   }

  // Initialize the transport object
  /*
   * Here we change all of the internal dimensions to be sufficient.
   * We get the object ready to do property evaluations.
   * A lot of the input required to do property evaluations is 
   * contained in the LiquidTransportParams class that is 
   * filled in TransportFactory. 
   *
   * @param tr  Transport parameters for all of the species
   *            in the phase.
   */
  bool LiquidTransport::initLiquid(LiquidTransportParams& tr) {

    int k;
    // constant substance attributes
    m_thermo = tr.thermo;
    m_velocityBasis = tr.velocityBasis_;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp, 0.0);
    copy(m_thermo->molecularWeights().begin(), 
	 m_thermo->molecularWeights().end(), m_mw.begin());

    /*
     *  Get the input Viscosities
     */
    m_viscSpecies.resize(m_nsp, 0.0);
    m_viscTempDep_Ns.resize(m_nsp, 0);
    //for each species, assign viscosity model and coefficients
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      m_viscTempDep_Ns[k] =  ltd.viscosity;
    }

    /*
     *  Get the input Thermal Conductivities
     */
    m_lambdaSpecies.resize(m_nsp, 0.0);
    m_lambdaTempDep_Ns.resize(m_nsp, 0);
    //for each species, assign thermal conductivity model 
    for (k = 0; k < m_nsp; k++) {
      Cantera::LiquidTransportData &ltd = tr.LTData[k];
      m_lambdaTempDep_Ns[k] =  ltd.thermalCond;
    }

    /*
     *  Get the input Hydrodynamic Radii
     */
    m_hydrodynamic_radius.resize(m_nsp, 0.0);
    m_radiusTempDep_Ns.resize(m_nsp, 0);
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
    m_diffTempDep_Ns.resize(m_nsp, 0);
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
    //m_radiusMixModel = tr.hydroRadius;
    m_diffMixModel = tr.speciesDiffusivity;
    m_bdiff.resize(m_nsp,m_nsp, 0.0);
    //Don't really need to update this here.  
    //It is updated in updateDiff_T()
    m_diffMixModel->getMatrixTransProp( m_bdiff ); 

    m_mode       = tr.mode_;

    m_massfracs.resize(m_nsp, 0.0);
    m_massfracs_tran.resize(m_nsp, 0.0);
    m_molefracs.resize(m_nsp, 0.0);
    m_molefracs_tran.resize(m_nsp, 0.0);
    m_concentrations.resize(m_nsp, 0.0);
    m_actCoeff.resize(m_nsp, 0.0);
    m_chargeSpecies.resize(m_nsp, 0.0);
    for ( int i = 0; i < m_nsp; i++ )
      m_chargeSpecies[i] = m_thermo->charge( i );
    m_volume_spec.resize(m_nsp, 0.0);
    m_Grad_lnAC.resize(m_nsp, 0.0); 
    m_spwork.resize(m_nsp, 0.0);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);
    m_Grad_mu.resize(m_nDim * m_nsp, 0.0);

    m_flux.resize(m_nsp, m_nDim, 0.0);
    m_Vdiff.resize(m_nsp, m_nDim, 0.0);


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

  // Returns the viscosity of the solution
  /*
   *  The viscosity calculation is handled by subclasses of 
   *  LiquidTranInteraction as specified in the input file.  
   *  These in turn employ subclasses of LTPspecies to 
   *  determine the individual species viscosities.
   */ 
  doublereal LiquidTransport::viscosity() {
        
    update_T();
    update_C();

    if (m_visc_mix_ok) return m_viscmix;
  
    ////// LiquidTranInteraction method
    m_viscmix = m_viscMixModel->getMixTransProp( m_viscTempDep_Ns );

    return m_viscmix;

    /*
    // update m_viscSpecies[] if necessary
    if (!m_visc_temp_ok) {
      updateViscosity_T();
    }

    if (!m_visc_conc_ok) {
      updateViscosities_C();
    }
    */
  }

  // Returns the pure species viscosities for all species
  /*
   *  The pure species viscosities are evaluated using the 
   *  appropriate subclasses of LTPspecies as specified in the 
   *  input file.
   *
   * @param visc  array of length "number of species"
   *              to hold returned viscosities.
   */
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
   *  The species hydrodynamic radii are evaluated using the 
   *  appropriate subclasses of LTPspecies as specified in the 
   *  input file.
   *
   * @param radius  array of length "number of species"
   *                to hold returned radii.
   */
  void LiquidTransport::getSpeciesHydrodynamicRadius(doublereal* const radius) {
    update_T();
    if (!m_radi_temp_ok) {
      updateHydrodynamicRadius_T();
    }
    copy(m_hydrodynamic_radius.begin(), m_hydrodynamic_radius.end(), radius); 

  }

 //================================================================

  // Return the thermal conductivity of the solution
  /*
   *  The thermal conductivity calculation is handled by subclasses of 
   *  LiquidTranInteraction as specified in the input file.  
   *  These in turn employ subclasses of LTPspecies to 
   *  determine the individual species thermal condictivities.
   */ 
  doublereal LiquidTransport::thermalConductivity() {
   
    update_T();
    update_C();

    if (!m_cond_mix_ok) {
     m_lambda = m_lambdaMixModel->getMixTransProp( m_lambdaTempDep_Ns );
     m_cond_mix_ok = true;
    } 

    return m_lambda;
  }


  /****************** thermal diffusion coefficients ************/

  //! Return the thermal diffusion coefficients
  /*!
   *  These are all zero for this simple implementaion
   *
   *  @param dt thermal diffusion coefficients
   */
  void LiquidTransport::getThermalDiffCoeffs(doublereal* const dt) {
    for (int k = 0; k < m_nsp; k++) {
      dt[k] = 0.0;
    }
  }

  /******************* binary diffusion coefficients **************/


  // Returns the binary diffusion coefficients
  /*
   *   The binary diffusion coefficients are specified in the input
   *   file through the LiquidTransportInteractions class.  These
   *   are the binary interaction coefficients employed in the 
   *   Stefan-Maxwell equation.
   *   
   *   @param ld  number of species in system
   *   @param d   vector of binary diffusion coefficients
   *          units = m2 s-1. length = ld*ld = (number of species)^2
   */
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
  // Get the Electrical mobilities (m^2/V/s).
  /*
   *  The electrical mobilities are not well defined 
   *  in the context of LiquidTransport because the Stefan Maxwell 
   *  equation is solved.  Here the electrical mobilities 
   *  are calculated from the mixture-averaged
   *  diffusion coefficients through a call to getMixDiffCoeffs() 
   *  using the Einstein relation
   *
   *     \f[ 
   *          \mu^e_k = \frac{F D_k}{R T}
   *     \f]
   *
   *  Note that this call to getMixDiffCoeffs() requires
   *  a solve of the Stefan Maxwell equation making this  
   *  determination of the mixture averaged diffusion coefficients 
   *  a {\em slow} method for obtaining diffusion coefficients.  
   * 
   *  Also note that the Stefan Maxwell solve will be based upon 
   *  the thermodynamic state (including gradients) most recently
   *  set.  Gradients can be set specifically using set_Grad_V,
   *  set_Grad_X and set_Grad_T or through calls to 
   *  getSpeciesFluxes, getSpeciesFluxesES, getSpeciesVdiff, 
   *  getSpeciesVdiffES, etc.
   *
   * @param mobil_e  Returns the electrical mobilities of
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
  // Get the fluid mobilities (s kmol/kg).
  /*
   *  The fluid mobilities are not well defined 
   *  in the context of LiquidTransport because the Stefan Maxwell 
   *  equation is solved.  Here the fluid mobilities 
   *  are calculated from the mixture-averaged
   *  diffusion coefficients through a call to getMixDiffCoeffs() 
   *  using the Einstein relation
   *
   *     \f[ 
   *          \mu^f_k = \frac{D_k}{R T}
   *     \f]
   *
   *  Note that this call to getMixDiffCoeffs() requires
   *  a solve of the Stefan Maxwell equation making this  
   *  determination of the mixture averaged diffusion coefficients 
   *  a {\em slow} method for obtaining diffusion coefficients.  
   * 
   *  Also note that the Stefan Maxwell solve will be based upon 
   *  the thermodynamic state (including gradients) most recently
   *  set.  Gradients can be set specifically using set_Grad_V,
   *  set_Grad_X and set_Grad_T or through calls to 
   *  getSpeciesFluxes, getSpeciesFluxesES, getSpeciesVdiff, 
   *  getSpeciesVdiffES, etc.
   *
   * @param mobil_f  Returns the fluid mobilities of
   *               the species in array \c mobil_f. The array must be
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
  //! Specify the value of the gradient of the temperature
  /*!
   * @param grad_T Gradient of the temperature (length num dimensions);
   */
  void LiquidTransport::set_Grad_T(const doublereal* const grad_T) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_T[a] = grad_T[a];
    }
  }
  //==============================================================
  //! Specify the value of the gradient of the voltage
  /*!
   *
   * @param grad_V Gradient of the voltage (length num dimensions);
   */
  void LiquidTransport::set_Grad_V(const doublereal* const grad_V) {
    for (int a = 0; a < m_nDim; a++) {
      m_Grad_V[a] = grad_V[a];
    }
  }
  //==============================================================
  //! Specify the value of the gradient of the MoleFractions
  /*!
   *
   * @param grad_X Gradient of the mole fractions(length nsp * num dimensions);
   */
  void LiquidTransport::set_Grad_X(const doublereal* const grad_X) {
    int itop = m_nDim * m_nsp;
    for (int i = 0; i < itop; i++) {
      m_Grad_X[i] = grad_X[i];
    }
  }
  //==============================================================

  // Compute the mixture electrical conductivity from 
  // the Stefan-Maxwell equation.
  /*
   *  To compute the mixture electrical conductance, the Stefan
   *  Maxwell equation is solved for zero species gradients and 
   *  for unit potential gradient, \f$ \nabla V \f$.  
   *  The species fluxes are converted to current by summing over 
   *  the charge-weighted fluxes according to 
   *  \f[
   *      \vec{i} = \sum_{i} z_i F \rho \vec{V_i} / W_i 
   *  \f]
   *  where \f$ z_i \f$ is the charge on species i,
   *  \f$ F \f$ is Faradays constant,  \f$ \rho \f$  is the density,
   *  \f$ W_i \f$ is the molecular mass of species i.
   *  The conductance, \f$ \kappa \f$ is obtained from 
   *  \f[
   *      \kappa = \vec{i} / \nabla V.
   *  \f]
   */
  doublereal LiquidTransport::getElectricConduct( ) {
    doublereal gradT = 0.0;
    doublereal gradX[m_nDim * m_nsp];
    doublereal gradV[m_nDim];
    for (int i = 0; i < m_nDim; i++) {
      for (int k = 0; k < m_nsp; k++) 
	gradX[ i*m_nDim + k] = 0.0;
      gradV[i] = 1.0;
    }

    set_Grad_T(&gradT);
    set_Grad_X(gradX);
    set_Grad_V(gradV);

    doublereal *fluxes = new doublereal( m_nsp * m_nDim );
    doublereal current;

    getSpeciesFluxesExt(m_nDim, fluxes);

    //sum over species charges, fluxes, Faraday to get current
    // Since we want the scalar conductivity, we need only consider one-dim
    for (int i = 0; i < 1; i++) {
      current = 0.0;
      for (int k = 0; k < m_nsp; k++) 
	current += m_chargeSpecies[k] *  Faraday  * fluxes[k] / m_mw[k];
      //divide by unit potential gradient
      current /= - gradV[i];
    }
    delete fluxes;
    return current;

  }

  //! Compute the electric current density in A/m^2
  /**
   *  The electric current is computed first by computing the 
   *  species diffusive fluxes using  the Stefan Maxwell solution
   *  and then the current, \f$ \vec{i} \f$ by summing over 
   *  the charge-weighted fluxes according to 
   *  \f[
   *      \vec{i} = \sum_{i} z_i F \rho \vec{V_i} / W_i 
   *  \f]
   *  where \f$ z_i \f$ is the charge on species i,
   *  \f$ F \f$ is Faradays constant,  \f$ \rho \f$  is the density,
   *  \f$ W_i \f$ is the molecular mass of species i.
   * 
   * @param ndim The number of spatial dimensions (1, 2, or 3).
   * @param grad_T The temperature gradient (ignored in this model).
   * @param ldx  Leading dimension of the grad_X array.
   * @param grad_T The temperature gradient (ignored in this model).
   * @param ldf  Leading dimension of the grad_V and current vectors.
   * @param grad_V The electrostatic potential gradient.
   * @param current The electric current in A/m^2.
   */
  void LiquidTransport::getElectricCurrent(int ndim, 
					  const doublereal* grad_T, 
					  int ldx, 
					  const doublereal* grad_X, 
					  int ldf, 
					  const doublereal* grad_V, 
					  doublereal* current) {
    
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_V);

    doublereal *fluxes = new doublereal( m_nsp * m_nDim );

    getSpeciesFluxesExt(ldf, fluxes);

    //sum over species charges, fluxes, Faraday to get current
    for (int i = 0; i < m_nDim; i++) {
      current[i] = 0.0;
      for (int k = 0; k < m_nsp; k++) 
	current[i] += m_chargeSpecies[k] * Faraday * fluxes[k] / m_mw[k];
      //divide by unit potential gradient
    }
    delete fluxes;

  }

  // Get the species diffusive velocities wrt to 
  // the averaged velocity, 
  // given the gradients in mole fraction and temperature
  /* 
   * The average velocity can be computed on a mole-weighted 
   * or mass-weighted basis, or the diffusion velocities may 
   * be specified as relative to a specific species (i.e. a 
   * solvent) all according to the velocityBasis input parameter.
   * 
   *  Units for the returned fluxes are kg m-2 s-1.
   * 
   *  @param ndim Number of dimensions in the flux expressions
   *  @param grad_T Gradient of the temperature
   *                 (length = ndim)
   *  @param ldx  Leading dimension of the grad_X array 
   *              (usually equal to m_nsp but not always)
   *  @param grad_X Gradients of the mole fraction
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
   *  @param ldf  Leading dimension of the fluxes array 
   *              (usually equal to m_nsp but not always)
   *  @param Vdiff  Output of the diffusive velocities.
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
   */
  void LiquidTransport::getSpeciesVdiff(int ndim, 
					 const doublereal* grad_T, 
					 int ldx, const doublereal* grad_X, 
					 int ldf, doublereal* Vdiff) {
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesVdiffExt(ldf, Vdiff);
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
  void LiquidTransport::getSpeciesVdiffES(int ndim, 
					   const doublereal* grad_T, 
					   int ldx, 
					   const doublereal* grad_X, 
					   int ldf, 
					   const doublereal* grad_V, 
					   doublereal* Vdiff) {
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    set_Grad_V(grad_V);
    getSpeciesVdiffExt(ldf, Vdiff);
  }

  //  Return the species diffusive mass fluxes wrt to
  //  the averaged velocity in [kmol/m^2/s].
  /*
   *
   * The diffusive mass flux of species \e k is computed 
   * using the Stefan-Maxwell equation
   * \f[
   *     X_i \nabla \mu_i 
   *                     = RT \sum_i \frac{X_i X_j}{D_{ij}} 
   *                           ( \vec{V}_j - \vec{V}_i )
   * \f]
   * to determine the diffusion velocity and 
   * \f[
   *      \vec{N}_i = C_T X_i \vec{V}_i
   * \f]
   * to determine the diffusion flux.  Here \f$ C_T \f$ is the 
   * total concentration of the mixture [kmol/m^3], \f$ D_{ij} \f$
   * are the Stefa-Maxwell interaction parameters in [m^2/s],
   * \f$ \vec{V}_{i} \f$ is the diffusion velocity of species \e i,
   * \f$ \mu_i \f$ is the electrochemical potential of species \e i.
   *
   * Note that for this method, there is no argument for the 
   * gradient of the electric potential (voltage).  Electric 
   * potential gradients can be set with set_Grad_V() or
   * method getSpeciesFluxesES() can be called.x
   *
   * The diffusion velocity is relative to an average velocity 
   * that can be computed on a mole-weighted 
   * or mass-weighted basis, or the diffusion velocities may 
   * be specified as relative to a specific species (i.e. a 
   * solvent) all according to the \verbatim <velocityBasis> 
   * \endverbatim input parameter.
   
   * @param ndim The number of spatial dimensions (1, 2, or 3).
   * @param grad_T The temperature gradient (ignored in this model).
   *                 (length = ndim)
   * @param ldx  Leading dimension of the grad_X array.
   *              (usually equal to m_nsp but not always)
   * @param grad_X Gradients of the mole fraction
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
   * @param ldf  Leading dimension of the fluxes array 
   *              (usually equal to m_nsp but not always)
   * @param grad_Phi Gradients of the electrostatic potential
   *             length = ndim
   * @param fluxes  Output of the diffusive mass fluxes
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
   */
  void LiquidTransport::getSpeciesFluxes(int ndim, 
					 const doublereal* grad_T, 
					 int ldx, const doublereal* grad_X, 
					 int ldf, doublereal* fluxes) {
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesFluxesExt(ldf, fluxes);
  }

  //  Return the species diffusive mass fluxes wrt to
  //  the averaged velocity in [kmol/m^2/s].
  /*
   *
   * The diffusive mass flux of species \e k is computed 
   * using the Stefan-Maxwell equation
   * \f[
   *     X_i \nabla \mu_i 
   *                     = RT \sum_i \frac{X_i X_j}{D_{ij}} 
   *                           ( \vec{V}_j - \vec{V}_i )
   * \f]
   * to determine the diffusion velocity and 
   * \f[
   *      \vec{N}_i = C_T X_i \vec{V}_i
   * \f]
   * to determine the diffusion flux.  Here \f$ C_T \f$ is the 
   * total concentration of the mixture [kmol/m^3], \f$ D_{ij} \f$
   * are the Stefa-Maxwell interaction parameters in [m^2/s],
   * \f$ \vec{V}_{i} \f$ is the diffusion velocity of species \e i,
   * \f$ \mu_i \f$ is the electrochemical potential of species \e i.
   *
   * The diffusion velocity is relative to an average velocity 
   * that can be computed on a mole-weighted 
   * or mass-weighted basis, or the diffusion velocities may 
   * be specified as relative to a specific species (i.e. a 
   * solvent) all according to the \verbatim <velocityBasis> 
   * \endverbatim input parameter.
   
   * @param ndim The number of spatial dimensions (1, 2, or 3).
   * @param grad_T The temperature gradient (ignored in this model).
   *                 (length = ndim)
   * @param ldx  Leading dimension of the grad_X array.
   *              (usually equal to m_nsp but not always)
   * @param grad_X Gradients of the mole fraction
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
   * @param ldf  Leading dimension of the fluxes array 
   *              (usually equal to m_nsp but not always)
   * @param grad_Phi Gradients of the electrostatic potential
   *             length = ndim
   * @param fluxes  Output of the diffusive mass fluxes
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
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

  //  Return the species diffusive velocities relative to 
  //  the averaged velocity.  
  /*
   * This method acts similarly to getSpeciesVdiffES() but
   * requires all gradients to be preset using methods 
   * set_Grad_X(), set_Grad_V(), set_Grad_T().  
   * See the documentation of getSpeciesVdiffES() for details.
   *      
   *  @param ldf  Leading dimension of the Vdiff array.
   *  @param Vdiff  Output of the diffusive velocities.
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
   */
  void LiquidTransport::getSpeciesVdiffExt(int ldf, doublereal* Vdiff) {
    int n, k;

    update_T();
    update_C();

    update_Grad_lnAC();

    stefan_maxwell_solve();

    for (n = 0; n < m_nDim; n++) {
      for (k = 0; k < m_nsp; k++) {
	Vdiff[n*ldf + k] = m_Vdiff(k,n);
      }
    }
  }

  //  Return the species diffusive fluxes relative to 
  //  the averaged velocity.  
  /*
   * This method acts similarly to getSpeciesFluxesES() but
   * requires all gradients to be preset using methods 
   * set_Grad_X(), set_Grad_V(), set_Grad_T().  
   * See the documentation of getSpeciesFluxesES() for details.
   *      
   *  units = kg/m2/s
   *
   *  @param ldf  Leading dimension of the Vdiff array.
   *  @param fluxes  Output of the diffusive fluxes.
   *             Flat vector with the m_nsp in the inner loop.
   *             length = ldx * ndim
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
  }

  // Get the Mixture diffusion coefficients  [m^2/s]
  /*
   *  The mixture diffusion coefficients are not well defined 
   *  in the context of LiquidTransport because the Stefan Maxwell 
   *  equation is solved.  Here the mixture diffusion coefficients 
   *  are defined according to Ficks law: 
   *  \f[
   *     X_i \vec{V_i} = -D_i \nabla X_i. 
   *  \f]
   *  Solving Ficks Law for \f$ D_i \f$ gives a mixture diffusion 
   *  coefficient
   *  \f[
   *     D_i = - X_i \vec{V_i} / ( \nabla X_i ). 
   *  \f]
   *  If \f$ \nabla X_i = 0 \f$ this is undefined and the 
   *  nonsensical value -1 is returned.  
   *  
   *  Note that this evaluation of \f$ \vec{V_i} \f$  requires 
   *  a solve of the Stefan Maxwell equation making this  
   *  determination of the mixture averaged diffusion coefficients 
   *  a {\em slow} method for obtaining diffusion coefficients.  
   * 
   *  Also note that the Stefan Maxwell solve will be based upon 
   *  the thermodynamic state (including gradients) most recently
   *  set.  Gradients can be set specifically using set_Grad_V,
   *  set_Grad_X and set_Grad_T or through calls to 
   *  getSpeciesFluxes, getSpeciesFluxesES, getSpeciesVdiff, 
   *  getSpeciesVdiffES, etc.
   *  
   *  @param d vector of mixture diffusion coefficients
   *          units = m2 s-1. length = number of species
   */
  void LiquidTransport::getMixDiffCoeffs(doublereal* const d) {

    update_T();
    update_C();

    update_Grad_lnAC();

    stefan_maxwell_solve();

    for ( int n = 0; n < m_nDim; n++) {
      for (int k = 0; k < m_nsp; k++) {
	if ( m_Grad_X[n*m_nsp + k] != 0.0 ) {
	  d[n*m_nsp + k] = - m_Vdiff(k,n) * m_molefracs[k] 
	    / m_Grad_X[n*m_nsp + k];
	} else {
	  //avoid divide by zero with nonsensical response
	  d[n*m_nsp + k] = - 1.0; 
	}
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
      m_thermo->getMassFractions(DATA_PTR(m_massfracs));
      m_thermo->getMoleFractions(DATA_PTR(m_molefracs));
      m_thermo->getConcentrations(DATA_PTR(m_concentrations));
      concTot_ = 0.0;
      concTot_tran_ = 0.0;
      for (int k = 0; k < m_nsp; k++) {
	m_molefracs[k] = fmaxx(0.0, m_molefracs[k]);
	m_molefracs_tran[k] = fmaxx(MIN_X, m_molefracs[k]);
	m_massfracs_tran[k] = fmaxx(MIN_X, m_massfracs[k]);
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

  /*************************************************************************
   *
   *    methods to update species temperature-dependent properties 
   *
   *************************************************************************/

  /**
   * Update the temperature-dependent parts of the species
   * thermal conductivity internally using calls to the
   * appropriate LTPspecies subclass. 
   */
  void LiquidTransport::updateCond_T() {

    int k;

    for (k = 0; k < m_nsp; k++) {
      m_lambdaSpecies[k] = m_lambdaTempDep_Ns[k]->getSpeciesTransProp() ;
    }
    m_cond_temp_ok = true;
    m_cond_mix_ok = false;
  }


  //! Update the binary Stefan-Maxwell diffusion coefficients 
  //! wrt T using calls to the appropriate LTPspecies subclass
  void LiquidTransport::updateDiff_T() {

    m_diffMixModel->getMatrixTransProp( m_bdiff );
    m_diff_temp_ok = true;
    m_diff_mix_ok = false;
  }


  //! Update the pure-species viscosities functional dependence on concentration.
  void LiquidTransport::updateViscosities_C() {
    m_visc_conc_ok = true;
  }


  /**
   * Updates the array of pure species viscosities internally 
   * using calls to the appropriate LTPspecies subclass.
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


  //!  Update the temperature-dependent hydrodynamic radius terms
  //!  for each species internally  using calls to the
  //!  appropriate LTPspecies subclass
  void LiquidTransport::updateHydrodynamicRadius_T() {
    int k;

    for (k = 0; k < m_nsp; k++) {
      m_hydrodynamic_radius[k] = m_radiusTempDep_Ns[k]->getSpeciesTransProp() ;
    } 
    m_radi_temp_ok = true;
    m_radi_mix_ok = false;
  }

  //!  Updates the internal value of the gradient of the  
  //!  logarithm of the activity coefficients, which is 
  //!  used in the gradient of the chemical potential. 
  /**
   * Evaluate the gradients of the activity coefficients 
   * as they alter the diffusion coefficient.  
   *
   *  The gradient of the chemical potential can be written in terms of 
   *  gradient of the logarithm of the mole fraction times a correction
   *  associated with the gradient of the activity coefficient relative to 
   *  that of the mole fraction.  Specifically, the gradients of the 
   *  logarithms of each are involved according to the formula 
   
   *  \f[
   *      \nabla \mu_k = RT \nabla ( \ln X_k ) 
   *      \left[ 1 + \nabla ( \ln \gamma_k ) / \nabla ( \ln X_k ) \right]
   *  \f]
   *  
   * The required quantity is the derivitive of the logarithm of the 
   * activity coefficient with respect to the derivative of the 
   * logarithm of the mole fraction (or whatever concentration 
   * variable we are using to express chemical potential.
   *
   * Updates the vector over species i:
   * \[
   *    \partial \left[ \ln ( \gamma_i ) \right] 
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


  /*
   *
   *    Solve for the diffusional velocities in the Stefan-Maxwell equations
   *
   */
  //! Solve the stefan_maxell equations for the diffusive fluxes.
  /*
   * The diffusive mass flux of species \e k is computed 
   * using the Stefan-Maxwell equation
   * \f[
   *     X_i \nabla \mu_i 
   *                     = RT \sum_i \frac{X_i X_j}{D_{ij}} 
   *                           ( \vec{V}_j - \vec{V}_i )
   * \f]
   * to determine the diffusion velocity and 
   * \f[
   *      \vec{N}_i = C_T X_i \vec{V}_i
   * \f]
   * to determine the diffusion flux.  Here \f$ C_T \f$ is the 
   * total concentration of the mixture [kmol/m^3], \f$ D_{ij} \f$
   * are the Stefa-Maxwell interaction parameters in [m^2/s],
   * \f$ \vec{V}_{i} \f$ is the diffusion velocity of species \e i,
   * \f$ \mu_i \f$ is the electrochemical potential of species \e i.
   *
   * The diffusion velocity is relative to an average velocity 
   * that can be computed on a mole-weighted 
   * or mass-weighted basis, or the diffusion velocities may 
   * be specified as relative to a specific species (i.e. a 
   * solvent) all according to the \verbatim <velocityBasis> 
   * \endverbatim input parameter.
   *
   * One of the Stefan Maxwell equations is replaced by the appropriate
   * definition of the mass-averaged velocity, the mole-averaged velocity 
   * or the specification that velocities are relative to that 
   * of one species.
   */
  void LiquidTransport::stefan_maxwell_solve() {
    int i, j, a;
    doublereal tmp;
    m_B.resize(m_nsp, m_nDim, 0.0);
    m_A.resize(m_nsp, m_nsp, 0.0);

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
      //equation for the reference velocity 
      for (j = 0; j < m_nsp; j++) {
	if ( m_velocityBasis == VB_MOLEAVG )
	  m_A(0,j) = m_molefracs_tran[j];
	else if ( m_velocityBasis == VB_MASSAVG )
	  m_A(0,j) = m_massfracs_tran[j];
	else if ( ( m_velocityBasis >= 0 ) 
		  && ( m_velocityBasis < m_nsp ) )
	  // use species number m_velocityBasis as reference velocity
	  if ( m_velocityBasis == j ) m_A(0,j) = 1.0; 
	else 
	  throw CanteraError("LiquidTransport::stefan_maxwell_solve",
			     "Unknown reference velocity provided.");
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
	    m_A(i,i) -=   tmp;
	    m_A(i,j)  = + tmp;
	  }
	}
      }

      //! invert and solve the system  Ax = b. Answer is in m_B
      solve(m_A, m_B);
  	
      break;
    case 2:  /* 2-D approximation */
      m_B(0,0) = 0.0;
      m_B(0,1) = 0.0;
      //equation for the reference velocity 
      for (j = 0; j < m_nsp; j++) {
	if ( m_velocityBasis == VB_MOLEAVG )
	  m_A(0,j) = m_molefracs_tran[j];
	else if ( m_velocityBasis == VB_MASSAVG )
	  m_A(0,j) = m_massfracs_tran[j];
	else if ( ( m_velocityBasis >= 0 ) 
		  && ( m_velocityBasis < m_nsp ) )
	  // use species number m_velocityBasis as reference velocity
	  if ( m_velocityBasis == j ) m_A(0,j) = 1.0; 
	else 
	  throw CanteraError("LiquidTransport::stefan_maxwell_solve",
			     "Unknown reference velocity provided.");
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
	    m_A(i,i) -=   tmp;
	    m_A(i,j)  = + tmp;
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
      //equation for the reference velocity 
      for (j = 0; j < m_nsp; j++) {
	if ( m_velocityBasis == VB_MOLEAVG )
	  m_A(0,j) = m_molefracs_tran[j];
	else if ( m_velocityBasis == VB_MASSAVG )
	  m_A(0,j) = m_massfracs_tran[j];
	else if ( ( m_velocityBasis >= 0 ) 
		  && ( m_velocityBasis < m_nsp ) )
	  // use species number m_velocityBasis as reference velocity
	  if ( m_velocityBasis == j ) m_A(0,j) = 1.0; 
	else 
	  throw CanteraError("LiquidTransport::stefan_maxwell_solve",
			     "Unknown reference velocity provided.");
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
	    m_A(i,i) -=   tmp;
	    m_A(i,j)  = + tmp;
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
