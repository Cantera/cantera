/**
 *  @file LiquidTransport.h
 *   Header file defining class LiquidTransport
 */
/*
 * $Revision$
 * $Date$
 */

#ifndef CT_LIQUIDTRAN_H
#define CT_LIQUIDTRAN_H



// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

using namespace std;

// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"
#include "TransportParams.h"
#include "LiquidTransportParams.h"

namespace Cantera {

  class LiquidTransportParams;

    
  //! Class LiquidTransport implements models for transport
  //! properties for liquid phases.  
  /*!
   *  Liquid Transport is set up with some flexibility in 
   *  this class.  Transport properties like viscostiy
   *  and thermal conductivity are allowed flexibility within
   *  the constraints of the LiquidTransportProperty and 
   *  LiquidTransportInteractions classes. For species 
   *  diffusion, the LiquidTransport class focuses on 
   *  the Stefan-Maxwell equation to determine the diffusion 
   *  velocities.  Other options for liquid diffusion include 
   *  solvent-dominated diffusion, and a class SolventTransport 
   *  should be forthcoming.  
   *
   *  The class LiquidTransport has several roles.  
   *  -# It brings together the individual species transport 
   *     properties, expressed as subclasses of LTPspecies 
   *     (Liquid Transport Properties of Species), with 
   *     models for the composition dependence of liquid 
   *     transport properties expressed as subclasses of 
   *     LiquidTranInteraction.
   * 
   *  -# It calculates the bulk velocity \f$ \vec{v} \f$ and  
   *     individual species diffusion velocities, \f$ \vec{V_i} \f$ 
   *     using the Stefan-Maxwell equations.  It is 
   *     possible to set a flag to calculate relative to a 
   *     mass-averaged bulk velocity, relative to a mole-averaged 
   *     bulk velocity or relative to a single species velocity 
   *     using the <velocityBasis basis="mass"> keyword.  
   *     Mass-averaged velocities are the default for which the 
   *     diffusion velocities satisfy 
   *     \f[
   *        \sum_{i} Y_i \vec{V_i} = 0
   *     \f] 
   *     for mass fraction \f$ Y_i \f$.  For mole-averaged velocities
   *     \f[
   *        \sum_{i} X_i \vec{V_i} = 0
   *     \f] 
   *     for mole fraction \f$ X_i \f$.
   * 
   *  -# It provides acccess to a number of derived quantities
   *     related to transport properties as described in the 
   *     various methods below.
   *     
   *  
   *  Within LiquidTransport, the state is presumed to be 
   *  defined in terms of the  species mole fraction, 
   *  temperature and pressure.  Charged species are expected 
   *  and quantities like the electric current are computed
   *  based on a combined electrochemcial potential.
   *  
   *
   *  @ingroup tranprops
   */
  class LiquidTransport : public Transport {
  public:

    typedef  vector_fp Coeff_T_;
 

    //! Default constructor.  
    /*!
     * This requires call to initLiquid(LiquidTransportParams& tr)
     * after filling LiquidTransportParams to complete instantiation.
     * The filling of LiquidTransportParams is currently carried out 
     * in the TransportFactory class, but might be moved at some point.
     * 
     * @param thermo  ThermoPhase object holding species information.
     * @param ndim    Number of spatial dimensions.
     */
    LiquidTransport(thermo_t* thermo = 0, int ndim = 1);

    //! Copy Constructor for the %LiquidThermo object.
    /*!
     * @param right  %LiquidTransport to be copied
     */
    LiquidTransport(const LiquidTransport &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %LiquidTransport object to be copied 
     *                 into the current one.
     */
    LiquidTransport&  operator=(const  LiquidTransport& right);
    
    //! Duplication routine for objects which inherit from
    //! %Transport
    /*!
     *  This virtual routine can be used to duplicate %Transport objects
     *  inherited from %Transport even if the application only has
     *  a pointer to %Transport to work with.
     *
     *  These routines are basically wrappers around the derived copy
     *  constructor.
     */
    virtual Transport *duplMyselfAsTransport() const;


    //! virtual destructor
    virtual ~LiquidTransport();

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     * A lot of the input required to do property evaluations is 
     * contained in the LiquidTransportParams class that is 
     * filled in TransportFactory. 
     *
     * @param tr  Transport parameters for all of the species
     *            in the phase.
     */
    virtual bool initLiquid(LiquidTransportParams& tr);

    friend class TransportFactory;


    //! Return the model id for this transport parameterization
    virtual int model() const {
      return cLiquidTransport; 
    }


    //! Returns the viscosity of the solution
    /*!
     *  The viscosity calculation is handled by subclasses of 
     *  LiquidTranInteraction as specified in the input file.  
     *  These in turn employ subclasses of LTPspecies to 
     *  determine the individual species viscosities.
     */ 
    virtual doublereal viscosity();

    //! Returns the pure species viscosities for all species
    /*!
     *  The pure species viscosities are evaluated using the 
     *  appropriate subclasses of LTPspecies as specified in the 
     *  input file.
     *
     * @param visc  array of length "number of species"
     *              to hold returned viscosities.
     */
    virtual void getSpeciesViscosities(doublereal* const visc);

    //! Returns the hydrodynamic radius for all species 
    /*!
     *  The species hydrodynamic radii are evaluated using the 
     *  appropriate subclasses of LTPspecies as specified in the 
     *  input file.
     *
     * @param radius  array of length "number of species"
     *                to hold returned radii.
     */
    virtual void getSpeciesHydrodynamicRadius(doublereal* const radius);

    //! Returns the binary diffusion coefficients
    /*!
     *   The binary diffusion coefficients are specified in the input
     *   file through the LiquidTransportInteractions class.  These
     *   are the binary interaction coefficients employed in the 
     *   Stefan-Maxwell equation.
     *   
     *   @param ld  number of species in system
     *   @param d   vector of binary diffusion coefficients
     *          units = m2 s-1. length = ld*ld = (number of species)^2
     */
    virtual void getBinaryDiffCoeffs(const int ld, doublereal* const d);

    //! Get the Mixture diffusion coefficients
    /*!
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
     *  a \e slow method for obtaining diffusion coefficients.  
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
    virtual void getMixDiffCoeffs(doublereal* const d);


    //! Return the thermal diffusion coefficients
    /*!
     *  These are all zero for this simple implementaion
     *
     *  @param dt thermal diffusion coefficients
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Return the thermal conductivity of the solution
    /*!
     *  The thermal conductivity calculation is handled by subclasses of 
     *  LiquidTranInteraction as specified in the input file.  
     *  These in turn employ subclasses of LTPspecies to 
     *  determine the individual species thermal condictivities.
     */ 
    virtual doublereal thermalConductivity();

    //! Get the Electrical mobilities (m^2/V/s).
    /*!
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
     *  a \e slow method for obtaining diffusion coefficients.  
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
    virtual void getMobilities(doublereal* const mobil_e);

    //! Get the fluid mobilities (s kmol/kg).
    /*!
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
     *  a \e slow method for obtaining diffusion coefficients.  
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
    virtual void getFluidMobilities(doublereal* const mobil_f);

    //! Specify the value of the gradient of the voltage
    /*!
     *
     * @param grad_V Gradient of the voltage (length num dimensions);
     */
    virtual void set_Grad_V(const doublereal* const grad_V);

    //! Specify the value of the gradient of the temperature
    /*!
     * @param grad_T Gradient of the temperature (length num dimensions);
     */
    virtual void set_Grad_T(const doublereal* const grad_T);

    //! Specify the value of the gradient of the MoleFractions
    /*!
     *
     * @param grad_X Gradient of the mole fractions(length nsp * num dimensions);
     */
    virtual void set_Grad_X(const doublereal* const grad_X);

     //! Compute the mixture electrical conductivity from 
     //! the Stefan-Maxwell equation.
     /**
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
      * 
      */
     doublereal getElectricConduct( );

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
     void getElectricCurrent(int ndim, 
			     const doublereal* grad_T, 
			     int ldx, 
			     const doublereal* grad_X, 
			     int ldf, 
			     const doublereal* grad_V, 
			     doublereal* current) ;
     


    //! Get the species diffusive velocities wrt to 
    //! the averaged velocity, 
    //! given the gradients in mole fraction and temperature
    /*!
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
     * @param ldx  Leading dimension of the grad_X array 
     *              (usually equal to m_nsp but not always)
     * @param grad_X Gradients of the mole fraction
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     * @param ldf  Leading dimension of the fluxes array 
     *              (usually equal to m_nsp but not always)
     * @param Vdiff  Output of the diffusive velocities.
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     */
    virtual void getSpeciesVdiff(int ndim, 
				 const doublereal* grad_T, 
				 int ldx, 
				 const doublereal* grad_X,
				 int ldf, 
				 doublereal* Vdiff);

    //! Get the species diffusive mass fluxes wrt to 
    //! the averaged velocity, 
    //! given the gradients in mole fraction, temperature 
    //! and electrostatic potential.
    /*!
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
     * @param ldx  Leading dimension of the grad_X array 
     *              (usually equal to m_nsp but not always)
     * @param grad_X Gradients of the mole fraction
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     * @param ldf  Leading dimension of the fluxes array 
     *              (usually equal to m_nsp but not always)
     * @param grad_Phi Gradients of the electrostatic potential
     *                 (length = ndim)
     * @param fluxes  Output of the diffusive mass fluxes
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     */
    virtual void getSpeciesVdiffES(int ndim, 
				    const doublereal* grad_T, 
				    int ldx, 
				    const doublereal* grad_X,
				    int ldf, 
				    const doublereal* grad_Phi,
				    doublereal* Vdiff) ;



     //!  Return the species diffusive mass fluxes wrt to
     //!  the averaged velocity in [kmol/m^2/s].
     /**
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
     virtual void getSpeciesFluxes(int ndim, 
				  const doublereal* grad_T, 
				  int ldx, const doublereal* grad_X, 
				  int ldf, doublereal* fluxes);

     //!  Return the species diffusive mass fluxes wrt to
     //!  the averaged velocity in [kmol/m^2/s].
     /**
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
     virtual void getSpeciesFluxesES(int ndim, 
				     const doublereal* grad_T, 
				     int ldx, 
				     const doublereal* grad_X, 
				     int ldf, 
				     const doublereal* grad_Phi,
				     doublereal* fluxes);

     //!  Return the species diffusive velocities relative to 
     //!  the averaged velocity.  
     /**
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
    virtual void getSpeciesVdiffExt(int ldf, doublereal* Vdiff);

     //!  Return the species diffusive fluxes relative to 
     //!  the averaged velocity.  
     /**
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
    virtual void getSpeciesFluxesExt(int ldf, doublereal* fluxes);

  protected:
    //! Returns true if temperature has changed, 
    //! in which case flags are set to recompute transport properties.
    /*!
     *  This is called whenever a transport property is requested.
     *  The first task is to check whether the temperature has changed
     *  since the last call to update_T().
     *  If it hasn't then an immediate return is carried out.
     *
     *
     *   Note this should be a lightweight function since it's
     *   part of all of the interfaces.
     *
     *     @internal
     *
     * @return  Returns true if the temperature has changed, and false otherwise
     */
    virtual bool update_T();

    //! Returns true if mixture composition has changed, 
    //! in which case flags are set to recompute transport properties.
    /*!
     *   This is called for every interface call to check whether
     *   the concentrations have changed. Concentrations change
     *   whenever the pressure or the mole fraction has changed.
     *   If it has changed, the recalculations should be done.
     *
     *   Note this should be a lightweight function since it's
     *   part of all of the interfaces.
     *
     *   @internal
     *
     * @return  Returns true if the mixture composition has changed, and false otherwise.
     */ 
    virtual bool update_C();

    
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
     virtual void update_Grad_lnAC();


    //! Solve the stefan_maxell equations for the diffusive fluxes.
    /**
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
    void stefan_maxwell_solve();

    //!  Update the temperature-dependent viscosity terms.
    //!  Updates the array of pure species viscosities, and the 
    //!  weighting functions in the viscosity mixture rule.
    /*!
     * The flag m_visc_ok is set to true.
     */
    void updateViscosity_T();

    //!  Update the temperature-dependent hydrodynamic radius terms
    //!  for each species 
    /*!
     * The flag m_radi_temp_ok is set to true.
     */
    void updateHydrodynamicRadius_T();

    //! Update the temperature-dependent parts of the mixture-averaged 
    //! thermal conductivity.     
    void updateCond_T();

    //! Update the concentration parts of the viscosities
    /*!
     *  Internal routine is run whenever the update_boolean
     *  m_visc_conc_ok is false. Currently there is no concentration 
     *  dependence for the pure species viscosities.
     *
     * @internal
     */
    void updateViscosities_C();
 
    //! Update the concentration dependence of the hydrodynamic radius
    /*!
     *  Internal routine is run whenever the update_boolean
     *  m_radi_conc_ok is false. Currently there is no concentration 
     *  dependence for the hydrodynamic radius.
     *
     * @internal
     */
    void updateHydrodynamicRadius_C();
 
    //! Update the binary diffusion coefficients wrt T.
    /*!
     *   These are evaluated
     *   from the polynomial fits at unit pressure (1 Pa).
     */
    void updateDiff_T();


  private:

    //! Number of species in the mixture
    int m_nsp;

    //! Minimum temperature applicable to the transport property eval
    doublereal m_tmin;

    //! Maximum temperature applicable to the transport property evaluator
    doublereal m_tmax;

    //! Local Copy of the molecular weights of the species
    /*!
     *  Length is Equal to the number of species in the mechanism.
     */
    vector_fp  m_mw;

    //! Viscosity temperature dependence type
    /*!
     *  Types of temperature dependencies:
     *     0  - Independent of temperature (only one implemented so far)
     *     1  - extended arrhenius form
     *     2  - polynomial in temperature form
     */
    std::vector<LTPspecies*> m_viscTempDep_Ns;

    //! Viscosity mixing model type
    /*!
     *  Types of mixing models supported:
     *     2  - Mole fraction weighting of species viscosities
     *     3  - Mass fraction weighting of species viscosities
     *     4  - Mole fraction weighting of logarithms of species viscosities
     */
    LiquidTranInteraction *m_viscMixModel;

    //! Thermal conductivity temperature dependence type
    /*!
     *  Types of temperature dependencies:
     *     0  - Independent of temperature (only one implemented so far)
     *     1  - extended arrhenius form
     *     2  - polynomial in temperature form
     */
    std::vector<LTPspecies*> m_lambdaTempDep_Ns;

    //! Thermal conductivity mixing model type
    /*!
     *  Types of mixing models supported:
     *     2  - Mole fraction weighting of species viscosities
     *     3  - Mass fraction weighting of species viscosities
     */
    LiquidTranInteraction *m_lambdaMixModel;
 
   //! Diffusion coefficient temperature dependence type
    /*!
     *  Types of temperature dependencies:
     *     0  - Independent of temperature (only one implemented so far)
     *     1  - extended arrhenius form
     *     2  - polynomial in temperature form
     */
    std::vector<LTPspecies*> m_diffTempDep_Ns;

    //! Species diffusivity mixing model type
    /*!
     *  Types of mixing models supported:
     *     5  - Pairwise interactions -- Setfan-Maxwell diffusion coefficients
     */
    LiquidTranInteraction *m_diffMixModel;

    //! Setfan-Maxwell diffusion coefficients
    DenseMatrix m_diff_Dij;


    std::vector<bool> useHydroRadius_;

   //!Hydrodynamic radius temperature dependence type
    /*!
     *  Types of temperature dependencies:
     *     0  - Independent of temperature
     *     1  - extended arrhenius form
     *     2  - polynomial in temperature form
     */
    std::vector<LTPspecies*> m_radiusTempDep_Ns;

    //! Species hydrodynamic radius
    vector_fp  m_hydrodynamic_radius;

    //! Hydrodynamic radius mixing model type
    /*!
     *  Types of mixing models supported:
     *     0  - No mixing model allowed
     */
    LiquidTranInteraction *m_radiusMixModel;


    //! Polynomial coefficients of the binary diffusion coefficients
    /*!
     * These express the temperature dependendence of the
     * binary diffusivities. An overall pressure dependence is then
     * added.
     */
    /*
    std::vector<vector_fp>            m_diffcoeffs;
    */


    //! Internal value of the gradient of the mole fraction vector
    /*!
     *  Note, this is the only gradient value that can and perhaps
     *  should reflect the true state of the mole fractions in the
     *  application solution vector. In other words no cropping or
     *  massaging of the values to make sure they are above zero
     *  should occur. - developing ....
     *
     *  m_nsp is the number of species in the fluid
     *  k is the species index
     *  n is the dimensional index (x, y, or z). It has a length
     *    equal to m_nDim
     *  
     *    m_Grad_X[n*m_nsp + k]
     */
    vector_fp m_Grad_X;

    //! Gradient of the logarithm of the activity coefficients 
    //! with respect to the logarithm of the mole fraction, plus one.
    /*!
     *  This quantity appears in the gradient of the chemical potential.  
     *  It multiplies the gradient of the mole fraction, and in this way 
     * serves to "modify" the diffusion coefficient.  
     *
     *    m_Grad_X[k] = 1 + \partial \left[ \ln ( \gamma_i ) \right] 
     *                  / \partial \left[ \ln ( \X_i  ) \right] 
     * 
     * Note that where "molefraction is used here, whatever 
     * concentration-related variable applies, so that if 
     * molality is the concentration variable, the gradient of the 
     * activity coefficient should be with respect to the molality. 
     *  
     */
    vector_fp m_Grad_lnAC;

    //! Internal value of the gradient of the Temperature vector
    /*!
     *  Generally, if a transport property needs this 
     *  in its evaluation it will look to this place
     *  to get it. 
     *
     *  No internal property is precalculated based on gradients.
     *  Gradients are assumed to be freshly updated before
     *  every property call.
     */
    vector_fp m_Grad_T;

    //! Internal value of the gradient of the Pressure vector
    /*!
     *  Generally, if a transport property needs this 
     *  in its evaluation it will look to this place
     *  to get it. 
     *
     *  No internal property is precalculated based on gradients.
     *  Gradients are assumed to be freshly updated before
     *  every property call.
     */
    vector_fp m_Grad_P;

    //! Internal value of the gradient of the Electric Voltage
    /*!
     *  Generally, if a transport property needs this 
     *  in its evaluation it will look to this place
     *  to get it. 
     *
     *  No internal property is precalculated based on gradients.
     *  Gradients are assumed to be freshly updated before
     *  every property call.
     */
    vector_fp m_Grad_V;

    //! Gradient of the electrochemical potential
    /*!
     *  m_nsp is the number of species in the fluid
     *  k is the species index
     *  n is the dimensional index (x, y, or z)
     *  
     *  \f[
     *     m_Grad_mu[n*m_nsp + k]
     *  \f]
     */
    vector_fp m_Grad_mu;

    // property values

    //! Array of Binary Diffusivities
    /*!
     *   Depends on the temperature. We have set the pressure dependence
     *   to zero for this liquid phase constituitve model
     *
     *  This has a size equal to nsp x nsp
     *  It is a symmetric matrix.
     *  D_ii is the self diffusion coefficient. D_ii is not
     *  needed except for when there is one species in the mixture.
     *
     * units m2/sec
     */
    DenseMatrix  m_bdiff;

    //! Species viscosities and their logarithm
    /*!
     *  Viscosity of the species and its logarithm
     *   Length = number of species
     *
     *   Depends on the temperature. We have set the pressure dependence
     *   to zero for this liquid phase constituitve model
     *
     * controlling update boolean -> m_visc_temp_ok
     */
    vector_fp m_viscSpecies;

    //! Internal value of the species individual thermal conductivities
    /*!
     * Then a mixture rule is applied to get the solution conductivities
     *
     * Depends on the temperature and perhaps pressure, but
     * not the species concentrations
     *
     * controlling update boolean -> m_cond_temp_ok
     */
    vector_fp  m_lambdaSpecies;

    //! State of the mole fraction vector.
    int m_iStateMF;

    //! Local copy of the mass fractions of the species in the phase
    /*!
     *  The mass fraction vector comes from the ThermoPhase object. 
     *
     * length = m_nsp
     */
    vector_fp m_massfracs;

    //! Local copy of the mass fractions of the species in the phase
    /** 
     * This version of the mass fraction vector is adjusted to a
     * minimum lower bound of MIN_X for use in transport calculations.
     */ 
    vector_fp m_massfracs_tran;

    //! Local copy of the mole fractions of the species in the phase
    /*!
     *  The mole fractions here are assumed to be bounded by 0.0 and 1.0
     *  and they are assumed to add up to one exactly. This mole
     *  fraction vector comes from the ThermoPhase object. Derivative
     *  quantities from this are referred to as bounded.
     *
     * Update info?
     * length = m_nsp
     */
    vector_fp m_molefracs;

    //! Non-zero mole fraction vector used in transport property calculations
    /*!
     *  The mole fractions here are assumed to be bounded by MIN_X and 1.0
     *  and they may not be assumed to add up to one. This
     *  mole fraction vector is created from the ThermoPhase object.
     *  Derivative quantities of this use the _tran suffix.
     *
     * Update info?
     * length = m_nsp
     */
    vector_fp m_molefracs_tran;

    vector_fp Xdelta_;

    //! Local copy of the concentrations of the species in the phase
    /*!
     *  The concentrations are consistent with the m_molefracs
     *  vector which is bounded and sums to one.
     *
     * Update info?
     * length = m_nsp
     */
    vector_fp m_concentrations;

    //! Local copy of the total concentration.
    /*!
     * This is consistent with the m_concentrations[] and
     *  m_molefracs[] vector.
     */
    doublereal concTot_;

    //! Local copy of the total concentration.
    /*!
     *  This is consistent with the x_molefracs_tran vector and
     *  with the concTot_ number;
     */
    doublereal concTot_tran_;

    //! Mean molecular mass
    doublereal meanMolecularWeight_;

    //! Density
    doublereal dens_;

    //! Local copy of the charge of each species
    /*!
     *  Contains the charge of each species (length m_nsp)
     */
    vector_fp m_chargeSpecies;

    //! Specific volume for each species.  Local copy from thermo object.
    vector_fp m_volume_spec;

    vector_fp m_actCoeff;

    //! RHS to the stefan-maxwell equation
    DenseMatrix   m_B;

    //! Matrix for the stefan maxwell equation.
    DenseMatrix m_A;

    //! Current Temperature -> locally storred
    /*!
     * This is used to test whether new temperature computations
     * should be performed.
     */
    doublereal m_temp;

    //! Current value of the pressure
    doublereal m_press;

    //! Solution of the flux system
    /*!
     *  This is the mass flux of species k
     *  in units of kg m-3 s-1.
     */
    Array2D m_flux;

    //! Solution of the Stefan Maxwell equation
    /*!
     *  This is the diffusion velocity of species k
     *  in units of m/s and relative to the mole-averaged velocity.
     */
    Array2D m_Vdiff;

    //! Saved value of the mixture thermal conductivity
    doublereal m_lambda;

    //! Saved value of the mixture viscosity
    doublereal m_viscmix;

    //! work space
    /*!
     *   Length is equal to m_nsp
     */
    vector_fp  m_spwork;



  private:    
    //! Boolean indicating that the top-level mixture viscosity is current
    /*!
     *  This is turned false for every change in T, P, or C.
     */
    bool m_visc_mix_ok;

    //! Boolean indicating that weight factors wrt viscosity is current
    bool m_visc_temp_ok;
 
    //! Flag to indicate that the pure species viscosities
    //! are current wrt the concentration
    bool m_visc_conc_ok;

    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_radi_mix_ok;

    //! Boolean indicating that temperature dependence of
    //! hydrodynamic radius is current 
    bool m_radi_temp_ok;
 
    //! Flag to indicate that the hydrodynamic radius is current
    //! is current wrt the concentration
    bool m_radi_conc_ok;

    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_diff_mix_ok;

    //! Boolean indicating that binary diffusion coeffs are current
    bool m_diff_temp_ok;

    //! Flag to indicate that the pure species conductivities
    //! are current wrt the temperature
    bool m_cond_temp_ok;

    //! Boolean indicating that mixture conductivity is current
    bool m_cond_mix_ok;

    //! Mode indicator for transport models -- currently unused.
    int m_mode;

    //! Debugging flags
    /*!
     *  Turn on to get debugging information
     */
    bool m_debug;

    //! Number of dimensions 
    /*!
     * Either 1, 2, or 3
     */
    int m_nDim;

    //! Throw an exception if this method is invoked. 
    /*!
     * This probably indicates something is not yet implemented.
     *
     * @pram msg   Indicates the member function which is not implemented
     */
    doublereal err(std::string msg) const;

  };
}
#endif






