/**
 *
 *  @file SimpleTransport.h
 *   Header file defining class SimpleTransport
 */
#ifndef CT_SIMPLETRAN_H
#define CT_SIMPLETRAN_H



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

    
  //! Class LiquidTransport implements mixture-averaged transport
  //! properties for liquid phases.
  /*!
   *  The model is based on that
   *  described by Newman, Electrochemical Systems
   *
   *  The velocity of species i may be described by the
   *  following equation p. 297 (12.1)
   *
   *   \f[
   *     c_i \nabla \mu_i = R T \sum_j \frac{c_i c_j}{c_T D_{ij}}
   *         (\mathbf{v}_j - \mathbf{v}_i)
   *   \f]
   *
   * This as written is degenerate by 1 dof.
   *
   * To fix this we must add in the definition of the mass averaged
   * velocity of the solution. We will call the simple bold-faced 
   *  \f$\mathbf{v} \f$
   * symbol the mass-averaged velocity. Then, the relation
   * between \f$\mathbf{v}\f$ and the individual species velocities is
   * \f$\mathbf{v}_i\f$
   *
   *  \f[
   *      \rho_i \mathbf{v}_i =  \rho_i \mathbf{v} + \mathbf{j}_i
   *  \f]
   *   where \f$\mathbf{j}_i\f$ are the diffusional fluxes of species i
   *   with respect to the mass averaged velocity and
   *
   *  \f[
   *       \sum_i \mathbf{j}_i = 0
   *  \f]
   * 
   *  and
   *
   *  \f[
   *     \sum_i \rho_i \mathbf{v}_i =  \rho \mathbf{v} 
   *  \f]
   * 
   * Using these definitions, we can write
   *
   *  \f[
   *      \mathbf{v}_i =  \mathbf{v} + \frac{\mathbf{j}_i}{\rho_i}
   *  \f]
   *
   *
   *  \f[
   *     c_i \nabla \mu_i = R T \sum_j \frac{c_i c_j}{c_T D_{ij}}
   *         (\frac{\mathbf{j}_j}{\rho_j} - \frac{\mathbf{j}_i}{\rho_i})
   *         =  R T \sum_j \frac{1}{D_{ij}}
   *         (\frac{x_i \mathbf{j}_j}{M_j} - \frac{x_j \mathbf{j}_i}{M_i})
   *  \f]
   *
   * The equations that we actually solve are
   *
   *  \f[
   *     c_i \nabla \mu_i =
   *         =  R T \sum_j \frac{1}{D_{ij}}
   *         (\frac{x_i \mathbf{j}_j}{M_j} - \frac{x_j \mathbf{j}_i}{M_i})
   *  \f]
   *  and we replace the 0th equation with the following:
   *
   *  \f[
   *       \sum_i \mathbf{j}_i = 0
   *  \f]
   *
   *  When there are charged species, we replace the rhs with the
   *  gradient of the electrochemical potential to obtain the
   *  modified equation
   *  
   *  \f[
   *     c_i \nabla \mu_i + c_i F z_i \nabla \Phi 
   *         =  R T \sum_j \frac{1}{D_{ij}}
   *         (\frac{x_i \mathbf{j}_j}{M_j} - \frac{x_j \mathbf{j}_i}{M_i})
   *  \f]
   *
   * With this formulation we may solve for the diffusion velocities,
   * without having to worry about what the mass averaged velocity
   * is.
   *
   *  <H2> Viscosity Calculation  </H2>
   * 
   *  The viscosity calculation may be broken down into two parts.
   *  In the first part, the viscosity of the pure species are calculated
   *  In the second part, a mixing rule is applied, based on the
   *  Wilkes correlation, to yield the mixture viscosity.
   *  
   *
   *
   */
  class SimpleTransport : public Transport {
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
    SimpleTransport(thermo_t* thermo = 0, int ndim = 1);

    //!Copy Constructor for the %LiquidThermo object.
    /*!
     * @param right  %LiquidTransport to be copied
     */
    SimpleTransport(const SimpleTransport &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %LiquidTransport object to be copied 
     *                 into the current one.
     */
    SimpleTransport& operator=(const  SimpleTransport& right);
    
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
    virtual ~SimpleTransport() {}

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     *
     * @param tr  Transport parameters for all of the species
     *            in the phase.
     */
    virtual bool initLiquid(LiquidTransportParams& tr);

    friend class TransportFactory;


    //! Return the model id for this transport parameterization
    virtual int model() const {
      return cSimpleTransport; 
    }

    //! overloaded base class methods

    //! Returns the mixture viscosity of the solution
    /*!
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
    virtual doublereal viscosity();

    //! Returns the pure species viscosities
    /*!
     *  The pure species viscosities are to be given in an Arrhenius 
     * form in accordance with activated-jump-process dominated transport.
     */
    virtual void getSpeciesViscosities(doublereal* const visc);

    //! Returns the binary diffusion coefficients
    /*!
     *   @param ld 
     *   @param d   
     */
    virtual void getBinaryDiffCoeffs(const int ld, doublereal* const d);

    //! Get the Mixture diffusion coefficients
    /*!
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


    //! Returns the mixture thermal conductivity of the solution
    /*!
     * The thermal is computed using the general mixture rules
     * specified in the variable compositionDepType_.
     *
     *  Controlling update boolean = m_condmix_ok
     *
     *  Units are in W/m/K  or equivalently kg m / s3 / K
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

    virtual doublereal thermalConductivity();

    //! Get the electrical Mobilities (m^2/V/s).
    /*!
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
    virtual void getMobilities(doublereal* const mobil_e);

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
    virtual void getFluidMobilities(doublereal* const mobil_f);

    //! Specify the valpdaue of the gradient of the voltage
    /*!
     *
     * @param grad_V Gradient of the voltage (length num dimensions);
     */
    virtual void set_Grad_V(const doublereal* const grad_V);

    //! Specify the value of the gradient of the temperature
    /*!
     *
     * @param grad_V Gradient of the temperature (length num dimensions);
     */
    virtual void set_Grad_T(const doublereal* const grad_T);

    //! Specify the value of the gradient of the MoleFractions
    /*!
     *
     * @param grad_X Gradient of the mole fractions(length nsp * num dimensions);
     */
    virtual void set_Grad_X(const doublereal* const grad_X);


    /**
     * @param ndim The number of spatial dimensions (1, 2, or 3).
     * @param grad_T The temperature gradient (ignored in this model).
     * @param ldx  Leading dimension of the grad_X array.
     * The diffusive mass flux of species \e k is computed from
     *
     *
     */
     virtual void getSpeciesFluxes(int ndim, 
				  const doublereal* grad_T, 
				  int ldx, const doublereal* grad_X, 
				  int ldf, doublereal* fluxes);

    //!  Return the species diffusive mass fluxes wrt to
    //!  the mass averaged velocity,
    /*!
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
    virtual void getSpeciesFluxesExt(int ldf, doublereal* fluxes);

  protected:

    //! Handles the effects of changes in the Temperature, internally
    //! within the object.
    /*!
     *  This is called whenever a transport property is requested.
     *  The first task is to check whether the temperature has changed
     *  since the last call to update_T().
     *  If it hasn't then an immediate return is carried out.
     *
     *     @internal
     *
     * @return  Returns true if the temperature has changed, and false otherwise
     */
    virtual bool update_T();

    //! Handles the effects of changes in the mixture concentration
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
     */ 
    virtual bool update_C();

    //!  Update the temperature-dependent viscosity terms.
    //!  Updates the array of pure species viscosities, and the 
    //!  weighting functions in the viscosity mixture rule.
    /*!
     * The flag m_visc_temp_ok is set to true.
     */
    void updateViscosity_T();

    //! Update the temperature-dependent parts of the mixture-averaged 
    //! thermal conductivity.     
    void updateCond_T();

    //! Update the concentration parts of the viscosities
    /*!
     *  Internal routine is run whenever the update_boolean
     *  is false. This routine will calculate
     *  internal values for the species viscosities.
     *
     * @internal
     */
    void updateViscosities_C();
 
    //! Update the binary diffusion coefficients wrt T.
    /*!
     *   These are evaluated
     *   from the polynomial fits at unit pressure (1 Pa).
     */
    void updateDiff_T();


  private:

    //! Number of species in the mixture
    int m_nsp;

    //! Temperature dependence type
    /*!
     *     The following coefficients are allowed to have simple
     *     temperature dependencies:
     *         mixture viscosity
     *         mixture thermal conductivity
     *         diffusitivy
     *
     *  Types of temperature dependencies:
     *     0  - Independent of temperature (only one implemented so far)
     *     1  - extended arrhenius form
     *     2  - power law form
     */
    int tempDepType_;

    //! Composition dependence of the transport properties
    /*!
     *   The following coefficients are allowed to have simple
     *   composition dependencies
     *       mixture viscosity
     *       mixture thermal conductivity
     *       
     *
     *   Types of composition dependencies
     *    0 - Solvent values (i.e., species 0) contributes only
     *    1 - linear combination of mole fractions; 
     */
    int compositionDepType_;

    bool useHydroRadius_;

    //! Boolean indicating whether electro-migration term should be
    //! added
    /*!
     *
     */
    bool doMigration_;

    //! Minimum temperature applicable to the transport property eval
    doublereal m_tmin;

    //! Maximum temperature applicable to the transport property evaluator
    doublereal m_tmax;

    //! Local Copy of the molecular weights of the species
    /*!
     *  Length is Equal to the number of species in the mechanism.
     */
    vector_fp  m_mw;

    //! Pure species viscosities in Arrhenius temperature-dependent form.
    std::vector<Coeff_T_>  m_coeffVisc_Ns; 
  
    //! Pure species thermal conductivities in Arrhenius temperature-dependent form.
    /*!
     *
     */
    std::vector<Coeff_T_>  m_coeffLambda_Ns; 
  

    //! Pure species viscosities in Arrhenius temperature-dependent form.
    std::vector<Coeff_T_>  m_coeffDiff_Ns; 

    
    std::vector<Coeff_T_>  m_coeffHydroRadius_Ns; 
  

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


    // property values

  

    //! Vector of Species Diffusivities
    /*!
     *   Depends on the temperature. We have set the pressure dependence
     *   to zero for this liquid phase constituitve model
     *
     *   units m2/s
     */
    vector_fp m_diffSpecies;

    //! Species viscosities
    /*!
     *   Viscosity of the species 
     *   Length = number of species
     *
     *   Depends on the temperature. We have set the pressure dependence
     *   to zero for this model
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
    vector_fp  m_condSpecies;

    //! State of the mole fraction vector.
    int m_iStateMF;

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

    //! Mean molecular weight
    doublereal meanMolecularWeight_;

    //! Density
    doublereal dens_;

    //! Local copy of the charge of each species
    /*!
     *  Contains the charge of each species (length m_nsp)
     */
    vector_fp m_chargeSpecies;
  
    //! Current Temperature -> locally storred
    /*!
     * This is used to test whether new temperature computations
     * should be performed.
     */
    doublereal m_temp;


    //! Current value of the pressure
    doublereal m_press;

 
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
 
    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_diff_mix_ok;

    //! Boolean indicating that binary diffusion coeffs are current
    bool m_diff_temp_ok;

    //! Flag to indicate that the pure species conductivities
    //! are current wrt the temperature
    bool m_cond_temp_ok;

    //! Boolean indicating that mixture conductivity is current
    bool m_cond_mix_ok;


    //! Number of dimensions 
    /*!
     * Either 1, 2, or 3
     */
    int m_nDim;

  private:
    
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






