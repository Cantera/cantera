/**
 *
 *  @file LiquidTransport.h
 *   Header file defining class LiquidTransport
 */
/*
 * $Revision: 1.9 $
 * $Date: 2009/03/27 18:24:39 $
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
#include "LiquidTransportParams.h"

namespace Cantera {

  const int LVISC_CONSTANT     = 0;
  const int LVISC_WILKES       = 1;
  const int LVISC_MIXTUREAVG   = 2;

  const int LDIFF_MIXDIFF_UNCORRECTED     = 0;
  const int LDIFF_MIXDIFF_FLUXCORRECTED  = 1;
  const int LDIFF_MULTICOMP_STEFANMAXWELL  = 2;



  class TransportParams;

    
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
  class LiquidTransport : public Transport {
  public:

    //! default constructor
    LiquidTransport(thermo_t* thermo = 0, int ndim = 1);

    //!Copy Constructor for the %LiquidThermo object.
    /*!
     * @param right  ThermoPhase to be copied
     */
    LiquidTransport(const LiquidTransport &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %ThermoPhase object to be copied into the
     *                 current one.
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
    virtual ~LiquidTransport() {}

    //! Return the model id for this transport parameterization
    virtual int model() {
      return cLiquidTransport; 
    }

    //! overloaded base class methods

    //! Returns the viscosity of the solution
    /*!
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
     *
     * Controlling update boolean m_viscmix_ok
     */ 
    virtual doublereal viscosity();

    //! Returns the pure species viscosities
    /*!
     *
     * 
     */
    virtual void getSpeciesViscosities(doublereal* const visc);

    virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Return the thermal conductivity of the solution
    /*!
     * The thermal conductivity is computed from the following mixture rule:
     *   \f[
     *    \lambda = 0.5 \left( \sum_k X_k \lambda_k 
     *     + \frac{1}{\sum_k X_k/\lambda_k}\right)
     *   \f]
     *
     *  Controlling update boolean = m_condmix_ok
     */
    virtual doublereal thermalConductivity();

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


    //! Get the Mobilities
    /*!
     * @param mobil
     */
    virtual void getMobilities(doublereal* const mobil);

    //! Specify the value of the gradient of the voltage
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

   virtual void update_Grad_lnAC();

  protected:
    //! Handles the effects of changes in the Temperature, internally
    //! within the object.
    /*!
     *  This is called whenever a transport property is
     *  requested.  
     *  The first task is to check whether the temperature has changed
     *  since the last call to update_T().
     *  If it hasn't then an immediate return is carried out.
     *
     *
     *   Note this should be a lightweight function since it's
     *   part of all of the interfaces.
     *
     *     @internal
     */ 
    virtual void update_temp();

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
    virtual void update_conc();

  public:
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

    //!  Return the species diffusive mass fluxes
    /*!
     *
     *
     *
     * @param ndim The number of spatial dimensions (1, 2, or 3).
     * @param grad_T The temperature gradient (ignored in this model).
     * @param ldx  Leading dimension of the grad_X array.
     * The diffusive mass flux of species \e k is computed from
     *
     *
     */
    virtual void getSpeciesDiffusiveMassFluxes(doublereal* const fluxes);
    
    /**
     * @param ndim The number of spatial dimensions (1, 2, or 3).
     * @param grad_T The temperature gradient (ignored in this model).
     * @param ldx  Leading dimension of the grad_X array.
     * The diffusive mass flux of species \e k is computed from
     *
     *
     */
    virtual void getSpeciesFluxesExt(int ldf, doublereal* fluxes);


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



    //! Solve the stefan_maxell equations for the diffusive fluxes.
    void stefan_maxwell_solve();


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

    // polynomial fits
    vector<vector<int> >         m_poly;

    //! Polynomial coefficients of the viscosity
    /*!
     * These express the temperature dependendence of the pures
     * species viscosities.
     */
    std::vector<vector_fp> viscCoeffsVector_;

    //! Polynomial coefficients of the conductivities
    /*!
     * These express the temperature dependendence of the pures
     * species conductivities
     */
    vector<vector_fp>            m_condcoeffs;

    //! Polynomial coefficients of the binary diffusion coefficients
    /*!
     * These express the temperature dependendence of the
     * binary diffusivities. An overall pressure dependence is then
     * added.
     */
    vector<vector_fp>            m_diffcoeffs;
 

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
     *  ck  m_Grad_mu[n*m_nsp + k]
     */
    vector_fp m_ck_Grad_mu;

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

    //! Species viscosities
    /*!
     *  Viscosity of the species
     *   Length = number of species
     *
     *   Depends on the temperature. We have set the pressure dependence
     *   to zero for this liquid phase constituitve model
     *
     * controlling update boolean -> m_visc_temp_ok
     */
    vector_fp viscSpecies_;

    //! Sqrt of the species viscosities
    /*!
     * The sqrt(visc) is used in the mixing formulas
     * Length = m_nsp
     *
     * Depends on the temperature and perhaps pressure, but
     * not the species concentrations
     *
     * controlling update boolean m_visc_temp_ok
     */
    vector_fp m_sqvisc;

    //! Internal value of the species individual thermal conductivities
    /*!
     * Then a mixture rule is applied to get the solution conductivities
     *
     * Depends on the temperature and perhaps pressure, but
     * not the species concentrations
     *
     * controlling update boolean -> m_cond_temp_ok
     */
    vector_fp  m_cond;

    //! Polynomials of the log of the temperature
    vector_fp m_polytempvec;

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

    doublereal meanMolecularWeight_;
    doublereal dens_;

    //! Local copy of the charge of each species
    /*!
     *  Contains the charge of each species (length m_nsp)
     */
    vector_fp m_chargeSpecies;

  
    vector_fp volume_specPM_;

    vector_fp actCoeffMolar_;

    vector_fp lnActCoeffMolarDelta_;

    //! Stefan-Maxwell Diffusion Coefficients at T, P and C
    /*!
     *   These diffusion coefficients are considered to be
     *  a function of Temperature, Pressure, and Concentration.
     */
    DenseMatrix m_DiffCoeff_StefMax;

    //! Viscosity model
    /*!
     *
     */
    int viscosityModel_;

    //! viscosity weighting functions
    DenseMatrix m_phi;         

    //! Matrix of the ratios of the species molecular weights
    /*!
     * m_wratjk(i,j) = (m_mw[j]/m_mw[k])**0.25
     */
    DenseMatrix m_wratjk;

    //! Matrix of the ratios of the species molecular weights
    /*!
     * m_wratkj1(i,j) = (1.0 + m_mw[k]/m_mw[j])**0.5
     */
    DenseMatrix m_wratkj1;

    //! RHS to the stefan-maxwell equation
    DenseMatrix   m_B;

    //! Matrix for the stefan maxwell equation.
    DenseMatrix m_A;

    //! Internal storage for the species LJ well depth
    vector_fp   m_eps;

    //! Internal storage for species polarizability
    vector_fp   m_alpha;

    

    //! Current Temperature -> locally storred
    /*!
     * This is used to test whether new temperature computations
     * should be performed.
     */
    doublereal m_temp;

    //! Current log(T)
    doublereal m_logt;

    //! Current value of kT
    doublereal m_kbt;

    //! Current Temperature **0.5
    doublereal m_sqrt_t;

    //! Current Temperature **0.25
    doublereal m_t14;

    //! Current Temperature **1.5
    doublereal m_t32;

    //! Current temperature function
    /*!
     *  This is equal to sqrt(Boltzmann * T)
     */
    doublereal m_sqrt_kbt;

    //! Current value of the pressure
    doublereal m_press;

    //! Solution of the flux system
    /*!
     *  This is the mass flux of species k
     *  in units of kg m-3 s-1.
     */
    Array2D m_flux;

    //! Saved value of the mixture thermal conductivity
    doublereal m_lambda;

    //! Saved value of the mixture viscosity
    doublereal m_viscmix;

    // work space
    vector_fp  m_spwork;

    //! Internal Function
  protected:
    //!  Update the temperature-dependent viscosity terms.
    //!  Updates the array of pure species viscosities, and the 
    //!  weighting functions in the viscosity mixture rule.
    /*!
     * The flag m_visc_ok is set to true.
     */
    void updateViscosity_temp();

    //! Update the temperature-dependent parts of the mixture-averaged 
    //! thermal conductivity.     
    void updateCond_temp();

    //! Update the concentration parts of the viscosities
    /*!
     *  Internal routine is run whenever the update_boolean
     *  m_visc_conc_ok is false. This routine will calculate
     *  internal values for the species viscosities.
     *
     * @internal
     */
    void updateViscosities_conc();
 
    //! Update the binary diffusion coefficients wrt T.
    /*!
     *   These are evaluated
     *   from the polynomial fits at unit pressure (1 Pa).
     */
    void updateDiff_temp();

  private:    
    //! Boolean indicating that the top-level mixture viscosity is current
    /*!
     *  This is turned false for every change in T, P, or C.
     */
    bool m_visc_mix_ok;

    //! Boolean indicating that weight factors wrt viscosity is current
    bool m_visc_temp_ok;
 
    //! Flag to indicate that the pure species viscosities
    //! are current wrt the temperature
    bool m_visc_conc_ok;

    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_diff_mix_ok;

    //! Boolean indicating that binary diffusion coeffs are current
    bool m_diff_temp_ok;

    //! Flag to indicate that the pure species conductivities
    //! are current wrt the temperature
    bool m_cond_temp_ok;

    //! Boolean indicating that mixture conductivity is current
    bool m_cond_mix_ok;

    //! Mode for fitting the species viscosities
    /*!
     * Either its CK_Mode or its cantera mode
     * in CK_Mode visc is fitted to a polynomial
     * in Cantera mode sqrt(visc) is fitted.
     */
    int m_mode;

    //! Internal storage for the diameter - diameter
    //! species interactions
    DenseMatrix m_diam;

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
  };
}
#endif






