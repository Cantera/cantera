/**
 *  @file MixTransport.h
 *    Headers for the MixTransport object, which models transport properties
 *    in ideal gas solutions using a mixture averaged approximation
 *    (see \ref tranprops and \link Cantera::MixTransport MixTransport \endlink) .
 *
 */
/* $Author$
 * $Revision$
 * $Date$
 */
// Copyright 2001  California Institute of Technology


#ifndef CT_MIXTRAN_H
#define CT_MIXTRAN_H

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

// STL includes
#include <vector>
#include <string>
#include <map>
#include <numeric>
#include <algorithm>

// Cantera includes
#include "TransportBase.h"
#include "DenseMatrix.h"

namespace Cantera {


  class GasTransportParams;

  
  //! Class MixTransport implements mixture-averaged transport properties for ideal gas mixtures.
  /*!
   *    The model is based on that described by Kee, Coltrin, and Glarborg, "Theoretical and
   *    Practical Aspects of Chemically Reacting Flow Modeling."
   *
   *    
   * The viscosity is computed using the Wilke mixture rule (kg /m /s)
   *
   *    \f[
   *        \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
   *    \f]
   *
   *     Here \f$ \mu_k \f$ is the viscosity of pure species \e k, and
   * 
   *    \f[ 
   *        \Phi_{k,j} = \frac{\left[1 
   *                     + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
   *                     {\sqrt{8}\sqrt{1 + M_k/M_j}}
   *    \f]
   *
   *
   * The thermal conductivity is computed from the following mixture rule:
   *   \f[
   *          \lambda = 0.5 \left( \sum_k X_k \lambda_k  + \frac{1}{\sum_k X_k/\lambda_k} \right)
   *   \f]
   *
   *  It's used to compute the flux of energy due to a thermal gradient
   *
   *   \f[
   *          j_T =  - \lambda  \nabla T
   *   \f]
   *   
   *  The flux of energy has units of energy (kg m2 /s2) per second per area.
   *
   *  The units of lambda are W / m K which is equivalent to kg m / s^3 K.
   *
   *
   */
  class MixTransport : public Transport {

  protected:

    //! Default constructor.  
    /*!
     *
     */
    MixTransport();

  public:

    //!Copy Constructor for the %MixTransport object.
    /*!
     * @param right  %LiquidTransport to be copied
     */
    MixTransport(const MixTransport &right);

    //! Assignment operator
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %LiquidTransport object to be copied 
     *                 into the current one.
     */
    MixTransport& operator=(const  MixTransport& right);
    
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


    //! Destructor
    virtual ~MixTransport();

    //! Return the model id for transport
    /*!
     * @return cMixtureAverage
     */
    virtual int model() const {
      return cMixtureAveraged;
    }

    //! Viscosity of the mixture  (kg /m /s)
    /*!
     * The viscosity is computed using the Wilke mixture rule (kg /m /s)
     *
     *    \f[
     *        \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
     *    \f]
     *
     *     Here \f$ \mu_k \f$ is the viscosity of pure species \e k, and
     * 
     *    \f[ 
     *        \Phi_{k,j} = \frac{\left[1 
     *                     + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
     *                     {\sqrt{8}\sqrt{1 + M_k/M_j}}
     *    \f]
     * 
     *  @return   Returns the viscosity of the mixture  ( units =  Pa s  = kg /m /s)
     *
     * @see updateViscosity_T();
     */
    virtual doublereal viscosity();

    virtual void getSpeciesViscosities(doublereal* visc)
    { update_T();  updateViscosity_T(); copy(m_visc.begin(), m_visc.end(), visc); }

    //! Return the thermal diffusion coefficients
    /*!
     * For this approximation, these are all zero.
     *
     *  Eqns. (12.168) shows how they are used in an expression for the species flux.
     *
     * @param dt  Vector of thermal diffusion coefficients. Units = kg/m/s
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Returns the mixture thermal conductivity (W/m /K)
    /*!
     * The thermal conductivity is computed from the following mixture rule:
     *   \f[
     *          \lambda = 0.5 \left( \sum_k X_k \lambda_k  + \frac{1}{\sum_k X_k/\lambda_k} \right)
     *   \f]
     *
     *  It's used to compute the flux of energy due to a thermal gradient
     *
     *   \f[
     *          j_T =  - \lambda  \nabla T
     *   \f]
     *   
     *  The flux of energy has units of energy (kg m2 /s2) per second per area.
     *
     *  The units of lambda are W / m K which is equivalent to kg m / s^3 K.
     *
     * @return Returns the mixture thermal conductivity, with units of W/m/K
     */
    virtual doublereal thermalConductivity();

    virtual void getBinaryDiffCoeffs(const int ld, doublereal* const d);

    
    //! Returns the Mixture-averaged diffusion coefficients [m^2/s]. 
    /*!
     * Returns the mixture averaged diffusion coefficients for a gas, appropriate for calculating the
     * mass averged diffusive flux with respect to the mass averaged velocity using gradients of the
     * mole fraction. 
     * Note, for the single species case or the pure fluid case the routine returns the self-diffusion coefficient.
     * This is need to avoid a Nan result in the formula below.
     *
     *  This is Eqn. 12.180 from "Chemicaly Reacting Flow"
     *
     *   \f[
     *       D_{km}' = \frac{\left( \bar{M} - X_k M_k \right)}{ \bar{\qquad M \qquad } }  {\left( \sum_{j \ne k} \frac{X_j}{D_{kj}} \right) }^{-1} 
     *   \f]
     *
     *
     *
     *  @param d  Output Vector of mixture diffusion coefficients, \f$  D_{km}'  \f$  ,  for each species (m^2/s).
     *            length m_nsp
     */
    virtual void getMixDiffCoeffs(doublereal* const d);


    virtual void getMobilities(doublereal* const mobil);
    virtual void update_T();
    virtual void update_C();

    //! Get the species diffusive mass fluxes wrt to the mass averaged velocity, 
    //! given the gradients in mole fraction and temperature
    /*!
     *  Units for the returned fluxes are kg m-2 s-1.
     *
     *
     * The diffusive mass flux of species \e k is computed from
     * \f[
     *          \vec{j}_k = -n M_k D_k \nabla X_k.
     * \f]
     *
     *  @param ndim      Number of dimensions in the flux expressions
     *  @param grad_T    Gradient of the temperature
     *                    (length = ndim)
     * @param ldx        Leading dimension of the grad_X array 
     *                   (usually equal to m_nsp but not always)
     * @param grad_X     Gradients of the mole fraction
     *                   Flat vector with the m_nsp in the inner loop.
     *                   length = ldx * ndim
     * @param ldf  Leading dimension of the fluxes array 
     *              (usually equal to m_nsp but not always)
     * @param fluxes  Output of the diffusive mass fluxes
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     */
    virtual void getSpeciesFluxes(int ndim, 
				  const doublereal* grad_T,
				  int ldx,
				  const doublereal* grad_X, 
				  int ldf, doublereal* fluxes);

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     *
     * @param tr  Transport parameters for all of the species
     *            in the phase.
     */
    virtual bool initGas( GasTransportParams& tr );

    friend class TransportFactory;

    
    //! Return a structure containing all of the pertinent parameters about a species that was
    //! used to construct the Transport properties in this object.
    /*!
     * @param kspec Species number to obtain the properties from.
     *
     * @return GasTransportData  returned structure.
     */
    struct GasTransportData getGasTransportData(int kspec) const;

 

  private:

    //! Calculate the pressure from the ideal gas law
    doublereal pressure_ig() const {
      return (m_thermo->molarDensity() * GasConstant *
	      m_thermo->temperature());
    }
    void updateThermal_T();
    void updateViscosity_T();
    void updateCond_T();
    void updateSpeciesViscosities();
    void updateDiff_T();
    void correctBinDiffCoeffs();



    // --------- Member Data -------------

    // mixture attributes
    int m_nsp;

    //! Minimum value of the temperature that this transport parameterization is valid
    doublereal m_tmin;

    //! Maximum value of the temperature that this transport parameterization is valid
    doublereal m_tmax;

    //! Local copy of the species molecular weights.
    vector_fp  m_mw;

    // polynomial fits
    std::vector<vector_fp>            m_visccoeffs;
    std::vector<vector_fp>            m_condcoeffs;
    std::vector<vector_fp>            m_diffcoeffs;
    vector_fp                    m_polytempvec;

    // property values
    DenseMatrix m_bdiff;

    //! vector of species viscosities (kg /m /s)
    /*!
     *  These are used in wilke's rule to calculate the viscosity of the solution
     *  length = m_kk
     */
    vector_fp m_visc;

    //! vector of square root of species viscosities sqrt(kg /m /s)
    /*!
     *  These are used in wilke's rule to calculate the viscosity of the solution
     *  length = m_kk
     */
    vector_fp m_sqvisc;

    //! vector of species thermal conductivities (W/m /K)
    /*!
     *  These are used in wilke's rule to calculate the viscosity of the solution
     *  units = W /m /K = kg m /s^3 /K.
     *  length = m_kk
     */
    vector_fp m_cond;

    //! Vector of species molefractions
    /*!
     *  These are processed so that all mole fractions are >= MIN_X
     *  Length = m_kk
     */
    vector_fp m_molefracs;

    std::vector<std::vector<int> > m_poly;
 
    //! Viscosity Weighting Functions
    DenseMatrix m_phi;  
    DenseMatrix m_wratjk;
    DenseMatrix m_wratkj1;

    //! Rotational relaxation number for the species in the current phase
    /*!
     *  Not used in this routine -> just a passthrough
     *
     * length is the number of species in the phase
     * units are dimensionless
     */
    vector_fp m_zrot;

    //! Dimensionless rotational heat capacity of the species in the current phase
    /*!
     *  These values are 0, 1 and 1.5 for single-molecule, linear, and nonlinear species respectively
     *  length is the number of species in the pahse
     *  units are dimensionless  (Cr / R)
     */
    vector_fp m_crot;
    vector_fp m_cinternal;

    //! Lennard-Jones well-depth of the species in the current phase
    /*!
     * length is the number of species in the phase
     * Units are Joules (Note this is not Joules/kmol) (note, no kmol -> this is a per molecule amount)
     */
    vector_fp m_eps;
    vector_fp m_alpha;
    vector_fp m_dipoleDiag;

    //! Current value of the temperature at which the properties in this object are calculated (Kelvin)
    doublereal m_temp;

    //! Current value of the log of the temperature
    doublereal m_logt;

    //! Current value of Boltzman's constant times the temperature (Joules)
    doublereal m_kbt;

    //! Current value of temperature to 1/4 power
    doublereal m_t14;

    //! Current value of temperature to the 3/2 power
    doublereal m_t32;

    doublereal m_sqrt_kbt;
    doublereal m_sqrt_t;

    vector_fp  m_sqrt_eps_k;
    DenseMatrix m_log_eps_k;
    vector_fp  m_frot_298;
    vector_fp  m_rotrelax;

    doublereal m_lambda;
    doublereal m_viscmix;

    // work space
    vector_fp  m_spwork;

  
    bool m_viscmix_ok;
    bool m_viscwt_ok;
    bool m_spvisc_ok;
    bool m_diffmix_ok;
    bool m_bindiff_ok;
    bool m_abc_ok;
    bool m_spcond_ok;
    bool m_condmix_ok;

    int m_mode;

    DenseMatrix m_epsilon;
    DenseMatrix m_diam;
    DenseMatrix incl;

    //! Debug flag - turns on more printing
    bool m_debug;
  };
}
#endif
