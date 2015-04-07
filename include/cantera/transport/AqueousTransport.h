/**
 *  @file AqueousTransport.h
 *   Header file defining class AqueousTransport
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_AQUEOUSTRAN_H
#define CT_AQUEOUSTRAN_H

// Cantera includes
#include "TransportBase.h"
#include "cantera/numerics/DenseMatrix.h"
#include "TransportParams.h"
#include "LiquidTransportParams.h"

namespace Cantera
{

class LiquidTransportParams;

//! Class AqueousTransport implements mixture-averaged transport
//! properties for brine phases.
/*!
 *  The model is based on that described by Newman, Electrochemical Systems
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
 * With this formulation we may solve for the diffusion velocities, without
 * having to worry about what the mass averaged velocity is.
 *
 *  <H2> Viscosity Calculation  </H2>
 *
 *  The viscosity calculation may be broken down into two parts.
 *  In the first part, the viscosity of the pure species are calculated
 *  In the second part, a mixing rule is applied, based on the
 *  Wilkes correlation, to yield the mixture viscosity.
 * @ingroup tranprops
 */
class AqueousTransport : public Transport
{
public:
    AqueousTransport();

    virtual int model() const {
        return cAqueousTransport;
    }

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
     * Controlling update boolean = m_viscwt_ok
     *
     *  @param visc     Vector of species viscosities
     */
    virtual void getSpeciesViscosities(doublereal* const visc);

    //! Return a vector of Thermal diffusion coefficients [kg/m/sec].
    /*!
     * The thermal diffusion coefficient \f$ D^T_k \f$ is defined
     * so that the diffusive mass flux of species <I>k</I> induced by the
     * local temperature gradient is given by the following formula
     *
     *    \f[
     *         M_k J_k = -D^T_k \nabla \ln T.
     *    \f]
     *
     *   The thermal diffusion coefficient can be either positive or negative.
     *
     *  In this method we set it to zero.
     *
     * @param dt On return, dt will contain the species thermal
     *           diffusion coefficients.  Dimension dt at least as large as
     *           the number of species. Units are kg/m/s.
     */
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

    virtual void getBinaryDiffCoeffs(const size_t ld, doublereal* const d);

    //! Get the Mixture diffusion coefficients
    /*!
     * For the single species case or the pure fluid case the routine returns
     * the self-diffusion coefficient. This is need to avoid a NaN result.
     *  @param d vector of mixture diffusion coefficients
     *          units = m2 s-1. length = number of species
     */
    virtual void getMixDiffCoeffs(doublereal* const d);

    virtual void getMobilities(doublereal* const mobil_e);

    virtual void getFluidMobilities(doublereal* const mobil_f);

    //! Specify the value of the gradient of the voltage
    /*!
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
     * @param grad_X Gradient of the mole fractions(length nsp * num dimensions);
     */
    virtual void set_Grad_X(const doublereal* const grad_X);

    //! Handles the effects of changes in the Temperature, internally
    //! within the object.
    /*!
     *  This is called whenever a transport property is requested. The first
     *  task is to check whether the temperature has changed since the last
     *  call to update_T(). If it hasn't then an immediate return is carried
     *  out.
     */
    virtual void update_T();

    //! Handles the effects of changes in the mixture concentration
    /*!
     *  This is called the first time any transport property is requested from
     *  Mixture after the concentrations have changed.
     */
    virtual void update_C();

    virtual void getSpeciesFluxes(size_t ndim,  const doublereal* const grad_T,
                                  size_t ldx, const doublereal* const grad_X,
                                  size_t ldf, doublereal* const fluxes);

    //!  Return the species diffusive mass fluxes wrt to the specified averaged velocity,
    /*!
     *   This method acts similarly to getSpeciesFluxesES() but
     *   requires all gradients to be preset using methods set_Grad_X(), set_Grad_V(), set_Grad_T().
     *   See the documentation of getSpeciesFluxesES() for details.
     *
     *  units = kg/m2/s
     *
     * Internally, gradients in the in mole fraction, temperature
     * and electrostatic potential contribute to the diffusive flux
     *
     * The diffusive mass flux of species \e k is computed from the following formula
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
     *  @param ldf     Stride of the fluxes array. Must be equal to or greater than the number of species.
     *  @param fluxes  Output of the diffusive fluxes. Flat vector with the m_nsp in the inner loop.
     *                   length = ldx * ndim
     */
    virtual void getSpeciesFluxesExt(size_t ldf, doublereal* const fluxes);

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     *
     * @param tr  Transport parameters for all of the species in the phase.
     */
    virtual bool initLiquid(LiquidTransportParams& tr);

    friend class TransportFactory;

    /**
     * Return a structure containing all of the pertinent parameters
     * about a species that was used to construct the Transport
     * properties in this object.
     *
     * @param k Species number to obtain the properties about.
     */
    class LiquidTransportData getLiquidTransportData(int k);

    //! Solve the Stefan-Maxwell equations for the diffusive fluxes.
    void stefan_maxwell_solve();

private:
    //! Local Copy of the molecular weights of the species
    /*!
     *  Length is Equal to the number of species in the mechanism.
     */
    vector_fp  m_mw;

    //! Polynomial coefficients of the viscosity
    /*!
     * These express the temperature dependence of the pure species viscosities.
     */
    std::vector<vector_fp>            m_visccoeffs;

    //! Polynomial coefficients of the conductivities
    /*!
     * These express the temperature dependence of the pure species conductivities
     */
    std::vector<vector_fp>            m_condcoeffs;

    //! Polynomial coefficients of the binary diffusion coefficients
    /*!
     * These express the temperature dependence of the binary diffusivities.
     * An overall pressure dependence is then added.
     */
    std::vector<vector_fp>            m_diffcoeffs;

    //! Internal value of the gradient of the mole fraction vector
    /*!
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
     *  Generally, if a transport property needs this in its evaluation it
     *  will look to this place to get it.
     *
     *  No internal property is precalculated based on gradients. Gradients
     *  are assumed to be freshly updated before every property call.
     */
    vector_fp m_Grad_T;

    //! Internal value of the gradient of the Electric Voltage
    /*!
     *  Generally, if a transport property needs this in its evaluation it
     *  will look to this place to get it.
     *
     *  No internal property is precalculated based on gradients. Gradients
     *  are assumed to be freshly updated before every property call.
     */
    vector_fp m_Grad_V;

    //! Gradient of the electrochemical potential
    /*!
     *  m_nsp is the number of species in the fluid
     *  k is the species index
     *  n is the dimensional index (x, y, or z)
     *
     *    m_Grad_mu[n*m_nsp + k]
     */
    vector_fp m_Grad_mu;

    // property values

    //! Array of Binary Diffusivities
    /*!
     * This has a size equal to nsp x nsp
     * It is a symmetric matrix.
     * D_ii is undefined.
     *
     * units m2/sec
     */
    DenseMatrix                  m_bdiff;

    //! Species viscosities
    /*!
     *  Viscosity of the species
     *   Length = number of species
     *
     * Depends on the temperature and perhaps pressure, but
     * not the species concentrations
     *
     * controlling update boolean -> m_spvisc_ok
     */
    vector_fp m_visc;

    //! Sqrt of the species viscosities
    /*!
     * The sqrt(visc) is used in the mixing formulas
     * Length = m_nsp
     *
     * Depends on the temperature and perhaps pressure, but
     * not the species concentrations
     *
     * controlling update boolean m_spvisc_ok
     */
    vector_fp m_sqvisc;

    //! Internal value of the species individual thermal conductivities
    /*!
     * Then a mixture rule is applied to get the solution conductivities
     *
     * Depends on the temperature and perhaps pressure, but
     * not the species concentrations
     *
     * controlling update boolean -> m_spcond_ok
     */
    vector_fp  m_cond;

    //! Polynomials of the log of the temperature
    vector_fp m_polytempvec;

    //! State of the mole fraction vector.
    int m_iStateMF;

    //! Local copy of the mole fractions of the species in the phase
    /*!
     * Update info?
     * length = m_nsp
     */
    vector_fp m_molefracs;

    //! Local copy of the concentrations of the species in the phase
    /*!
     * Update info?
     * length = m_nsp
     */
    vector_fp m_concentrations;

    //! Local copy of the charge of each species
    /*!
     *  Contains the charge of each species (length m_nsp)
     */
    vector_fp m_chargeSpecies;

    //! Stefan-Maxwell Diffusion Coefficients at T, P and C
    /*!
     *   These diffusion coefficients are considered to be
     *  a function of Temperature, Pressure, and Concentration.
     */
    DenseMatrix m_DiffCoeff_StefMax;

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
    Array2D   m_B;

    //! Matrix for the stefan maxwell equation.
    DenseMatrix m_A;

    //! Internal storage for the species LJ well depth
    vector_fp   m_eps;

    //! Internal storage for species polarizability
    vector_fp   m_alpha;

    //! Current Temperature -> locally stored
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
    Array2D   m_flux;

    //! saved value of the mixture thermal conductivity
    doublereal m_lambda;

    //! Saved value of the mixture viscosity
    doublereal m_viscmix;

    //! work space of size m_nsp
    vector_fp  m_spwork;

    //! Update the temperature-dependent viscosity terms.
    /*!
     *  Updates the array of pure species viscosities, and the weighting
     *  functions in the viscosity mixture rule. The flag m_visc_ok is set to
     *  true.
     */
    void updateViscosity_T();

    //! Update the temperature-dependent parts of the mixture-averaged
    //! thermal conductivity.
    void updateCond_T();

    //! Update the species viscosities
    /*!
     *  Internal routine is run whenever the update_boolean
     *  m_spvisc_ok is false. This routine will calculate
     *  internal values for the species viscosities.
     */
    void updateSpeciesViscosities();

    //! Update the binary diffusion coefficients wrt T.
    /*!
     *   These are evaluated from the polynomial fits at unit pressure (1 Pa).
     */
    void updateDiff_T();

    //! Boolean indicating that mixture viscosity is current
    bool m_viscmix_ok;

    //! Boolean indicating that weight factors wrt viscosity is current
    bool m_viscwt_ok;

    //! Flag to indicate that the pure species viscosities
    //! are current wrt the temperature
    bool m_spvisc_ok;

    //! Boolean indicating that mixture diffusion coeffs are current
    bool m_diffmix_ok;

    //! Boolean indicating that binary diffusion coeffs are current
    bool m_bindiff_ok;

    //! Flag to indicate that the pure species conductivities
    //! are current wrt the temperature
    bool m_spcond_ok;

    //! Boolean indicating that mixture conductivity is current
    bool m_condmix_ok;

    //! Mode for fitting the species viscosities
    /*!
     * Either it's CK_Mode or it's cantera mode
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
    size_t m_nDim;
};
}
#endif
