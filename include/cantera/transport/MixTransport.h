/**
 *  @file MixTransport.h
 *    Headers for the MixTransport object, which models transport properties
 *    in ideal gas solutions using a mixture averaged approximation
 *    (see \ref tranprops and \link Cantera::MixTransport MixTransport \endlink) .
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_MIXTRAN_H
#define CT_MIXTRAN_H

#include "GasTransport.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{

class GasTransportParams;

//! Class MixTransport implements mixture-averaged transport properties for ideal gas mixtures.
/*!
 * The model is based on that described by Kee, Coltrin, and Glarborg,
 * "Theoretical and Practical Aspects of Chemically Reacting Flow Modeling."
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
 *  @ingroup tranprops
 */
class MixTransport : public GasTransport
{
protected:

    //! Default constructor.
    MixTransport();

public:

    MixTransport(const MixTransport& right);
    MixTransport& operator=(const  MixTransport& right);
    virtual Transport* duplMyselfAsTransport() const;

    //! Return the model id for transport
    /*!
     * @return cMixtureAverage
     */
    virtual int model() const {
        return cMixtureAveraged;
    }

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

    //! Get the Electrical mobilities (m^2/V/s).
    /*!
     *   This function returns the mobilities. In some formulations
     *   this is equal to the normal mobility multiplied by Faraday's constant.
     *
     *   Here, the mobility is calculated from the diffusion coefficient using the Einstein relation
     *
     *     \f[
     *          \mu^e_k = \frac{F D_k}{R T}
     *     \f]
     *
     * @param mobil  Returns the mobilities of the species in array \c mobil. The array must be
     *               dimensioned at least as large as the number of species.
     */
    virtual void getMobilities(doublereal* const mobil);

    //! Update the internal parameters whenever the temperature has changed
    /*!
     *  This is called whenever a transport property is requested if
     *  the temperature has changed since the last call to update_T().
     */
    virtual void update_T();

    //! Update the internal parameters whenever the concentrations have changed
    /*!
     *  This is called whenever a transport property is requested if the
     *  concentrations have changed since the last call to update_C().
     */
    virtual void update_C();

    //! Get the species diffusive mass fluxes wrt to the mass averaged velocity,
    //! given the gradients in mole fraction and temperature
    /*!
     *  Units for the returned fluxes are kg m-2 s-1.
     *
     * The diffusive mass flux of species \e k is computed from
     * \f[
     *          \vec{j}_k = -n M_k D_k \nabla X_k.
     * \f]
     *
     *  @param ndim      Number of dimensions in the flux expressions
     *  @param grad_T    Gradient of the temperature
     *                    (length = ndim)
     *  @param ldx        Leading dimension of the grad_X array
     *                   (usually equal to m_nsp but not always)
     *  @param grad_X     Gradients of the mole fraction
     *                    Flat vector with the m_nsp in the inner loop.
     *                    length = ldx * ndim
     *  @param ldf        Leading dimension of the fluxes array
     *                   (usually equal to m_nsp but not always)
     *  @param fluxes     Output of the diffusive mass fluxes
     *                    Flat vector with the m_nsp in the inner loop.
     *                       length = ldx * ndim
     */
    virtual void getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                  size_t ldx, const doublereal* const grad_X,
                                  size_t ldf, doublereal* const fluxes);

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     *
     * @param tr  Transport parameters for all of the species
     *            in the phase.
     */
    virtual bool initGas(GasTransportParams& tr);

    friend class TransportFactory;

private:

    //! Calculate the pressure from the ideal gas law
    doublereal pressure_ig() const {
        return (m_thermo->molarDensity() * GasConstant *
                m_thermo->temperature());
    }

    //! Update the temperature dependent parts of the species thermal conductivities
    /*!
     * These are evaluated from the polynomial fits of the temperature and are
     * assumed to be independent of pressure
     */
    void updateCond_T();

private:
    //! Polynomial fits to the thermal conductivity of each species
    /*!
     *  m_condcoeffs[k] is vector of polynomial coefficients for species k
     *  that fits the thermal conductivity
     */
    std::vector<vector_fp>            m_condcoeffs;

    //! vector of species thermal conductivities (W/m /K)
    /*!
     *  These are used in wilke's rule to calculate the viscosity of the
     *  solution. units = W /m /K = kg m /s^3 /K. length = m_kk.
     */
    vector_fp m_cond;

    //! Internal storage for the calculated mixture thermal conductivity
    /*!
     *  Units = W /m /K
     */
    doublereal m_lambda;

    //! Update boolean for the species thermal conductivities
    bool m_spcond_ok;

    //! Update boolean for the mixture rule for the mixture thermal conductivity
    bool m_condmix_ok;
public:
    vector_fp m_eps;
    vector_fp m_sigma;
    vector_fp m_alpha;
    DenseMatrix m_dipole;
    vector_fp m_zrot;
    vector_fp m_crot;
private:
    //! Debug flag - turns on more printing
    bool m_debug;
};
}
#endif
