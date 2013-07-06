/**
 *  @file PecosTransport.h
 *   Header file defining class PecosTransport
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_PECOSTRAN_H
#define CT_PECOSTRAN_H

#include "TransportBase.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{

class GasTransportParams;

/**
 * Class PecosTransport implements mixture-averaged transport
 * properties for ideal gas mixtures.
 */
class PecosTransport : public Transport
{

public:
    virtual int model() const {
        return cPecosTransport;
    }

    //! Viscosity of the mixture
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
     */
     virtual doublereal viscosity();

    virtual void getSpeciesViscosities(doublereal* const visc) {
        update_T();
        updateViscosity_T();
        copy(m_visc.begin(), m_visc.end(), visc);
    }

    //! Return the thermal diffusion coefficients
    /*!
     * For this approximation, these are all zero.
     */
    virtual void getThermalDiffCoeffs(doublereal* const dt);

    //! Returns the mixture thermal conductivity
    /*!
     * This is computed using the lumped model,
     * \f[
     *    k = k^{tr} + k^{ve}
     * \f]
     * where,
     * \f[
     *    k^{tr}= 5/2 \mu_s C_{v,s}^{trans} + \mu_s C_{v,s}^{rot}
     * \f]
     * and,
     * \f[
     *    k^{ve}= \mu_s C_{v,s}^{vib} + \mu_s C_{v,s}^{elec}
     * \f]
     *
     * The thermal conductivity is computed using the Wilke mixture rule.
     * \f[
     *     k = \sum_s \frac{k_s X_s}{\sum_j \Phi_{s,j} X_j}.
     * \f]
     * Here \f$ k_s \f$ is the conductivity of pure species \e s,
     * and
     * \f[
     * \Phi_{s,j} = \frac{\left[1
     * + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_s}}\right)}\right]^2}
     * {\sqrt{8}\sqrt{1 + M_s/M_j}}
     * \f]
     * @see updateCond_T();
     * @todo Reconcile these these formulas with the implementation
     */
    virtual doublereal thermalConductivity();

    //! binary diffusion coefficients
    /*!
     *  Using Ramshaw's self-consistent Effective Binary Diffusion
     *  (1990, J. Non-Equilib. Thermo)
     */
    virtual void getBinaryDiffCoeffs(const size_t ld, doublereal* const d);

    //! Mixture-averaged diffusion coefficients [m^2/s].
    /*!
    *   For the single species case or the pure fluid case the routine returns
    *   the self-diffusion coefficient. This is need to avoid a NaN result.
    */
    virtual void getMixDiffCoeffs(doublereal* const d);

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    //! These are the coefficients for calculating the molar diffusive fluxes
    //! from the species mole fraction gradients, computed according to
    //! Eq. 12.176 in "Chemically Reacting Flow":
    //!
    //! \f[  D_{km}^* = \frac{1-X_k}{\sum_{j \ne k}^K X_j/\mathcal{D}_{kj}} \f]
    //!
    //! @param[out] d vector of mixture-averaged diffusion coefficients for
    //!     each species, length m_nsp.
    void getMixDiffCoeffsMole(doublereal* const d);

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    //! These are the coefficients for calculating the diffusive mass fluxes
    //! from the species mass fraction gradients, computed according to
    //! Eq. 12.178 in "Chemically Reacting Flow":
    //!
    //! \f[  \frac{1}{D_{km}} = \sum_{j \ne k}^K \frac{X_j}{\mathcal{D}_{kj}} +
    //!     \frac{X_k}{1-Y_k} \sum_{j \ne k}^K \frac{Y_j}{\mathcal{D}_{kj}} \f]
    //!
    //! @param[out] d vector of mixture-averaged diffusion coefficients for
    //!     each species, length m_nsp.
    void getMixDiffCoeffsMass(doublereal* const d);

    virtual void getMobilities(doublereal* const mobil);
    virtual void update_T();

    /**
     *  This is called the first time any transport property is requested from
     *  Mixture after the concentrations have changed.
     */
    virtual void update_C();

    //! Get the species diffusive mass fluxes wrt to the mass averaged
    //! velocity, given the gradients in mole fraction and temperature
    /*!
     * The diffusive mass flux of species \e k is computed from
     * \f[
     * \vec{j}_k = -n M_k D_k \nabla X_k + \frac{\rho_k}{\rho} \sum_r n M_r D_r \nabla X_r
     * \f]
     * This neglects pressure, forced and thermal diffusion.
     * Units for the returned fluxes are kg m-2 s-1.
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
     * @param fluxes  Output of the diffusive mass fluxes
     *             Flat vector with the m_nsp in the inner loop.
     *             length = ldx * ndim
     */
    virtual void getSpeciesFluxes(size_t ndim,
                                  const doublereal* const grad_T,
                                  size_t ldx,
                                  const doublereal* const grad_X,
                                  size_t ldf, doublereal* const fluxes);

    //! Initialize the transport object
    /*!
     * Here we change all of the internal dimensions to be sufficient.
     * We get the object ready to do property evaluations.
     *
     * @param tr  Transport parameters for all of the species in the phase.
     */
    virtual bool initGas(GasTransportParams& tr);

    /**
     * Reads the transport table specified (currently defaults to internal file)
     *
     * Reads the user-specified transport table, appending new species
     * data and/or replacing default species data.
     */
    void read_blottner_transport_table();

    friend class TransportFactory;

protected:
    PecosTransport();

private:

    //! Calculate the pressure from the ideal gas law
    doublereal pressure_ig() const {
        return (m_thermo->molarDensity() * GasConstant *
                m_thermo->temperature());
    }

    // mixture attributes
    int m_nsp;
    vector_fp  m_mw;

    // polynomial fits
    std::vector<vector_fp> m_visccoeffs;
    std::vector<vector_fp> m_condcoeffs;
    std::vector<vector_fp> m_diffcoeffs;
    vector_fp  m_polytempvec;

    // blottner fits
    //int species = 20;
    double a[500], b[500], c[500];

    // property values
    DenseMatrix                  m_bdiff;
    vector_fp                    m_visc;
    vector_fp                    m_sqvisc;
    vector_fp                    m_cond;

    vector_fp                    m_molefracs;

    std::vector<std::vector<int> > m_poly;
    std::vector<vector_fp> m_astar_poly;
    std::vector<vector_fp> m_bstar_poly;
    std::vector<vector_fp> m_cstar_poly;
    std::vector<vector_fp> m_om22_poly;
    DenseMatrix          m_astar;
    DenseMatrix          m_bstar;
    DenseMatrix          m_cstar;
    DenseMatrix          m_om22;

    DenseMatrix m_phi;            // viscosity weighting functions
    DenseMatrix m_wratjk, m_wratkj1;

    vector_fp   m_zrot;
    vector_fp   m_crot;
    vector_fp   m_cinternal;
    vector_fp   m_eps;
    vector_fp   m_alpha;
    vector_fp   m_dipoleDiag;

    doublereal m_temp, m_logt, m_kbt, m_t14, m_t32;
    doublereal m_sqrt_kbt, m_sqrt_t;

    vector_fp  m_sqrt_eps_k;
    DenseMatrix m_log_eps_k;
    vector_fp  m_frot_298;
    vector_fp  m_rotrelax;

    doublereal m_lambda;
    doublereal m_viscmix;

    // work space
    vector_fp  m_spwork;

    void updateThermal_T();

    /**
     * Update the temperature-dependent viscosity terms. Updates the array of
     * pure species viscosities, and the weighting functions in the viscosity
     * mixture rule. The flag m_visc_ok is set to true.
     */
    void updateViscosity_T();

    /**
     * Update the temperature-dependent parts of the mixture-averaged
     * thermal conductivity.
     *
     * Calculated as,
     * \f[
     *    k= \mu_s (5/2 * C_{v,s}^{trans} + C_{v,s}^{rot} + C_{v,s}^{vib}
     * \f]
     */
    void updateCond_T();

    /**
     * Update the pure-species viscosities. (Pa-s) = (kg/m/sec)
     *
     * Using Blottner fit for viscosity. Defines kinematic viscosity
     * of the form
     * \f[
     *   \mu_s\left(T\right) = 0.10 \exp\left(A_s\left(\log T\right)^2 + B_s\log T + C_s\right)
     * \f]
     * where \f$ A_s \f$, \f$ B_s \f$, and \f$ C_s \f$ are constants.
     */
    void updateSpeciesViscosities();

    /**
     * Update the binary diffusion coefficients. These are evaluated
     * from the polynomial fits at unit pressure (1 Pa).
     */
    void updateDiff_T();
    void correctBinDiffCoeffs();
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
    bool m_debug;

    // specific heats
    vector_fp            cv_rot;
    vector_fp            cp_R;
    vector_fp            cv_int;

};
}
#endif
