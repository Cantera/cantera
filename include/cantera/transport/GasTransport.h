/**
 * @file GasTransport.h
 */

#ifndef CT_GAS_TRANSPORT_H
#define CT_GAS_TRANSPORT_H

#include "TransportBase.h"

namespace Cantera
{

//! Class GasTransport implements some functions and properties that are
//! shared by the MixTransport and MultiTransport classes.
//! @ingroup tranprops
class GasTransport : public Transport
{
public:
    GasTransport(const GasTransport& right);
    GasTransport& operator=(const GasTransport& right);

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

    //! Get the pure-species viscosities
    virtual void getSpeciesViscosities(doublereal* const visc) {
        update_T();
        updateViscosity_T();
        std::copy(m_visc.begin(), m_visc.end(), visc);
    }

    //! Returns the matrix of binary diffusion coefficients.
    /*!
     *        d[ld*j + i] = rp * m_bdiff(i,j);
     *
     * @param ld   offset of rows in the storage
     * @param d    output vector of diffusion coefficients. Units of m**2 / s
     */
    virtual void getBinaryDiffCoeffs(const size_t ld, doublereal* const d);

    //! Returns the Mixture-averaged diffusion coefficients [m^2/s].
    /*!
     * Returns the mixture averaged diffusion coefficients for a gas,
     * appropriate for calculating the mass averaged diffusive flux with respect
     * to the mass averaged velocity using gradients of the mole fraction.
     * Note, for the single species case or the pure fluid case the routine
     * returns the self-diffusion coefficient. This is needed to avoid a Nan
     * result in the formula below.
     *
     *  This is Eqn. 12.180 from "Chemically Reacting Flow"
     *
     *   \f[
     *       D_{km}' = \frac{\left( \bar{M} - X_k M_k \right)}{ \bar{\qquad M \qquad } }  {\left( \sum_{j \ne k} \frac{X_j}{D_{kj}} \right) }^{-1}
     *   \f]
     *
     *  @param[out] d  Vector of mixture diffusion coefficients, \f$ D_{km}' \f$ ,
     *      for each species (m^2/s). length m_nsp
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
    virtual void getMixDiffCoeffsMole(doublereal* const d);

    //! Returns the mixture-averaged diffusion coefficients [m^2/s].
    /*! 
     * These are the coefficients for calculating the diffusive mass fluxes
     * from the species mass fraction gradients, computed according to
     * Eq. 12.178 in "Chemically Reacting Flow":
     *
     * \f[
     *     \frac{1}{D_{km}} = \sum_{j \ne k}^K \frac{X_j}{\mathcal{D}_{kj}} +
     *     \frac{X_k}{1-Y_k} \sum_{j \ne k}^K \frac{Y_j}{\mathcal{D}_{kj}}
     * \f]
     *
     * @param[out] d vector of mixture-averaged diffusion coefficients for
     *     each species, length m_nsp.
     */
    virtual void getMixDiffCoeffsMass(doublereal* const d);

protected:
    GasTransport(ThermoPhase* thermo=0);

    virtual bool initGas(GasTransportParams& tr);
    virtual void update_T();
    virtual void update_C() = 0;

    //! Update the temperature-dependent viscosity terms.
    /**
     * Updates the array of pure species viscosities, and the weighting
     * functions in the viscosity mixture rule. The flag m_visc_ok is set to true.
     *
     * The formula for the weighting function is from Poling and Prausnitz,
     * Eq. (9-5.14):
     *  \f[
     *      \phi_{ij} = \frac{ \left[ 1 + \left( \mu_i / \mu_j \right)^{1/2} \left( M_j / M_i \right)^{1/4} \right]^2 }
     *                    {\left[ 8 \left( 1 + M_i / M_j \right) \right]^{1/2}}
     *  \f]
     */
    virtual void updateViscosity_T();

    //! Update the pure-species viscosities. These are evaluated from the
    //! polynomial fits of the temperature and are assumed to be independent
    //! of pressure.
    virtual void updateSpeciesViscosities();

    //! Update the binary diffusion coefficients
    /*!
     * These are evaluated from the polynomial fits of the temperature at the unit pressure of 1 Pa.
     */
    virtual void updateDiff_T();

    //! Vector of species mole fractions. These are processed so that all mole
    //! fractions are >= *Tiny*. Length = m_kk.
    vector_fp m_molefracs;

    //! Internal storage for the viscosity of the mixture  (kg /m /s)
    doublereal m_viscmix;

    //! Update boolean for mixture rule for the mixture viscosity
    bool m_visc_ok;

    //! Update boolean for the weighting factors for the mixture viscosity
    bool m_viscwt_ok;

    //! Update boolean for the species viscosities
    bool m_spvisc_ok;

    //! Update boolean for the binary diffusivities at unit pressure
    bool m_bindiff_ok;

    //! Type of the polynomial fits to temperature. CK_Mode means Chemkin mode.
    //! Currently CA_Mode is used which are different types of fits to temperature.
    int m_mode;

    //! m_phi is a Viscosity Weighting Function. size = m_nsp * n_nsp
    DenseMatrix m_phi;

    //! work space length = m_kk
    vector_fp m_spwork;

    //! vector of species viscosities (kg /m /s). These are used in Wilke's
    //! rule to calculate the viscosity of the solution. length = m_kk.
    vector_fp m_visc;

    //! Polynomial fits to the viscosity of each species. m_visccoeffs[k] is
    //! the vector of polynomial coefficients for species k that fits the
    //! viscosity as a function of temperature.
    std::vector<vector_fp> m_visccoeffs;

    //! Local copy of the species molecular weights.
    vector_fp m_mw;

    //! Holds square roots of molecular weight ratios
    /*!
     *  @code
     *  m_wratjk(j,k)  = sqrt(mw[j]/mw[k])        j < k
     *  m_wratjk(k,j)  = sqrt(sqrt(mw[j]/mw[k]))  j < k
     *  @endcode
     */
    DenseMatrix m_wratjk;

    //! Holds square roots of molecular weight ratios
    /*!
     *  `m_wratjk1(j,k)  = sqrt(1.0 + mw[k]/mw[j])        j < k`
     */
    DenseMatrix m_wratkj1;

    //! vector of square root of species viscosities sqrt(kg /m /s). These are
    //! used in Wilke's rule to calculate the viscosity of the solution.
    //! length = m_kk.
    vector_fp m_sqvisc;

    //! Powers of the ln temperature, up to fourth order
    vector_fp m_polytempvec;

    //! Current value of the temperature at which the properties in this object
    //! are calculated (Kelvin).
    doublereal m_temp;

    //! Current value of Boltzman's constant times the temperature (Joules)
    doublereal m_kbt;

    //! current value of  Boltzman's constant times the temperature.
    //! (Joules) to 1/2 power
    doublereal m_sqrt_kbt;

    //! current value of temperature to 1/2 power
    doublereal m_sqrt_t;

    //! Current value of the log of the temperature
    doublereal m_logt;

    //! Current value of temperature to 1/4 power
    doublereal m_t14;

    //! Current value of temperature to the 3/2 power
    doublereal m_t32;

    //! Polynomial fits to the binary diffusivity of each species
    /*!
     *  m_diffcoeff[ic] is vector of polynomial coefficients for species  i species  j
     *  that fits the binary diffusion coefficient. The relationship between i
     *  j and ic is determined from the following algorithm:
     *
     *      int ic = 0;
     *      for (i = 0; i < m_nsp; i++) {
     *         for (j = i; j < m_nsp; j++) {
     *           ic++;
     *         }
     *      }
     */
    std::vector<vector_fp> m_diffcoeffs;

    //! Matrix of binary diffusion coefficients at the reference pressure and
    //! the current temperature Size is nsp x nsp.
    DenseMatrix m_bdiff;
};

} // namespace Cantera

#endif
