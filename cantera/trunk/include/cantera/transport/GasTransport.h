/**
 * @file GasTransport.h
 */

#ifndef CT_GAS_TRANSPORT_H
#define CT_GAS_TRANSPORT_H

#include "TransportBase.h"

namespace Cantera {

//! Class GasTransport implements some functions and properties that are
//! shared by the MixTransport and MultiTransport classes.
class GasTransport : public Transport
{
public:
    virtual ~GasTransport() {}
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

    //! Vector of species mole fractions. These are processed so that all mole
    //! fractions are >= MIN_X. Length = m_kk.
    vector_fp m_molefracs;

    //! Internal storage for the viscosity of the mixture  (kg /m /s)
    doublereal m_viscmix;

    //! Update boolean for mixture rule for the mixture viscosity
    bool m_visc_ok;

    //! Update boolean for the weighting factors for the mixture viscosity
    bool m_viscwt_ok;

    //! Update boolean for the species viscosities
    bool m_spvisc_ok;

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
     *   m_wratjk(j,k)  = sqrt(mw[j]/mw[k])        j < k
     *   m_wratjk(k,j)  = sqrt(sqrt(mw[j]/mw[k]))  j < k
     */
    DenseMatrix m_wratjk;

    //! Holds square roots of molecular weight ratios
    /*!
     *   m_wratjk1(j,k)  = sqrt(1.0 + mw[k]/mw[j])        j < k
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
};

} // namespace Cantera

#endif
