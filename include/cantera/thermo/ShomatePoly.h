/**
 *  @file ShomatePoly.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the Shomate temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::ShomatePoly ShomatePoly\endlink and
 *   \link Cantera::ShomatePoly2 ShomatePoly2\endlink).
 *    Shomate polynomial expressions.
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_SHOMATEPOLY1_H
#define CT_SHOMATEPOLY1_H

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{
//! The Shomate polynomial parameterization for one temperature range for one
//! species
/*!
 * Seven coefficients \f$(A,\dots,G)\f$ are used to represent
 * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as
 * polynomials in the temperature, \f$ T \f$ :
 *
 * \f[
 * \tilde{c}_p^0(T) = A + B t + C t^2 + D t^3 + \frac{E}{t^2}
 * \f]
 * \f[
 * \tilde{h}^0(T) = A t + \frac{B t^2}{2} + \frac{C t^3}{3}
 *                + \frac{D t^4}{4}  - \frac{E}{t} + F.
 * \f]
 * \f[
 * \tilde{s}^0(T) = A\ln t + B t + \frac{C t^2}{2}
 *                + \frac{D t^3}{3} - \frac{E}{2t^2} + G.
 * \f]
 *
 * In the above expressions, the thermodynamic polynomials are expressed in
 * dimensional units, but the temperature,\f$ t \f$, is divided by 1000. The
 * following dimensions are assumed in the above expressions:
 *
 *    - \f$ \tilde{c}_p^0(T)\f$ = Heat Capacity (J/gmol*K)
 *    - \f$ \tilde{h}^0(T) \f$ = standard Enthalpy (kJ/gmol)
 *    - \f$ \tilde{s}^0(T) \f$= standard Entropy (J/gmol*K)
 *    - \f$ t \f$= temperature (K) / 1000.
 *
 * For more information about Shomate polynomials, see the NIST website,
 * http://webbook.nist.gov/
 *
 * Before being used within Cantera, the dimensions must be adjusted to those
 * used by Cantera (i.e., Joules and kmol).
 *
 * @ingroup spthermo
 */
class ShomatePoly : public SpeciesThermoInterpType
{
public:
    //! Empty constructor
    ShomatePoly() {}

    //! Normal constructor
    /*!
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients, [A,B,C,D,E,F,G], used to set
     *                     the parameters for the species standard state.
     *
     *  See the class description for the polynomial representation of the
     *  thermo functions in terms of \f$ A, \dots, G \f$.
     */
    ShomatePoly(double tlow, double thigh, double pref, const double* coeffs) :
        SpeciesThermoInterpType(tlow, thigh, pref),
        m_coeff(7)
    {
        for (size_t i = 0; i < 7; i++) {
            m_coeff[i] = coeffs[i] * 1000 / GasConstant;
        }
    }

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        return new ShomatePoly(*this);
    }

    virtual int reportType() const {
        return SHOMATE;
    }

    virtual size_t temperaturePolySize() const { return 6; }

    virtual void updateTemperaturePoly(double T, double* T_poly) const {
        doublereal tt = 1.e-3*T;
        T_poly[0] = tt;
        T_poly[1] = tt * tt;
        T_poly[2] = T_poly[1] * tt;
        T_poly[3] = 1.0/T_poly[1];
        T_poly[4] = std::log(tt);
        T_poly[5] = 1.0/tt;
    }

    /*!
     * @copydoc SpeciesThermoInterpType::updateProperties
     *
     * Form of the temperature polynomial:
     *   - `t` is T/1000.
     *   - `t[0] = t`
     *   - `t[1] = t*t`
     *   - `t[2] = t[1]*t`
     *   - `t[3] = 1.0/t[1]`
     *   - `t[4] = log(t)`
     *   - `t[5] = 1.0/t;
     */
    virtual void updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const {
        doublereal A = m_coeff[0];
        doublereal Bt = m_coeff[1]*tt[0];
        doublereal Ct2 = m_coeff[2]*tt[1];
        doublereal Dt3 = m_coeff[3]*tt[2];
        doublereal Etm2 = m_coeff[4]*tt[3];
        doublereal Ftm1 = m_coeff[5]*tt[5];
        doublereal G = m_coeff[6];

        *cp_R = A + Bt + Ct2 + Dt3 + Etm2;
        *h_RT = A + 0.5*Bt + 1.0/3.0*Ct2 + 0.25*Dt3 - Etm2 + Ftm1;
        *s_R = A*tt[4] + Bt + 0.5*Ct2 + 1.0/3.0*Dt3 - 0.5*Etm2 + G;
    }

    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R, doublereal* h_RT,
                                      doublereal* s_R) const {
        double tPoly[6];
        updateTemperaturePoly(temp, tPoly);
        updateProperties(tPoly, cp_R, h_RT, s_R);
    }

    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const {
        n = 0;
        type = SHOMATE;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        for (int i = 0; i < 7; i++) {
            coeffs[i] = m_coeff[i] * GasConstant / 1000;
        }
    }

    virtual void modifyParameters(doublereal* coeffs) {
        for (size_t i = 0; i < 7; i++) {
            m_coeff[i] = coeffs[i] * 1000 / GasConstant;
        }
    }

    virtual doublereal reportHf298(doublereal* const h298 = 0) const {
        double cp_R, h_RT, s_R;
        updatePropertiesTemp(298.15, &cp_R, &h_RT, &s_R);
        return h_RT * GasConstant * 298.15;
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New) {
        doublereal hnow = reportHf298();
        doublereal delH = Hf298New - hnow;
        m_coeff[5] += delH / (1e3 * GasConstant);
    }

protected:
    //! Array of coefficients
    vector_fp m_coeff;
};

//! The Shomate polynomial parameterization for two temperature ranges for one
//! species
/*!
 * Seven coefficients \f$(A,\dots,G)\f$ are used to represent
 * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as
 * polynomials in the temperature, \f$ T \f$, in one temperature region:
 *
 * \f[
 * \tilde{c}_p^0(T) = A + B t + C t^2 + D t^3 + \frac{E}{t^2}
 * \f]
 * \f[
 * \tilde{h}^0(T) = A t + \frac{B t^2}{2} + \frac{C t^3}{3}
 *                + \frac{D t^4}{4}  - \frac{E}{t}  + F.
 * \f]
 * \f[
 * \tilde{s}^0(T) = A\ln t + B t + \frac{C t^2}{2}
 *                + \frac{D t^3}{3} - \frac{E}{2t^2}  + G.
 * \f]
 *
 * In the above expressions, the thermodynamic polynomials are expressed
 * in dimensional units, but the temperature,\f$ t \f$, is divided by 1000. The
 * following dimensions are assumed in the above expressions:
 *
 *    - \f$ \tilde{c}_p^0(T)\f$ = Heat Capacity (J/gmol*K)
 *    - \f$ \tilde{h}^0(T) \f$ = standard Enthalpy (kJ/gmol)
 *    - \f$ \tilde{s}^0(T) \f$= standard Entropy (J/gmol*K)
 *    - \f$ t \f$= temperature (K) / 1000.
 *
 * For more information about Shomate polynomials, see the NIST website,
 * http://webbook.nist.gov/
 *
 * Before being used within Cantera, the dimensions must be adjusted to those
 * used by Cantera (i.e., Joules and kmol).
 *
 * This function uses two temperature regions, each with a Shomate polynomial
 * representation to represent the thermo functions. There are 15 coefficients,
 * therefore, in this representation. The first coefficient is the midrange
 * temperature.
 *
 * @ingroup spthermo
 */
class ShomatePoly2 : public SpeciesThermoInterpType
{
public:
    //! Empty constructor
    ShomatePoly2()
        : m_midT(0.0)
    {
        m_coeff.resize(15);
    }

    //! Normal constructor
    /*!
     * @param tlow    Minimum temperature
     * @param thigh   Maximum temperature
     * @param pref    reference pressure (Pa).
     * @param coeffs  Vector of coefficients used to set the parameters for the
     *                standard state. [Tmid, 7 low-T coeffs, 7 high-T coeffs]
     */
    ShomatePoly2(double tlow, double thigh, double pref, const double* coeffs) :
        SpeciesThermoInterpType(tlow, thigh, pref),
        m_midT(coeffs[0]),
        msp_low(tlow, coeffs[0], pref, coeffs+1),
        msp_high(coeffs[0], thigh, pref, coeffs+8),
        m_coeff(coeffs, coeffs + 15)
    {
    }

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        return new ShomatePoly2(*this);
    }

    virtual int reportType() const {
        return SHOMATE2;
    }

    virtual size_t temperaturePolySize() const { return 7; }

    virtual void updateTemperaturePoly(double T, double* T_poly) const {
        msp_low.updateTemperaturePoly(T, T_poly);
    }

    //! @copydoc ShomatePoly::updateProperties
    virtual void updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const {
        double T = 1000 * tt[0];
        if (T <= m_midT) {
            msp_low.updateProperties(tt, cp_R, h_RT, s_R);
        } else {
            msp_high.updateProperties(tt, cp_R, h_RT, s_R);
        }
    }

    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R,
                                      doublereal* h_RT,
                                      doublereal* s_R) const {
        if (temp <= m_midT) {
            msp_low.updatePropertiesTemp(temp, cp_R, h_RT, s_R);
        } else {
            msp_high.updatePropertiesTemp(temp, cp_R, h_RT, s_R);
        }
    }

    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const {
        n = 0;
        type = SHOMATE2;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        for (int i = 0; i < 15; i++) {
            coeffs[i] = m_coeff[i];
        }
    }

    //! Modify parameters for the standard state
    /*!
     * Here, we take the tact that we will just regenerate the object.
     *
     * @param coeffs   Vector of coefficients used to set the
     *                 parameters for the standard state.
     */
    virtual void modifyParameters(doublereal* coeffs) {
        std::copy(coeffs, coeffs + 15, m_coeff.begin());
        m_midT = coeffs[0];
        msp_low = ShomatePoly(m_lowT, m_midT, m_Pref, coeffs+1);
        msp_high = ShomatePoly(m_midT, m_highT, m_Pref, coeffs+8);
    }

    virtual doublereal reportHf298(doublereal* const h298 = 0) const {
        doublereal h;
        if (298.15 <= m_midT) {
            h = msp_low.reportHf298(h298);
        } else {
            h = msp_high.reportHf298(h298);
        }
        if (h298) {
            *h298 = h;
        }
        return h;
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New) {
        doublereal h298now = reportHf298(0);
        doublereal delH = Hf298New - h298now;
        double h = msp_low.reportHf298(0);
        double hnew = h + delH;
        msp_low.modifyOneHf298(k, hnew);
        h = msp_high.reportHf298(0);
        hnew = h + delH;
        msp_high.modifyOneHf298(k, hnew);
    }

protected:
    //! Midrange temperature (kelvin)
    doublereal m_midT;
    //! Shomate polynomial for the low temperature region.
    ShomatePoly msp_low;
    //! Shomate polynomial for the high temperature region.
    ShomatePoly msp_high;
    //! Array of the original coefficients.
    vector_fp m_coeff;
};
}

#endif
