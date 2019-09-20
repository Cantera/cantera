/**
 *  @file ShomatePoly.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the Shomate temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::ShomatePoly ShomatePoly\endlink and
 *   \link Cantera::ShomatePoly2 ShomatePoly2\endlink).
 *    Shomate polynomial expressions.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SHOMATEPOLY1_H
#define CT_SHOMATEPOLY1_H

#include "cantera/thermo/SpeciesThermoInterpType.h"

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
    ShomatePoly() : m_coeff(7), m_coeff5_orig(0.0) {}

    //! Constructor with all input data
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
        m_coeff5_orig = m_coeff[5];
    }

    //! Set array of 7 polynomial coefficients. Input values are assumed to be
    //! on a kJ/mol basis.
    void setParameters(const vector_fp& coeffs) {
        if (coeffs.size() != 7) {
            throw CanteraError("ShomatePoly::setParameters", "Array must "
                "contain 7 coefficients, but {} were given.", coeffs.size());
        }
        for (size_t i = 0; i < 7; i++) {
            m_coeff[i] = coeffs[i] * 1000 / GasConstant;
        }
        m_coeff5_orig = m_coeff[5];
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

    virtual void resetHf298() {
        m_coeff[5] = m_coeff5_orig;
    }

protected:
    //! Array of coefficients
    vector_fp m_coeff;
    double m_coeff5_orig;
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
    ShomatePoly2() : m_midT(0.0) {}

    //! Constructor with all input data
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
        msp_high(coeffs[0], thigh, pref, coeffs+8)
    {
    }

    virtual void setMinTemp(double Tmin) {
        SpeciesThermoInterpType::setMinTemp(Tmin);
        msp_low.setMinTemp(Tmin);
    }

    virtual void setMaxTemp(double Tmax) {
        SpeciesThermoInterpType::setMaxTemp(Tmax);
        msp_high.setMaxTemp(Tmax);
    }

    virtual void setRefPressure(double Pref) {
        SpeciesThermoInterpType::setRefPressure(Pref);
        msp_low.setRefPressure(Pref);
        msp_high.setRefPressure(Pref);
    }

    /*!
     * @param Tmid  Temperature [K] at the boundary between the low and high
     *              temperature polynomials
     * @param low   Vector of 7 coefficients for the low temperature polynomial
     * @param high  Vector of 7 coefficients for the high temperature polynomial
     */
    void setParameters(double Tmid, const vector_fp& low, const vector_fp& high) {
        m_midT = Tmid;
        msp_low.setMaxTemp(Tmid);
        msp_high.setMinTemp(Tmid);
        msp_low.setParameters(low);
        msp_high.setParameters(high);
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

    virtual size_t nCoeffs() const { return 15; }

    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const {
        msp_low.reportParameters(n, type, tlow, coeffs[0], pref, coeffs + 1);
        msp_high.reportParameters(n, type, coeffs[0], thigh, pref, coeffs + 8);
        type = SHOMATE2;
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

    virtual void resetHf298() {
        msp_low.resetHf298();
        msp_high.resetHf298();
    }

protected:
    //! Midrange temperature (kelvin)
    doublereal m_midT;
    //! Shomate polynomial for the low temperature region.
    ShomatePoly msp_low;
    //! Shomate polynomial for the high temperature region.
    ShomatePoly msp_high;
};
}

#endif
