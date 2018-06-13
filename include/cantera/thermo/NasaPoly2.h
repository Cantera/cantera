/**
 *  @file NasaPoly2.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the NASA temperature polynomial form applied to two temperature regions
 *  (see \ref spthermo and class \link Cantera::NasaPoly2 NasaPoly2\endlink).
 *
 *  Two zoned NASA polynomial parameterization
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_NASAPOLY2_H
#define CT_NASAPOLY2_H

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/NasaPoly1.h"

namespace Cantera
{
/**
 * The NASA polynomial parameterization for two temperature ranges. This
 * parameterization expresses the heat capacity as a fourth-order polynomial.
 * Note that this is the form used in the 1971 NASA equilibrium program and by
 * the Chemkin software package, but differs from the form used in the more
 * recent NASA equilibrium program.
 *
 * Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
 * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as
 * polynomials in \f$ T \f$ :
 * \f[
 * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 * \f[
 * \frac{h^0(T)}{RT} = a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2
 *                   + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4  + \frac{a_5}{T}.
 * \f]
 * \f[
 * \frac{s^0(T)}{R} = a_0\ln T + a_1 T + \frac{a_2}{2} T^2
 *                  + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4  + a_6.
 * \f]
 *
 * This class is designed specifically for use by the class MultiSpeciesThermo.
 *
 * @ingroup spthermo
 */
class NasaPoly2 : public SpeciesThermoInterpType
{
public:
    //! Full Constructor
    /*!
     * @param tlow      output - Minimum temperature
     * @param thigh     output - Maximum temperature
     * @param pref      output - reference pressure (Pa).
     * @param coeffs    Vector of coefficients used to set the parameters for
     *                  the standard state [Tmid, 7 high-T coeffs, 7 low-T
     *                  coeffs]. This is the coefficient order used in the
     *                  standard NASA format.
     */
    NasaPoly2(doublereal tlow, doublereal thigh, doublereal pref,
              const doublereal* coeffs) :
        SpeciesThermoInterpType(tlow, thigh, pref),
        m_midT(coeffs[0]),
        mnp_low(tlow, coeffs[0], pref, coeffs + 8),
        mnp_high(coeffs[0], thigh, pref, coeffs + 1),
        m_coeff(coeffs, coeffs + 15) {
    }

    virtual int reportType() const {
        return NASA2;
    }

    virtual size_t temperaturePolySize() const { return 6; }

    virtual void updateTemperaturePoly(double T, double* T_poly) const {
        mnp_low.updateTemperaturePoly(T, T_poly);
    }

    //! @copydoc NasaPoly1::updateProperties
    void updateProperties(const doublereal* tt,
                          doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {
        if (tt[0] <= m_midT) {
            mnp_low.updateProperties(tt, cp_R, h_RT, s_R);
        } else {
            mnp_high.updateProperties(tt, cp_R, h_RT, s_R);
        }
    }

    void updatePropertiesTemp(const doublereal temp,
                              doublereal* cp_R,
                              doublereal* h_RT,
                              doublereal* s_R) const {
        if (temp <= m_midT) {
            mnp_low.updatePropertiesTemp(temp, cp_R, h_RT, s_R);
        } else {
            mnp_high.updatePropertiesTemp(temp, cp_R, h_RT, s_R);
        }
    }

    void reportParameters(size_t& n, int& type,
                          doublereal& tlow, doublereal& thigh,
                          doublereal& pref,
                          doublereal* const coeffs) const {
        n = 0;
        type = NASA2;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        for (int i = 0; i < 15; i++) {
            coeffs[i] = m_coeff[i];
        }
    }

    doublereal reportHf298(doublereal* const h298 = 0) const {
        double h;
        if (298.15 <= m_midT) {
            h = mnp_low.reportHf298(0);
        } else {
            h = mnp_high.reportHf298(0);
        }
        if (h298) {
            *h298 = h;
        }
        return h;
    }

    void resetHf298() {
        mnp_low.resetHf298();
        mnp_high.resetHf298();
    }

    void modifyOneHf298(const size_t k, const doublereal Hf298New) {
        doublereal h298now = reportHf298(0);
        doublereal delH = Hf298New - h298now;
        double h = mnp_low.reportHf298(0);
        double hnew = h + delH;
        mnp_low.modifyOneHf298(k, hnew);
        h = mnp_high.reportHf298(0);
        hnew = h + delH;
        mnp_high.modifyOneHf298(k, hnew);
    }

    void validate(const std::string& name);

protected:
    //! Midrange temperature
    doublereal m_midT;
    //! NasaPoly1 object for the low temperature region.
    NasaPoly1 mnp_low;
    //! NasaPoly1 object for the high temperature region.
    NasaPoly1 mnp_high;
    //! array of polynomial coefficients
    vector_fp m_coeff;
};

}
#endif
