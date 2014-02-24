
/**
 *  @file NasaPoly1.h
 *  Header for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on the NASA temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::NasaPoly1 NasaPoly1\endlink).
 *
 *  This parameterization has one NASA temperature region.
 */

#ifndef CT_NASAPOLY1_H
#define CT_NASAPOLY1_H
// Copyright 2001  California Institute of Technology

#include "cantera/base/global.h"
#include "SpeciesThermoInterpType.h"
#include <iostream>

namespace Cantera
{
/**
 * The NASA polynomial parameterization for one temperature range.
 * This parameterization expresses the heat capacity as a
 * fourth-order polynomial. Note that this is the form used in the
 * 1971 NASA equilibrium program and by the Chemkin software
 * package, but differs from the form used in the more recent NASA
 * equilibrium program.
 *
 * Seven coefficients \f$(a_0,\dots,a_6)\f$ are used to represent
 * \f$ c_p^0(T)\f$, \f$ h^0(T)\f$, and \f$ s^0(T) \f$ as
 * polynomials in \f$ T \f$ :
 * \f[
 * \frac{c_p(T)}{R} = a_0 + a_1 T + a_2 T^2 + a_3 T^3 + a_4 T^4
 * \f]
 * \f[
 * \frac{h^0(T)}{RT} = a_0 + \frac{a_1}{2} T + \frac{a_2}{3} T^2
 * + \frac{a_3}{4} T^3 + \frac{a_4}{5} T^4  + \frac{a_5}{T}.
 * \f]
 * \f[
 * \frac{s^0(T)}{R} = a_0\ln T + a_1 T + \frac{a_2}{2} T^2
 + \frac{a_3}{3} T^3 + \frac{a_4}{4} T^4  + a_6.
 * \f]
 *
 * This class is designed specifically for use by class NasaThermo.
 * @ingroup spthermo
 */
class NasaPoly1 : public SpeciesThermoInterpType
{
public:
    //! Empty constructor
    NasaPoly1()
        : m_coeff(7, 0.0) {}

    //! constructor used in templated instantiations
    /*!
     * @param n            Species index
     * @param tlow         Minimum temperature
     * @param thigh        Maximum temperature
     * @param pref         reference pressure (Pa).
     * @param coeffs       Vector of coefficients used to set the
     *                     parameters for the standard state.
     */
    NasaPoly1(size_t n, doublereal tlow, doublereal thigh, doublereal pref,
              const doublereal* coeffs) :
        SpeciesThermoInterpType(n, tlow, thigh, pref),
        m_coeff(vector_fp(7)) {
        std::copy(coeffs, coeffs + 7, m_coeff.begin());
    }

    //! copy constructor
    /*!
     * @param b object to be copied
     */
    NasaPoly1(const NasaPoly1& b) :
        SpeciesThermoInterpType(b),
        m_coeff(vector_fp(7))
    {
        std::copy(b.m_coeff.begin(),
                  b.m_coeff.begin() + 7,
                  m_coeff.begin());
    }

    //! assignment operator
    /*!
     * @param b object to be copied
     */
    NasaPoly1& operator=(const NasaPoly1& b) {
        if (&b != this) {
            SpeciesThermoInterpType::operator=(b);
            std::copy(b.m_coeff.begin(),
                      b.m_coeff.begin() + 7,
                      m_coeff.begin());
        }
        return *this;
    }

    virtual SpeciesThermoInterpType*
    duplMyselfAsSpeciesThermoInterpType() const {
        NasaPoly1* np = new NasaPoly1(*this);
        return (SpeciesThermoInterpType*) np;
    }

    virtual int reportType() const {
        return NASA1;
    }

    //! Update the properties for this species, given a temperature polynomial
    /*!
     * This method is called with a pointer to an array containing the
     * functions of temperature needed by this  parameterization, and three
     * pointers to arrays where the computed property values should be
     * written. This method updates only one value in each array.
     *
     * Temperature Polynomial:
     *  tt[0] = t;
     *  tt[1] = t*t;
     *  tt[2] = m_t[1]*t;
     *  tt[3] = m_t[2]*t;
     *  tt[4] = 1.0/t;
     *  tt[5] = std::log(t);
     *
     * @param tt      vector of temperature polynomials
     * @param cp_R    Vector of Dimensionless heat capacities. (length m_kk).
     * @param h_RT    Vector of Dimensionless enthalpies. (length m_kk).
     * @param s_R     Vector of Dimensionless entropies. (length m_kk).
     */
    virtual void updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT, doublereal* s_R) const {
        doublereal ct0 = m_coeff[2];          // a0
        doublereal ct1 = m_coeff[3]*tt[0];    // a1 * T
        doublereal ct2 = m_coeff[4]*tt[1];    // a2 * T^2
        doublereal ct3 = m_coeff[5]*tt[2];    // a3 * T^3
        doublereal ct4 = m_coeff[6]*tt[3];    // a4 * T^4

        doublereal cp, h, s;
        cp = ct0 + ct1 + ct2 + ct3 + ct4;
        h = ct0 + 0.5*ct1 + OneThird*ct2 + 0.25*ct3 + 0.2*ct4
            + m_coeff[0]*tt[4];               // last term is a5/T
        s = ct0*tt[5] + ct1 + 0.5*ct2 + OneThird*ct3
            +0.25*ct4 + m_coeff[1];           // last term is a6

        // return the computed properties in the location in the output
        // arrays for this species
        cp_R[m_index] = cp;
        h_RT[m_index] = h;
        s_R[m_index] = s;
        //writelog("NASA1: for species "+int2str(m_index)+", h_RT = "+
        //    fp2str(h)+"\n");
    }

    virtual void updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R, doublereal* h_RT,
                                      doublereal* s_R) const {
        double tPoly[6];
        tPoly[0]  = temp;
        tPoly[1]  = temp * temp;
        tPoly[2]  = tPoly[1] * temp;
        tPoly[3]  = tPoly[2] * temp;
        tPoly[4]  = 1.0 / temp;
        tPoly[5]  = std::log(temp);
        updateProperties(tPoly, cp_R, h_RT, s_R);
    }

    //! @deprecated
    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const {
        warn_deprecated("NasaPoly1::reportParameters");
        n = m_index;
        type = NASA1;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        coeffs[5] = m_coeff[0];
        coeffs[6] = m_coeff[1];
        for (int i = 2; i < 7; i++) {
            coeffs[i-2] = m_coeff[i];
        }
    }

    //! Modify parameters for the standard state
    /*!
     * @param coeffs   Vector of coefficients used to set the
     *                 parameters for the standard state.
     * @deprecated
     */
    virtual void modifyParameters(doublereal* coeffs) {
        warn_deprecated("NasaPoly1::modifyParameters");
        m_coeff[0] = coeffs[5];
        m_coeff[1] = coeffs[6];
        for (int i = 0; i < 5; i++) {
            m_coeff[i+2] = coeffs[i];
        }
    }

#ifdef H298MODIFY_CAPABILITY

    virtual doublereal reportHf298(doublereal* const h298 = 0) const {
        double tt[6];
        double temp = 298.15;
        tt[0]  = temp;
        tt[1]  = temp * temp;
        tt[2]  = tt[1] * temp;
        tt[3]  = tt[2] * temp;
        tt[4]  = 1.0 / temp;
        //tt[5]  = std::log(temp);
        doublereal ct0 = m_coeff[2];          // a0
        doublereal ct1 = m_coeff[3]*tt[0];    // a1 * T
        doublereal ct2 = m_coeff[4]*tt[1];    // a2 * T^2
        doublereal ct3 = m_coeff[5]*tt[2];    // a3 * T^3
        doublereal ct4 = m_coeff[6]*tt[3];    // a4 * T^4

        double h_RT = ct0 + 0.5*ct1 + OneThird*ct2 + 0.25*ct3 + 0.2*ct4
                      + m_coeff[0]*tt[4];               // last t

        double h = h_RT * GasConstant * temp;
        if (h298) {
            h298[m_index] = h;
        }
        return h;
    }

    virtual void modifyOneHf298(const size_t& k, const doublereal Hf298New) {
        if (k != m_index) {
            return;
        }
        double hcurr = reportHf298(0);
        double delH = Hf298New - hcurr;
        m_coeff[0] += (delH) / GasConstant;
    }

#endif

protected:
    //! array of polynomial coefficients
    vector_fp m_coeff;
};

}
#endif
