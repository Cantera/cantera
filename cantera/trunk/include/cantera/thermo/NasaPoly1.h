
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

#include "SpeciesThermoInterpType.h"

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
     * @param n       Species index
     * @param tlow    Minimum temperature
     * @param thigh   Maximum temperature
     * @param pref    reference pressure (Pa).
     * @param coeffs  Vector of coefficients used to set the parameters for the
     *                standard state, in the order [a5,a6,a0,a1,a2,a3,a4]
     * @deprecated  Use the constructor which does not take species index. To be
     *     removed after Cantera 2.2.
     */
    NasaPoly1(size_t n, doublereal tlow, doublereal thigh, doublereal pref,
              const doublereal* coeffs) :
        SpeciesThermoInterpType(n, tlow, thigh, pref),
        m_coeff(7)
    {
        for (size_t i = 0; i < 5; i++) {
            m_coeff[i] = coeffs[i+2];
        }
        m_coeff[5] = coeffs[0];
        m_coeff[6] = coeffs[1];
    }

    //! Normal constructor
    /*!
     * @param tlow    Minimum temperature
     * @param thigh   Maximum temperature
     * @param pref    reference pressure (Pa).
     * @param coeffs  Vector of coefficients used to set the parameters for the
     *                standard state, in the order [a0,a1,a2,a3,a4,a5,a6]
     */
    NasaPoly1(double tlow, double thigh, double pref, const double* coeffs)
        : SpeciesThermoInterpType(tlow, thigh, pref)
        , m_coeff(coeffs, coeffs+7)
    {
    }

    //! copy constructor
    /*!
     * @param b object to be copied
     */
    NasaPoly1(const NasaPoly1& b) :
        SpeciesThermoInterpType(b),
        m_coeff(b.m_coeff)
    {
    }

    //! assignment operator
    /*!
     * @param b object to be copied
     */
    NasaPoly1& operator=(const NasaPoly1& b) {
        if (&b != this) {
            SpeciesThermoInterpType::operator=(b);
            m_coeff = b.m_coeff;
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

    virtual size_t temperaturePolySize() const { return 6; }

    virtual void updateTemperaturePoly(double T, double* T_poly) const {
        T_poly[0] = T;
        T_poly[1] = T * T;
        T_poly[2] = T_poly[1] * T;
        T_poly[3] = T_poly[2] * T;
        T_poly[4] = 1.0 / T;
        T_poly[5] = std::log(T);
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
        doublereal ct0 = m_coeff[0];          // a0
        doublereal ct1 = m_coeff[1]*tt[0];    // a1 * T
        doublereal ct2 = m_coeff[2]*tt[1];    // a2 * T^2
        doublereal ct3 = m_coeff[3]*tt[2];    // a3 * T^3
        doublereal ct4 = m_coeff[4]*tt[3];    // a4 * T^4

        doublereal cp, h, s;
        cp = ct0 + ct1 + ct2 + ct3 + ct4;
        h = ct0 + 0.5*ct1 + 1.0/3.0*ct2 + 0.25*ct3 + 0.2*ct4
            + m_coeff[5]*tt[4];               // last term is a5/T
        s = ct0*tt[5] + ct1 + 0.5*ct2 + 1.0/3.0*ct3
            +0.25*ct4 + m_coeff[6];           // last term is a6

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
        updateTemperaturePoly(temp, tPoly);
        updateProperties(tPoly, cp_R, h_RT, s_R);
    }

    virtual void reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const {
        n = m_index;
        type = NASA1;
        tlow = m_lowT;
        thigh = m_highT;
        pref = m_Pref;
        std::copy(m_coeff.begin(), m_coeff.end(), coeffs);
    }

    //! Modify parameters for the standard state
    /*!
     * @param coeffs   Vector of coefficients used to set the
     *                 parameters for the standard state.
     */
    virtual void modifyParameters(doublereal* coeffs) {
        std::copy(coeffs, coeffs+7, m_coeff.begin());
    }

    virtual doublereal reportHf298(doublereal* const h298 = 0) const {
        double tt[6];
        double temp = 298.15;
        updateTemperaturePoly(temp, tt);
        doublereal ct0 = m_coeff[0];          // a0
        doublereal ct1 = m_coeff[1]*tt[0];    // a1 * T
        doublereal ct2 = m_coeff[2]*tt[1];    // a2 * T^2
        doublereal ct3 = m_coeff[3]*tt[2];    // a3 * T^3
        doublereal ct4 = m_coeff[4]*tt[3];    // a4 * T^4

        double h_RT = ct0 + 0.5*ct1 + 1.0/3.0*ct2 + 0.25*ct3 + 0.2*ct4
                      + m_coeff[5]*tt[4];               // last t

        double h = h_RT * GasConstant * temp;
        if (h298) {
            h298[m_index] = h;
        }
        return h;
    }

    virtual void modifyOneHf298(const size_t k, const doublereal Hf298New) {
        if (k != m_index) {
            return;
        }
        double hcurr = reportHf298(0);
        double delH = Hf298New - hcurr;
        m_coeff[5] += (delH) / GasConstant;
    }

protected:
    //! array of polynomial coefficients, stored in the order [a0, ..., a6]
    vector_fp m_coeff;
};

}
#endif
