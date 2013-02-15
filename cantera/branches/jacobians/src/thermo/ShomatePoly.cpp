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
#include "ShomatePoly.h"

namespace Cantera {

// Update the properties for this species, given a temperature polynomial
/*
 * This method is called with a pointer to an array containing the functions of
 * temperature needed by this  parameterization, and three pointers to arrays where the
 * computed property values should be written. This method updates only one value in
 * each array.
 *
 *  tt is T/1000.
 *  m_t[0] = tt;
 *  m_t[1] = tt*tt;
 *  m_t[2] = m_t[1]*tt;
 *  m_t[3] = 1.0/m_t[1];
 *  m_t[4] = log(tt);
 *  m_t[5] = 1.0/GasConstant;
 *  m_t[6] = 1.0/(GasConstant * T);
 *
 * @param tt      Vector of temperature polynomials
 * @param cp_R    Vector of Dimensionless heat capacities.
 *                (length m_kk).
 * @param h_RT    Vector of Dimensionless enthalpies.
 *                (length m_kk).
 * @param s_R     Vector of Dimensionless entropies.
 *                (length m_kk).
 */
template<typename ValAndDerivType>
void ShomatePoly<ValAndDerivType>::updateProperties(const ValAndDerivType* tt, ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                                     ValAndDerivType* s_R) const
{

    doublereal A = m_coeff[0];
    doublereal Bt = m_coeff[1] * tt[0];
    doublereal Ct2 = m_coeff[2] * tt[1];
    doublereal Dt3 = m_coeff[3] * tt[2];
    doublereal Etm2 = m_coeff[4] * tt[3];
    doublereal F = m_coeff[5];
    doublereal G = m_coeff[6];

    doublereal cp, h, s;
    cp = A + Bt + Ct2 + Dt3 + Etm2;
    h = tt[0] * (A + 0.5 * Bt + OneThird * Ct2 + 0.25 * Dt3 - Etm2) + F;
    s = A * tt[4] + Bt + 0.5 * Ct2 + OneThird * Dt3 - 0.5 * Etm2 + G;

    /*
     *  Shomate polynomials parameterizes assuming units of
     *  J/(gmol*K) for cp_r and s_R and kJ/(gmol) for h.
     *  However, Cantera assumes default MKS units of
     *  J/(kmol*K). This requires us to multiply cp and s
     *  by 1.e3 and h by 1.e6, before we then nondimensionlize
     *  the results by dividing by (GasConstant * T),
     *  where GasConstant has units of J/(kmol * K).
     */
    cp_R[m_index] = 1.e3 * cp * tt[5];
    h_RT[m_index] = 1.e6 * h * tt[6];
    s_R[m_index] = 1.e3 * s * tt[5];
}


template<>
void ShomatePoly<doubleFAD>::updateProperties(const doubleFAD* tt, doubleFAD* cp_R, doubleFAD* h_RT, doubleFAD* s_R) const
{

    doublereal A = m_coeff[0];
    doubleFAD Bt = m_coeff[1] * tt[0];
    doubleFAD Ct2 = m_coeff[2] * tt[1];
    doubleFAD Dt3 = m_coeff[3] * tt[2];
    doubleFAD Etm2 = m_coeff[4] * tt[3];
    doublereal F = m_coeff[5];
    doublereal G = m_coeff[6];

    doubleFAD cp, h, s;
    cp = A + Bt + Ct2 + Dt3 + Etm2;
    h = tt[0] * (A + 0.5 * Bt + OneThird * Ct2 + 0.25 * Dt3 - Etm2) + F;
    s = A * tt[4] + Bt + 0.5 * Ct2 + OneThird * Dt3 - 0.5 * Etm2 + G;

    /*
     *  Shomate polynomials parameterizes assuming units of
     *  J/(gmol*K) for cp_r and s_R and kJ/(gmol) for h.
     *  However, Cantera assumes default MKS units of
     *  J/(kmol*K). This requires us to multiply cp and s
     *  by 1.e3 and h by 1.e6, before we then nondimensionlize
     *  the results by dividing by (GasConstant * T),
     *  where GasConstant has units of J/(kmol * K).
     */
    cp_R[m_index] = 1.e3 * cp * tt[5];
    h_RT[m_index] = 1.e6 * h * tt[6];
    s_R[m_index] = 1.e3 * s * tt[5];
}

//==================================================================================================================================

// Update the properties for this species, given a temperature polynomial
/*
 * This method is called with a pointer to an array containing the functions of
 * temperature needed by this  parameterization, and three pointers to arrays where the
 * computed property values should be written. This method updates only one value in
 * each array.
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
 * @param cp_R    Vector of Dimensionless heat capacities.
 *                (length m_kk).
 * @param h_RT    Vector of Dimensionless enthalpies.
 *                (length m_kk).
 * @param s_R     Vector of Dimensionless entropies.
 *                (length m_kk).
 */
template<typename ValAndDerivType>
void ShomatePoly2<ValAndDerivType>::updateProperties(const ValAndDerivType* tt, ValAndDerivType* cp_R, ValAndDerivType* h_RT,
                                                     ValAndDerivType* s_R) const
{
    double T = 1000 * tt[0];
    if (T <= m_midT) {
        msp_low->updateProperties(tt, cp_R, h_RT, s_R);
    } else {
        msp_high->updateProperties(tt, cp_R, h_RT, s_R);
    }

}

template<>
void ShomatePoly2<doubleFAD>::updateProperties(const doubleFAD* tt, doubleFAD* cp_R, doubleFAD* h_RT, doubleFAD* s_R) const
{
    double T = 1000 * tt[0].val();
    if (T <= m_midT) {
        msp_low->updateProperties(tt, cp_R, h_RT, s_R);
    } else {
        msp_high->updateProperties(tt, cp_R, h_RT, s_R);
    }

}

template class ShomatePoly<doublereal> ;

#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class ShomatePoly<doubleFAD> ;
#endif
#endif


template class ShomatePoly2<doublereal> ;

#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class ShomatePoly2<doubleFAD> ;
#endif
#endif

}
