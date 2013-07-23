/**
 *  @file Nasa9Poly1.cpp
 *  Definitions for a single-species standard state object derived
 *  from
 *  \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink
 *  based
 *  on the NASA 9 coefficient temperature polynomial form applied to one temperature region
 *  (see \ref spthermo and class \link Cantera::Nasa9Poly1 Nasa9Poly1\endlink).
 *
 *  This parameterization has one NASA temperature region.
 */
// Copyright 2007  Sandia National Laboratories

#include "cantera/thermo/Nasa9Poly1.h"

namespace Cantera
{
Nasa9Poly1::Nasa9Poly1()
    : m_coeff(vector_fp(9))
{
    m_Pref = 1.0e5;
}

Nasa9Poly1::Nasa9Poly1(size_t n, doublereal tlow, doublereal thigh,
                       doublereal pref,
                       const doublereal* coeffs) :
    SpeciesThermoInterpType(n, tlow, thigh, pref),
    m_coeff(vector_fp(9))
{
    std::copy(coeffs, coeffs + 9, m_coeff.begin());
}

Nasa9Poly1::Nasa9Poly1(const Nasa9Poly1& b) :
    SpeciesThermoInterpType(b),
    m_coeff(vector_fp(9))
{
    std::copy(b.m_coeff.begin(),
              b.m_coeff.begin() + 9,
              m_coeff.begin());
}

Nasa9Poly1& Nasa9Poly1::operator=(const Nasa9Poly1& b)
{
    if (&b != this) {
        SpeciesThermoInterpType::operator=(b);
        std::copy(b.m_coeff.begin(),
                  b.m_coeff.begin() + 9,
                  m_coeff.begin());
    }
    return *this;
}

SpeciesThermoInterpType*
Nasa9Poly1::duplMyselfAsSpeciesThermoInterpType() const
{
    return new Nasa9Poly1(*this);
}

int Nasa9Poly1::reportType() const
{
    return NASA9;
}

void Nasa9Poly1::updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const
{

    doublereal ct0 = m_coeff[0] * tt[5];   // a0 / (T^2)
    doublereal ct1 = m_coeff[1] * tt[4];   // a1 / T
    doublereal ct2 = m_coeff[2];           // a2
    doublereal ct3 = m_coeff[3] * tt[0];   // a3 * T
    doublereal ct4 = m_coeff[4] * tt[1];   // a4 * T^2
    doublereal ct5 = m_coeff[5] * tt[2];   // a5 * T^3
    doublereal ct6 = m_coeff[6] * tt[3];   // a6 * T^4


    doublereal cpdivR = ct0 + ct1 + ct2 + ct3 + ct4 + ct5 + ct6;
    doublereal hdivRT = -ct0 + tt[6]*ct1  + ct2 + 0.5*ct3 + OneThird*ct4
                        + 0.25*ct5  + 0.2*ct6 + m_coeff[7] * tt[4];
    doublereal sdivR  = -0.5*ct0  - ct1 + tt[6]*ct2  + ct3  + 0.5*ct4
                        + OneThird*ct5 + 0.25*ct6 + m_coeff[8];

    // return the computed properties in the location in the output
    // arrays for this species
    cp_R[m_index] = cpdivR;
    h_RT[m_index] = hdivRT;
    s_R[m_index] = sdivR;
    //writelog("NASA9poly1: for species "+int2str(m_index)+", h_RT = "+
    //    fp2str(h)+"\n");
}

void Nasa9Poly1::updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R, doublereal* h_RT,
                                      doublereal* s_R) const
{
    double tPoly[7];
    tPoly[0]  = temp;
    tPoly[1]  = temp * temp;
    tPoly[2]  = tPoly[1] * temp;
    tPoly[3]  = tPoly[2] * temp;
    tPoly[4]  = 1.0 / temp;
    tPoly[5]  = tPoly[4] / temp;
    tPoly[6]  = std::log(temp);
    updateProperties(tPoly, cp_R, h_RT, s_R);
}

void Nasa9Poly1::reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const
{
    warn_deprecated("Nasa9Poly1::reportParameters");
    n = m_index;
    type = NASA9;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = 1;
    coeffs[1] = m_lowT;
    coeffs[2] = m_highT;
    for (int i = 0; i < 9; i++) {
        coeffs[i+3] = m_coeff[i];
    }
}

void Nasa9Poly1::modifyParameters(doublereal* coeffs)
{
    warn_deprecated("Nasa9Poly1::modifyParameters");
    for (int i = 0; i < 9; i++) {
        m_coeff[i] = coeffs[i];
    }
}

}
