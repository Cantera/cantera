/**
 * @file Nasa9Poly1.cpp Definitions for a single-species standard state object
 *     derived from
 * \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink based
 * on the NASA 9 coefficient temperature polynomial form applied to one
 * temperature region (see \ref spthermo and class \link Cantera::Nasa9Poly1
 * Nasa9Poly1\endlink).
 *
 * This parameterization has one NASA temperature region.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Nasa9Poly1.h"

namespace Cantera
{

Nasa9Poly1::Nasa9Poly1()
    : m_coeff(9)
{
}

Nasa9Poly1::Nasa9Poly1(double tlow, double thigh, double pref,
                       const double* coeffs) :
    SpeciesThermoInterpType(tlow, thigh, pref),
    m_coeff(coeffs, coeffs + 9)
{
}

void Nasa9Poly1::setParameters(const vector_fp &coeffs)
{
    if (coeffs.size() != 9) {
        throw CanteraError("Nasa9Poly1::setParameters", "Array must contain "
            "9 coefficients, but {} were given.", coeffs.size());
    }
    m_coeff = coeffs;
}

int Nasa9Poly1::reportType() const
{
    return NASA9;
}

void Nasa9Poly1::updateTemperaturePoly(double T, double* T_poly) const
{
    T_poly[0] = T;
    T_poly[1] = T * T;
    T_poly[2] = T_poly[1] * T;
    T_poly[3] = T_poly[2] * T;
    T_poly[4] = 1.0 / T;
    T_poly[5] = T_poly[4] / T;
    T_poly[6] = std::log(T);
}

void Nasa9Poly1::updateProperties(const doublereal* tt,
                                  doublereal* cp_R, doublereal* h_RT,
                                  doublereal* s_R) const
{

    doublereal ct0 = m_coeff[0] * tt[5]; // a0 / (T^2)
    doublereal ct1 = m_coeff[1] * tt[4]; // a1 / T
    doublereal ct2 = m_coeff[2]; // a2
    doublereal ct3 = m_coeff[3] * tt[0]; // a3 * T
    doublereal ct4 = m_coeff[4] * tt[1]; // a4 * T^2
    doublereal ct5 = m_coeff[5] * tt[2]; // a5 * T^3
    doublereal ct6 = m_coeff[6] * tt[3]; // a6 * T^4

    doublereal cpdivR = ct0 + ct1 + ct2 + ct3 + ct4 + ct5 + ct6;
    doublereal hdivRT = -ct0 + tt[6]*ct1 + ct2 + 0.5*ct3 + 1.0/3.0*ct4
                        + 0.25*ct5 + 0.2*ct6 + m_coeff[7] * tt[4];
    doublereal sdivR = -0.5*ct0 - ct1 + tt[6]*ct2 + ct3 + 0.5*ct4
                       + 1.0/3.0*ct5 + 0.25*ct6 + m_coeff[8];

    // return the computed properties for this species
    *cp_R = cpdivR;
    *h_RT = hdivRT;
    *s_R = sdivR;
}

void Nasa9Poly1::updatePropertiesTemp(const doublereal temp,
                                      doublereal* cp_R, doublereal* h_RT,
                                      doublereal* s_R) const
{
    double tPoly[7];
    updateTemperaturePoly(temp, tPoly);
    updateProperties(tPoly, cp_R, h_RT, s_R);
}

void Nasa9Poly1::reportParameters(size_t& n, int& type,
                                  doublereal& tlow, doublereal& thigh,
                                  doublereal& pref,
                                  doublereal* const coeffs) const
{
    n = 0;
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

}
