/**
 *  @file ConstCpPoly.cpp
 * Declarations for the \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType \endlink object that
 * employs a constant heat capacity assumption (see \ref spthermo and
 * \link Cantera::ConstCpPoly ConstCpPoly \endlink).
 */
// Copyright 2001  California Institute of Technology

#include "ConstCpPoly.h"

namespace Cantera
{
ConstCpPoly::ConstCpPoly()
    :  m_t0(0.0),
       m_cp0_R(0.0),
       m_h0_R(0.0),
       m_s0_R(0.0),
       m_logt0(0.0)
{
}

ConstCpPoly::ConstCpPoly(size_t n, doublereal tlow, doublereal thigh,
                         doublereal pref,
                         const doublereal* coeffs) :
    SpeciesThermoInterpType(n, tlow, thigh, pref)
{
    m_t0 = coeffs[0];
    m_h0_R = coeffs[1]  / GasConstant;
    m_s0_R = coeffs[2]  / GasConstant;
    m_cp0_R = coeffs[3] / GasConstant;
    m_logt0 = log(m_t0);
}

ConstCpPoly::ConstCpPoly(const ConstCpPoly& b) :
    SpeciesThermoInterpType(b),
    m_t0(b.m_t0),
    m_cp0_R(b.m_cp0_R),
    m_h0_R(b.m_h0_R),
    m_s0_R(b.m_s0_R),
    m_logt0(b.m_logt0)
{
}

ConstCpPoly& ConstCpPoly::operator=(const ConstCpPoly& b)
{
    if (&b != this) {
        m_t0    = b.m_t0;
        m_cp0_R = b.m_cp0_R;
        m_h0_R  = b.m_h0_R;
        m_s0_R  = b.m_s0_R;
        m_logt0 = b.m_logt0;
        SpeciesThermoInterpType::operator=(b);
    }
    return *this;
}

SpeciesThermoInterpType*
ConstCpPoly::duplMyselfAsSpeciesThermoInterpType() const
{
    return new ConstCpPoly(*this);
}

void ConstCpPoly::updateProperties(const doublereal* tt,
                                   doublereal* cp_R,
                                   doublereal* h_RT,
                                   doublereal* s_R) const
{
    double t = *tt;
    doublereal logt = log(t);
    doublereal rt = 1.0/t;
    cp_R[m_index] = m_cp0_R;
    h_RT[m_index] = rt*(m_h0_R + (t - m_t0) * m_cp0_R);
    s_R[m_index]  = m_s0_R + m_cp0_R * (logt - m_logt0);
}

void ConstCpPoly::updatePropertiesTemp(const doublereal temp,
                                       doublereal* cp_R,
                                       doublereal* h_RT,
                                       doublereal* s_R) const
{
    doublereal logt = log(temp);
    doublereal rt = 1.0/temp;
    cp_R[m_index] = m_cp0_R;
    h_RT[m_index] = rt*(m_h0_R + (temp - m_t0) * m_cp0_R);
    s_R[m_index]  = m_s0_R + m_cp0_R * (logt - m_logt0);
}

void ConstCpPoly::reportParameters(size_t& n, int& type,
                                   doublereal& tlow, doublereal& thigh,
                                   doublereal& pref,
                                   doublereal* const coeffs) const
{
    n = m_index;
    type = CONSTANT_CP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = m_t0;
    coeffs[1] = m_h0_R * GasConstant;
    coeffs[2] = m_s0_R * GasConstant;
    coeffs[3] = m_cp0_R * GasConstant;
}

void ConstCpPoly::modifyParameters(doublereal* coeffs)
{
    m_t0 = coeffs[0];
    m_h0_R = coeffs[1]  / GasConstant;
    m_s0_R = coeffs[2]  / GasConstant;
    m_cp0_R = coeffs[3] / GasConstant;
    m_logt0 = log(m_t0);
}

#ifdef H298MODIFY_CAPABILITY

doublereal ConstCpPoly::reportHf298(doublereal* const h298) const
{
    double temp = 298.15;
    doublereal h = GasConstant * (m_h0_R + (temp - m_t0) * m_cp0_R);
    if (h298) {
        h298[m_index] = h;
    }
    return h;
}

void ConstCpPoly::modifyOneHf298(const size_t& k, const doublereal Hf298New)
{
    if (k != m_index) {
        return;
    }
    doublereal hnow = reportHf298();
    doublereal delH = Hf298New - hnow;
    m_h0_R += delH / GasConstant;
}

#endif

}
