/**
 *  @file Nasa9PolyMultiTempRegion.cpp
 *  Definitions for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType
 *    SpeciesThermoInterpType\endlink  based
 *  on the NASA 9 coefficient temperature polynomial form
 *  applied to one temperature region
 *  (see \ref spthermo and class
 *   \link Cantera::Nasa9Poly1 Nasa9Poly1\endlink).
 *
 *  This parameterization has one NASA temperature region.
 */
// Copyright 2007  Sandia National Laboratories

#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/Nasa9PolyMultiTempRegion.h"

using namespace std;

namespace Cantera
{
Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion() :
    m_currRegion(0)
{
}

Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion(vector<Nasa9Poly1*>& regionPts) :
    m_currRegion(0)
{
    // From now on, we own these pointers
    for (Nasa9Poly1* region : regionPts) {
        m_regionPts.emplace_back(region);
    }
    m_lowerTempBounds.resize(regionPts.size());
    m_lowT = m_regionPts[0]->minTemp();
    m_highT = m_regionPts[m_regionPts.size()-1]->maxTemp();
    m_Pref = m_regionPts[0]->refPressure();
    for (size_t i = 0; i < m_regionPts.size(); i++) {
        m_lowerTempBounds[i] = m_regionPts[i]->minTemp();
        if (fabs(m_regionPts[i]->refPressure() - m_Pref) > 0.0001) {
            throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                               "refPressure inconsistency");
        }
        if (i > 0) {
            if (m_lowerTempBounds[i-1] >= m_lowerTempBounds[i]) {
                throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                                   "minTemp bounds inconsistency");
            }
            if (fabs(m_regionPts[i-1]->maxTemp() - m_lowerTempBounds[i]) > 0.0001) {
                throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                                   "Temp bounds inconsistency");
            }
        }
    }
}

Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion(const Nasa9PolyMultiTempRegion& b) :
    SpeciesThermoInterpType(b),
    m_lowerTempBounds(b.m_lowerTempBounds),
    m_currRegion(b.m_currRegion)
{
    m_regionPts.resize(b.m_regionPts.size());
    for (size_t i = 0; i < m_regionPts.size(); i++) {
        m_regionPts[i].reset(new Nasa9Poly1(*b.m_regionPts[i]));
    }
}

Nasa9PolyMultiTempRegion&
Nasa9PolyMultiTempRegion::operator=(const Nasa9PolyMultiTempRegion& b)
{
    if (&b != this) {
        SpeciesThermoInterpType::operator=(b);
        m_lowerTempBounds = b.m_lowerTempBounds;
        m_currRegion = b.m_currRegion;
        m_regionPts.resize(b.m_regionPts.size());
        for (size_t i = 0; i < m_regionPts.size(); i++) {
            m_regionPts[i].reset(new Nasa9Poly1(*b.m_regionPts[i]));
        }
    }
    return *this;
}

Nasa9PolyMultiTempRegion::~Nasa9PolyMultiTempRegion()
{
}

SpeciesThermoInterpType*
Nasa9PolyMultiTempRegion::duplMyselfAsSpeciesThermoInterpType() const
{
    return new Nasa9PolyMultiTempRegion(*this);
}

int Nasa9PolyMultiTempRegion::reportType() const
{
    return NASA9MULTITEMP;
}

void Nasa9PolyMultiTempRegion::updateTemperaturePoly(double T, double* T_poly) const
{
    T_poly[0] = T;
    T_poly[1] = T * T;
    T_poly[2] = T_poly[1] * T;
    T_poly[3] = T_poly[2] * T;
    T_poly[4] = 1.0 / T;
    T_poly[5] = T_poly[4] / T;
    T_poly[6] = std::log(T);
}

void Nasa9PolyMultiTempRegion::updateProperties(const doublereal* tt,
        doublereal* cp_R,
        doublereal* h_RT,
        doublereal* s_R) const
{
    m_currRegion = 0;
    for (size_t i = 1; i < m_regionPts.size(); i++) {
        if (tt[0] < m_lowerTempBounds[i]) {
            break;
        }
        m_currRegion++;
    }

    m_regionPts[m_currRegion]->updateProperties(tt, cp_R, h_RT, s_R);
}

void Nasa9PolyMultiTempRegion::updatePropertiesTemp(const doublereal temp,
        doublereal* cp_R, doublereal* h_RT,
        doublereal* s_R) const
{
    // Now find the region
    m_currRegion = 0;
    for (size_t i = 1; i < m_regionPts.size(); i++) {
        if (temp < m_lowerTempBounds[i]) {
            break;
        }
        m_currRegion++;
    }

    m_regionPts[m_currRegion]->updatePropertiesTemp(temp, cp_R, h_RT, s_R);
}

void Nasa9PolyMultiTempRegion::reportParameters(size_t& n, int& type,
        doublereal& tlow, doublereal& thigh,
        doublereal& pref,
        doublereal* const coeffs) const
{
    n = 0;
    type = NASA9MULTITEMP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    double ctmp[12];
    coeffs[0] = double(m_regionPts.size());
    int index = 1;
    size_t n_tmp = 0;
    int type_tmp = 0;
    double pref_tmp = 0.0;
    for (size_t iReg = 0; iReg < m_regionPts.size(); iReg++) {
        m_regionPts[iReg]->reportParameters(n_tmp, type_tmp,
                                            coeffs[index], coeffs[index+1],
                                            pref_tmp, ctmp);
        for (int i = 0; i < 9; i++) {
            coeffs[index+2+i] = ctmp[3+i];
        }
        index += 11;
    }
}

void Nasa9PolyMultiTempRegion::modifyParameters(doublereal* coeffs)
{
    int index = 3;
    for (size_t iReg = 0; iReg < m_regionPts.size(); iReg++) {
        m_regionPts[iReg]->modifyParameters(coeffs + index);
        index += 11;
    }
}

}
