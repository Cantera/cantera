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

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"
#include "Nasa9PolyMultiTempRegion.h"

using namespace std;

namespace Cantera
{
Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion() :
    m_numTempRegions(0),
    m_currRegion(0)
{
}

Nasa9PolyMultiTempRegion::
Nasa9PolyMultiTempRegion(std::vector<Cantera::Nasa9Poly1*> &regionPts) :
    m_numTempRegions(0),
    m_currRegion(0)
{
    m_numTempRegions = regionPts.size();
    // Do a shallow copy of the pointers. From now on, we will
    // own these pointers and be responsible for deleting them.
    m_regionPts = regionPts;
    m_lowerTempBounds.resize(m_numTempRegions);
    m_lowT = m_regionPts[0]->minTemp();
    m_highT = m_regionPts[m_numTempRegions-1]->maxTemp();
    m_Pref = m_regionPts[0]->refPressure();
    m_index = m_regionPts[0]->speciesIndex();
    for (size_t i = 0; i < m_numTempRegions; i++) {
        m_lowerTempBounds[i] = m_regionPts[i]->minTemp();
        if (m_regionPts[i]->speciesIndex() != m_index) {
            throw CanteraError("Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion",
                               "m_index inconsistency");
        }
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

Nasa9PolyMultiTempRegion::
Nasa9PolyMultiTempRegion(const Nasa9PolyMultiTempRegion& b) :
    SpeciesThermoInterpType(b),
    m_numTempRegions(b.m_numTempRegions),
    m_lowerTempBounds(b.m_lowerTempBounds),
    m_currRegion(b.m_currRegion)
{
    m_regionPts.resize(m_numTempRegions);
    for (size_t i = 0; i < m_numTempRegions; i++) {
        Nasa9Poly1* dptr = b.m_regionPts[i];
        m_regionPts[i] = new Nasa9Poly1(*dptr);
    }
}

Nasa9PolyMultiTempRegion&
Nasa9PolyMultiTempRegion::operator=(const Nasa9PolyMultiTempRegion& b)
{
    if (&b != this) {
        SpeciesThermoInterpType::operator=(b);
        for (size_t i = 0; i < m_numTempRegions; i++) {
            delete m_regionPts[i];
            m_regionPts[i] = 0;
        }
        m_numTempRegions = b.m_numTempRegions;
        m_lowerTempBounds = b.m_lowerTempBounds;
        m_currRegion = b.m_currRegion;
        m_regionPts.resize(m_numTempRegions);
        for (size_t i = 0; i < m_numTempRegions; i++) {
            m_regionPts[i] = new Nasa9Poly1(*(b.m_regionPts[i]));
        }
    }
    return *this;
}

Nasa9PolyMultiTempRegion::~Nasa9PolyMultiTempRegion()
{
    for (size_t i = 0; i < m_numTempRegions; i++) {
        delete m_regionPts[i];
        m_regionPts[i] = 0;
    }
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

void Nasa9PolyMultiTempRegion::updateProperties(const doublereal* tt,
        doublereal* cp_R,
        doublereal* h_RT,
        doublereal* s_R) const
{
    // Let's put some additional debugging here.
    // This is an external routine
#ifdef DEBUG_HKM
    double temp = tt[0];
    if (temp < m_regionPts[m_currRegion]->minTemp()) {
        if (m_currRegion != 0) {
            throw CanteraError("Nasa9PolyMultiTempRegion::updateProperties",
                               "region problem");
        }
    }
    if (temp > m_regionPts[m_currRegion]->maxTemp()) {
        if (m_currRegion != m_numTempRegions - 1) {
            throw CanteraError("Nasa9PolyMultiTempRegion::updateProperties",
                               "region problem");
        }
    }
#endif
    (m_regionPts[m_currRegion])->updateProperties(tt, cp_R, h_RT, s_R);
}

void Nasa9PolyMultiTempRegion::updatePropertiesTemp(const doublereal temp,
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
    // Now find the region
    m_currRegion = 0;
    for (size_t i = 1; i < m_numTempRegions; i++) {
        if (temp < m_lowerTempBounds[i]) {
            break;
        }
        m_currRegion++;
    }

    updateProperties(tPoly, cp_R, h_RT, s_R);
}

void Nasa9PolyMultiTempRegion::reportParameters(size_t& n, int& type,
        doublereal& tlow, doublereal& thigh,
        doublereal& pref,
        doublereal* const coeffs) const
{
    warn_deprecated("Nasa9PolyMultiTempRegion::reportParameters");
    n = m_index;
    type = NASA9MULTITEMP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    double ctmp[12];
    coeffs[0] = double(m_numTempRegions);
    int index = 1;
    size_t n_tmp = 0;
    int type_tmp = 0;
    double pref_tmp = 0.0;
    for (size_t iReg = 0; iReg < m_numTempRegions; iReg++) {
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
    warn_deprecated("Nasa9PolyMultiTempRegion::modifyParameters");
    int index = 3;
    for (size_t iReg = 0; iReg < m_numTempRegions; iReg++) {
        m_regionPts[iReg]->modifyParameters(coeffs + index);
        index += 11;
    }
}

}
