/**
 *  @file Nasa9PolyMultiTempRegion.cpp
 *  Definitions for a single-species standard state object derived
 *  from @link Cantera::SpeciesThermoInterpType
 *    SpeciesThermoInterpType@endlink  based
 *  on the NASA 9 coefficient temperature polynomial form
 *  applied to one temperature region
 *  (see @ref spthermo and class
 *   @link Cantera::Nasa9Poly1 Nasa9Poly1@endlink).
 *
 *  This parameterization has one NASA temperature region.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/Nasa9PolyMultiTempRegion.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion(span<Nasa9Poly1*> regionPts)
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

Nasa9PolyMultiTempRegion::Nasa9PolyMultiTempRegion(double tlow, double thigh, double pref,
                                                   span<const double> coeffs)
    : SpeciesThermoInterpType(tlow, thigh, pref)
{
    size_t regions = static_cast<size_t>(coeffs[0]);

    for (size_t i=0; i<regions; i++) {
        Nasa9Poly1* poly = new Nasa9Poly1(coeffs[11*i+1], coeffs[11*i+2],
                                          pref, coeffs.subspan(11*i + 3, 9));
        m_regionPts.emplace_back(poly);
    }

    m_lowerTempBounds.resize(regions);
    for (size_t i = 0; i < m_regionPts.size(); i++) {
        m_lowerTempBounds[i] = m_regionPts[i]->minTemp();
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

void Nasa9PolyMultiTempRegion::setParameters(const map<double, vector<double>>& regions)
{
    m_regionPts.clear();
    m_lowerTempBounds.clear();
    for (const auto& [Tmin, coeffs] : regions) {
        m_lowerTempBounds.push_back(Tmin);
        Nasa9Poly1* poly = new Nasa9Poly1;
        poly->setRefPressure(refPressure());
        poly->setMinTemp(Tmin);
        poly->setParameters(coeffs);
        if (!m_regionPts.empty()) {
            m_regionPts.back()->setMaxTemp(Tmin);
        }
        m_regionPts.emplace_back(poly);
    }
    m_regionPts.back()->setMaxTemp(maxTemp());
}

int Nasa9PolyMultiTempRegion::reportType() const
{
    return NASA9MULTITEMP;
}

void Nasa9PolyMultiTempRegion::updateTemperaturePoly(double T, span<double> T_poly) const
{
    T_poly[0] = T;
    T_poly[1] = T * T;
    T_poly[2] = T_poly[1] * T;
    T_poly[3] = T_poly[2] * T;
    T_poly[4] = 1.0 / T;
    T_poly[5] = T_poly[4] / T;
    T_poly[6] = std::log(T);
}

void Nasa9PolyMultiTempRegion::updateProperties(span<const double> tt,
        double& cp_R, double& h_RT, double& s_R) const
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

void Nasa9PolyMultiTempRegion::updatePropertiesTemp(const double temp,
        double& cp_R, double& h_RT, double& s_R) const
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

size_t Nasa9PolyMultiTempRegion::nCoeffs() const
{
    return 11*m_regionPts.size() + 1;
}

void Nasa9PolyMultiTempRegion::reportParameters(size_t& n, int& type,
        double& tlow, double& thigh, double& pref, span<double> coeffs) const
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
        m_regionPts[iReg]->reportParameters(
            n_tmp, type_tmp, coeffs[index], coeffs[index+1], pref_tmp, ctmp);
        for (int i = 0; i < 9; i++) {
            coeffs[index+2+i] = ctmp[3+i];
        }
        index += 11;
    }
}

void Nasa9PolyMultiTempRegion::getParameters(AnyMap& thermo) const
{
    thermo["model"] = "NASA9";
    SpeciesThermoInterpType::getParameters(thermo);
    auto T_ranges = m_lowerTempBounds;
    T_ranges.push_back(m_highT);
    thermo["temperature-ranges"].setQuantity(T_ranges, "K");
    thermo["data"] = vector<vector<double>>();
    for (const auto& region : m_regionPts) {
        region->getParameters(thermo);
    }
}

}
