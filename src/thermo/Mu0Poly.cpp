/**
 *  @file Mu0Poly.cpp
 *  Definitions for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on a piecewise constant mu0 interpolation
 *  (see \ref spthermo and class \link Cantera::Mu0Poly Mu0Poly\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Mu0Poly.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
Mu0Poly::Mu0Poly()
    : m_numIntervals(0)
    , m_H298(0.0)
{
}

Mu0Poly::Mu0Poly(double tlow, double thigh, double pref, const double* coeffs) :
    SpeciesThermoInterpType(tlow, thigh, pref),
    m_numIntervals(0),
    m_H298(0.0)
{
    std::map<double, double> T_mu;
    size_t nPoints = (size_t) coeffs[0];
    for (size_t i = 0; i < nPoints; i++) {
        T_mu[coeffs[2*i+2]] = coeffs[2*i+3];
    }
    setParameters(coeffs[1], T_mu);
}

void Mu0Poly::setParameters(double h0, const std::map<double, double>& T_mu)
{
    size_t nPoints = T_mu.size();
    if (nPoints < 2) {
        throw CanteraError("Mu0Poly::setParameters", "nPoints must be >= 2");
    }
    m_numIntervals = nPoints - 1;
    m_H298 = h0 / GasConstant;

    // Distribute the data into the internal arrays, and find the index of the
    // point at 298.15 K.
    size_t iT298 = npos;
    for (const auto& row : T_mu) {
        double T1 = row.first;
        if (T1 == 298.15) {
            iT298 = m_t0_int.size();
        }
        m_t0_int.push_back(T1);
        m_mu0_R_int.push_back(row.second / GasConstant);
    }
    if (iT298 == npos) {
        throw CanteraError("Mu0Poly::setParameters",
                           "One temperature has to be 298.15");
    }

    // Resize according to the number of points
    m_h0_R_int.resize(nPoints);
    m_s0_R_int.resize(nPoints);
    m_cp0_R_int.resize(nPoints);

    // Starting from the interval with T298, we go up
    m_h0_R_int[iT298] = m_H298;
    m_s0_R_int[iT298] = - (m_mu0_R_int[iT298] - m_h0_R_int[iT298]) / m_t0_int[iT298];
    for (size_t i = iT298; i < m_numIntervals; i++) {
        double T1 = m_t0_int[i];
        double s1 = m_s0_R_int[i];
        double T2 = m_t0_int[i+1];
        double deltaMu = m_mu0_R_int[i+1] - m_mu0_R_int[i];
        double deltaT = T2 - T1;
        double cpi = (deltaMu - T1 * s1 + T2 * s1) / (deltaT - T2 * log(T2/T1));
        m_cp0_R_int[i] = cpi;
        m_h0_R_int[i+1] = m_h0_R_int[i] + cpi * deltaT;
        m_s0_R_int[i+1] = s1 + cpi * log(T2/T1);
        m_cp0_R_int[i+1] = cpi;
    }

    // Starting from the interval with T298, we go down
    if (iT298 != 0) {
        m_h0_R_int[iT298] = m_H298;
        m_s0_R_int[iT298] = - (m_mu0_R_int[iT298] - m_h0_R_int[iT298]) / m_t0_int[iT298];
        for (size_t i = iT298 - 1; i != npos; i--) {
            double T1 = m_t0_int[i];
            double T2 = m_t0_int[i+1];
            double s2 = m_s0_R_int[i+1];
            double deltaMu = m_mu0_R_int[i+1] - m_mu0_R_int[i];
            double deltaT = T2 - T1;
            double cpi = (deltaMu - T1 * s2 + T2 * s2) / (deltaT - T1 * log(T2/T1));
            m_cp0_R_int[i] = cpi;
            m_h0_R_int[i] = m_h0_R_int[i+1] - cpi * deltaT;
            m_s0_R_int[i] = s2 - cpi * log(T2/T1);
            if (i == (m_numIntervals-1)) {
                m_cp0_R_int[i+1] = cpi;
            }
        }
    }
}

void Mu0Poly::updateProperties(const doublereal* tt, doublereal* cp_R,
                               doublereal* h_RT, doublereal* s_R) const
{
    size_t j = m_numIntervals;
    double T = *tt;
    for (size_t i = 0; i < m_numIntervals; i++) {
        double T2 = m_t0_int[i+1];
        if (T <=T2) {
            j = i;
            break;
        }
    }
    double T1 = m_t0_int[j];
    double cp_Rj = m_cp0_R_int[j];
    *cp_R = cp_Rj;
    *h_RT = (m_h0_R_int[j] + (T - T1) * cp_Rj)/T;
    *s_R = m_s0_R_int[j] + cp_Rj * (log(T/T1));
}

void Mu0Poly::updatePropertiesTemp(const doublereal T,
                                   doublereal* cp_R,
                                   doublereal* h_RT,
                                   doublereal* s_R) const
{
    updateProperties(&T, cp_R, h_RT, s_R);
}

size_t Mu0Poly::nCoeffs() const
{
  return 2*m_numIntervals + 4;
}

void Mu0Poly::reportParameters(size_t& n, int& type,
                               doublereal& tlow, doublereal& thigh,
                               doublereal& pref,
                               doublereal* const coeffs) const
{
    n = 0;
    type = MU0_INTERP;
    tlow = m_lowT;
    thigh = m_highT;
    pref = m_Pref;
    coeffs[0] = int(m_numIntervals)+1;
    coeffs[1] = m_H298 * GasConstant;
    int j = 2;
    for (size_t i = 0; i < m_numIntervals+1; i++) {
        coeffs[j] = m_t0_int[i];
        coeffs[j+1] = m_mu0_R_int[i] * GasConstant;
        j += 2;
    }
}

Mu0Poly* newMu0ThermoFromXML(const XML_Node& Mu0Node)
{
    bool dimensionlessMu0Values = false;

    doublereal h298 = 0.0;
    if (Mu0Node.hasChild("H298")) {
        h298 = getFloat(Mu0Node, "H298", "actEnergy");
    }

    size_t numPoints = 1;
    if (Mu0Node.hasChild("numPoints")) {
        numPoints = getInteger(Mu0Node, "numPoints");
    }

    vector_fp cValues(numPoints);
    const XML_Node* valNode_ptr = getByTitle(Mu0Node, "Mu0Values");
    if (!valNode_ptr) {
        throw CanteraError("newMu0ThermoFromXML", "missing Mu0Values");
    }
    getFloatArray(*valNode_ptr, cValues, true, "actEnergy");

    // Check to see whether the Mu0's were input in a dimensionless form. If
    // they were, then the assumed temperature needs to be adjusted from the
    // assumed T = 273.15
    if (valNode_ptr->attrib("units") == "Dimensionless") {
        dimensionlessMu0Values = true;
    }
    if (cValues.size() != numPoints) {
        throw CanteraError("newMu0ThermoFromXML", "numPoints inconsistent");
    }

    vector_fp cTemperatures(numPoints);
    const XML_Node* tempNode_ptr = getByTitle(Mu0Node, "Mu0Temperatures");
    if (!tempNode_ptr) {
        throw CanteraError("newMu0ThermoFromXML", "missing Mu0Temperatures");
    }
    getFloatArray(*tempNode_ptr, cTemperatures, false);
    if (cTemperatures.size() != numPoints) {
        throw CanteraError("newMu0ThermoFromXML", "numPoints inconsistent");
    }

    // Fix up dimensionless Mu0 values if input
    if (dimensionlessMu0Values) {
        for (size_t i = 0; i < numPoints; i++) {
            cValues[i] *= cTemperatures[i] / 273.15;
        }
    }

    vector_fp c(2 + 2 * numPoints);
    c[0] = static_cast<double>(numPoints);
    c[1] = h298;
    for (size_t i = 0; i < numPoints; i++) {
        c[2+i*2] = cTemperatures[i];
        c[2+i*2+1] = cValues[i];
    }

    return new Mu0Poly(fpValue(Mu0Node["Tmin"]), fpValue(Mu0Node["Tmax"]),
                       fpValue(Mu0Node["Pref"]), &c[0]);
}

}
