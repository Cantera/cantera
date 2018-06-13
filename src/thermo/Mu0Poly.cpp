/**
 *  @file Mu0Poly.cpp
 *  Definitions for a single-species standard state object derived
 *  from \link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType\endlink  based
 *  on a piecewise constant mu0 interpolation
 *  (see \ref spthermo and class \link Cantera::Mu0Poly Mu0Poly\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Mu0Poly.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

Mu0Poly::Mu0Poly(double tlow, double thigh, double pref, const double* coeffs) :
    SpeciesThermoInterpType(tlow, thigh, pref),
    m_numIntervals(0),
    m_H298(0.0)
{
    processCoeffs(coeffs);
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
        throw CanteraError("installMu0ThermoFromXML", "missing Mu0Values");
    }
    getFloatArray(*valNode_ptr, cValues, true, "actEnergy");

    // Check to see whether the Mu0's were input in a dimensionless form. If
    // they were, then the assumed temperature needs to be adjusted from the
    // assumed T = 273.15
    if (valNode_ptr->attrib("units") == "Dimensionless") {
        dimensionlessMu0Values = true;
    }
    if (cValues.size() != numPoints) {
        throw CanteraError("installMu0ThermoFromXML", "numPoints inconsistent");
    }

    vector_fp cTemperatures(numPoints);
    const XML_Node* tempNode_ptr = getByTitle(Mu0Node, "Mu0Temperatures");
    if (!tempNode_ptr) {
        throw CanteraError("installMu0ThermoFromXML",
                           "missing Mu0Temperatures");
    }
    getFloatArray(*tempNode_ptr, cTemperatures, false);
    if (cTemperatures.size() != numPoints) {
        throw CanteraError("installMu0ThermoFromXML", "numPoints inconsistent");
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

void Mu0Poly::processCoeffs(const doublereal* coeffs)
{
    size_t nPoints = (size_t) coeffs[0];
    if (nPoints < 2) {
        throw CanteraError("Mu0Poly",
                           "nPoints must be >= 2");
    }
    m_numIntervals = nPoints - 1;
    m_H298 = coeffs[1] / GasConstant;
    size_t iT298 = 0;

    // Resize according to the number of points
    m_t0_int.resize(nPoints);
    m_h0_R_int.resize(nPoints);
    m_s0_R_int.resize(nPoints);
    m_cp0_R_int.resize(nPoints);
    m_mu0_R_int.resize(nPoints);

    // Calculate the T298 interval and make sure that the temperatures are
    // strictly monotonic. Also distribute the data into the internal arrays.
    bool ifound = false;
    for (size_t i = 0, iindex = 2; i < nPoints; i++) {
        double T1 = coeffs[iindex];
        m_t0_int[i] = T1;
        m_mu0_R_int[i] = coeffs[iindex+1] / GasConstant;
        if (T1 == 298.15) {
            iT298 = i;
            ifound = true;
        }
        if (i < nPoints - 1 && coeffs[iindex+2] <= T1) {
            throw CanteraError("Mu0Poly",
                               "Temperatures are not monotonic increasing");
        }
        iindex += 2;
    }
    if (!ifound) {
        throw CanteraError("Mu0Poly",
                           "One temperature has to be 298.15");
    }

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

}
