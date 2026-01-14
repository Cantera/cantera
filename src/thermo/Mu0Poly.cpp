/**
 *  @file Mu0Poly.cpp
 *  Definitions for a single-species standard state object derived
 *  from @link Cantera::SpeciesThermoInterpType SpeciesThermoInterpType@endlink  based
 *  on a piecewise constant mu0 interpolation
 *  (see @ref spthermo and class @link Cantera::Mu0Poly Mu0Poly@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Mu0Poly.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{
Mu0Poly::Mu0Poly()
    : SpeciesThermoInterpType(0.0, std::numeric_limits<double>::infinity(), 0.0)
{
}

Mu0Poly::Mu0Poly(double tlow, double thigh, double pref, span<const double> coeffs) :
    SpeciesThermoInterpType(tlow, thigh, pref),
    m_numIntervals(0),
    m_H298(0.0)
{
    map<double, double> T_mu;
    size_t nPoints = (size_t) coeffs[0];
    for (size_t i = 0; i < nPoints; i++) {
        T_mu[coeffs[2*i+2]] = coeffs[2*i+3];
    }
    setParameters(coeffs[1], T_mu);
}

void Mu0Poly::setParameters(double h0, const map<double, double>& T_mu)
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
    for (const auto& [T1, mu] : T_mu) {
        if (T1 == 298.15) {
            iT298 = m_t0_int.size();
        }
        m_t0_int.push_back(T1);
        m_mu0_R_int.push_back(mu / GasConstant);
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

void Mu0Poly::updateProperties(span<const double> tt, double& cp_R,
                               double& h_RT, double& s_R) const
{
    size_t j = m_numIntervals;
    double T = tt[0];
    for (size_t i = 0; i < m_numIntervals; i++) {
        double T2 = m_t0_int[i+1];
        if (T <=T2) {
            j = i;
            break;
        }
    }
    double T1 = m_t0_int[j];
    double cp_Rj = m_cp0_R_int[j];
    cp_R = cp_Rj;
    h_RT = (m_h0_R_int[j] + (T - T1) * cp_Rj)/T;
    s_R = m_s0_R_int[j] + cp_Rj * (log(T/T1));
}

void Mu0Poly::updatePropertiesTemp(const double T, double& cp_R,
                                   double& h_RT, double& s_R) const
{
    updateProperties(span<const double>(&T, 1), cp_R, h_RT, s_R);
}

size_t Mu0Poly::nCoeffs() const
{
  return 2*m_numIntervals + 4;
}

void Mu0Poly::reportParameters(size_t& n, int& type, double& tlow, double& thigh,
                               double& pref, span<double> coeffs) const
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

void Mu0Poly::getParameters(AnyMap& thermo) const
{
    SpeciesThermoInterpType::getParameters(thermo);
    thermo["model"] = "piecewise-Gibbs";
    thermo["h0"].setQuantity(m_H298 * GasConstant, "J/kmol");
    AnyMap data;
    bool dimensionless = m_input.getBool("dimensionless", false);
    if (dimensionless) {
        thermo["dimensionless"] = true;
    }
    for (size_t i = 0; i < m_numIntervals+1; i++) {
        if (dimensionless) {
            data[fmt::format("{}", m_t0_int[i])] = m_mu0_R_int[i] / m_t0_int[i];
        } else {
            data[fmt::format("{}", m_t0_int[i])].setQuantity(
                m_mu0_R_int[i] * GasConstant, "J/kmol");
        }
    }
    thermo["data"] = std::move(data);
}

}
