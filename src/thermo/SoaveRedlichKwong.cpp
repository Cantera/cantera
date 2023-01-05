//! @file SoaveRedlichKwong.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SoaveRedlichKwong.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"

#include <boost/math/tools/roots.hpp>

namespace bmt = boost::math::tools;

namespace Cantera
{
const double SoaveRedlichKwong::omega_a = 4.27480233540E-01;
const double SoaveRedlichKwong::omega_b = 8.66403499650E-02;
const double SoaveRedlichKwong::omega_vc = 3.33333333333333E-01;



double SoaveRedlichKwong::dpdVCalc(double T, double molarVol, double& presCalc) const
{
    double denom = molarVol * (molarVol + m_b);
    double vmb = molarVol - m_b;
    return -GasConstant * T / (vmb * vmb) + m_aAlpha_mix * (2 * molarVol + m_b) / (denom * denom);
}

void SoaveRedlichKwong::calculatePressureDerivatives() const
{
    double T = temperature();
    double mv = molarVolume();
    double pres;

    m_dpdV = dpdVCalc(T, mv, pres);
    m_dpdT = GasConstant / (mv - m_b) - daAlpha_dT() / (mv * (mv + m_b));
}

double SoaveRedlichKwong::daAlpha_dT() const
{
    // Same as PengRobinson::daAlpha_dT()
    double daAlphadT = 0.0, k, Tc, sqtTr, coeff1, coeff2;
    for (size_t i = 0; i < m_kk; i++) {
        // Calculate first derivative of alpha for individual species
        Tc = speciesCritTemperature(m_a_coeffs(i,i), m_b_coeffs[i]);
        sqtTr = sqrt(temperature() / Tc);
        coeff1 = 1 / (Tc*sqtTr);
        coeff2 = sqtTr - 1;
        k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
    }
    // Calculate mixture derivative
    for (size_t i = 0; i < m_kk; i++) {
        for (size_t j = 0; j < m_kk; j++) {
            daAlphadT += moleFractions_[i] * moleFractions_[j] * 0.5
                         * m_aAlpha_binary(i, j)
                         * (m_dalphadT[i] / m_alpha[i] + m_dalphadT[j] / m_alpha[j]);
        }
    }
    return daAlphadT;
}

double SoaveRedlichKwong::d2aAlpha_dT2() const
{
    // Same as PengRobinson::d2aAlpha_dT2()
    for (size_t i = 0; i < m_kk; i++) {
        double Tcrit_i = speciesCritTemperature(m_a_coeffs(i, i), m_b_coeffs[i]);
        double sqt_Tr = sqrt(temperature() / Tcrit_i);
        double coeff1 = 1 / (Tcrit_i*sqt_Tr);
        double coeff2 = sqt_Tr - 1;
        // Calculate first and second derivatives of alpha for individual species
        double k = m_kappa[i];
        m_dalphadT[i] = coeff1 * (k*k*coeff2 - k);
        m_d2alphadT2[i] = (k*k + k) * coeff1 / (2*sqt_Tr*sqt_Tr*Tcrit_i);
    }

    // Calculate mixture derivative
    double d2aAlphadT2 = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        double alphai = m_alpha[i];
        for (size_t j = 0; j < m_kk; j++) {
            double alphaj = m_alpha[j];
            double alphaij = alphai * alphaj;
            double term1 = m_d2alphadT2[i] / alphai + m_d2alphadT2[j] / alphaj;
            double term2 = 2 * m_dalphadT[i] * m_dalphadT[j] / alphaij;
            double term3 = m_dalphadT[i] / alphai + m_dalphadT[j] / alphaj;
            d2aAlphadT2 += 0.5 * moleFractions_[i] * moleFractions_[j]
                           * m_aAlpha_binary(i, j)
                           * (term1 + term2 - 0.5 * term3 * term3);
        }
    }
    return d2aAlphadT2;
}


int SoaveRedlichKwong::solveCubic(double T, double pres, double a, double b, double aAlpha,
                                  double Vroot[3]) const
{
    double an = 1.0;
    double bn = - GasConstant * T / pres;
    double cn = (aAlpha - b * GasConstant * T) / pres - b * b;
    double dn = aAlpha * b / pres;

    double tc = a * omega_b / (b * omega_a * GasConstant);
    double pc = omega_b * R * tc / b;
    double vc = omega_vc * GasConstant * tc / pc;

    return MixtureFugacityTP::solveCubic(T, pres, a, b, aAlpha, Vroot,
                                         an, bn, cn, dn, tc, vc);
}

}
