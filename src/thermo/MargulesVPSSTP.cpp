/**
 *  @file MargulesVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess Gibbs free energy formulations related to Margules
 *   expansions (see @ref thermoprops
 *    and class @link Cantera::MargulesVPSSTP MargulesVPSSTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{
MargulesVPSSTP::MargulesVPSSTP(const string& inputFile, const string& id_)
{
    initThermoFile(inputFile, id_);
}

// -- Activities, Standard States, Activity Concentrations -----------

void MargulesVPSSTP::getLnActivityCoefficients(double* lnac) const
{
    // Update the activity coefficients
    s_update_lnActCoeff();

    // take the exp of the internally stored coefficients.
    for (size_t k = 0; k < m_kk; k++) {
        lnac[k] = lnActCoeff_Scaled_[k];
    }
}

// ------------ Partial Molar Properties of the Solution ------------

void MargulesVPSSTP::getChemPotentials(double* mu) const
{
    // First get the standard chemical potentials in molar form. This requires
    // updates of standard state as a function of T and P
    getStandardChemPotentials(mu);

    // Update the activity coefficients
    s_update_lnActCoeff();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        mu[k] += RT() * (log(xx) + lnActCoeff_Scaled_[k]);
    }
}

double MargulesVPSSTP::cv_mole() const
{
    return cp_mole() - GasConstant;
}

void MargulesVPSSTP::getPartialMolarEnthalpies(double* hbar) const
{
    // Get the nondimensional standard state enthalpies
    getEnthalpy_RT(hbar);

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RT() * temperature() * dlnActCoeffdT_Scaled_[k];
    }
}

void MargulesVPSSTP::getPartialMolarCp(double* cpbar) const
{
    // Get the nondimensional standard state entropies
    getCp_R(cpbar);
    double T = temperature();

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] -= 2 * T * dlnActCoeffdT_Scaled_[k] + T * T * d2lnActCoeffdT2_Scaled_[k];
    }
    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void MargulesVPSSTP::getPartialMolarEntropies(double* sbar) const
{
    // Get the nondimensional standard state entropies
    getEntropy_R(sbar);
    double T = temperature();

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void MargulesVPSSTP::getPartialMolarVolumes(double* vbar) const
{
    double T = temperature();

    // Get the standard state values in m^3 kmol-1
    getStandardVolumes(vbar);

    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        double g0 = (m_VHE_b_ij[i] - T * m_VSE_b_ij[i]);
        double g1 = (m_VHE_c_ij[i] - T * m_VSE_c_ij[i]);
        const double temp1 = g0 + g1 * XB;
        const double all = -1.0*XA*XB*temp1 - XA*XB*XB*g1;

        for (size_t iK = 0; iK < m_kk; iK++) {
            vbar[iK] += all;
        }
        vbar[iA] += XB * temp1;
        vbar[iB] += XA * temp1 + XA*XB*g1;
    }
}

void MargulesVPSSTP::initThermo()
{
    initLengths();
    if (m_input.hasKey("interactions")) {
        for (auto& item : m_input["interactions"].asVector<AnyMap>()) {
            auto& species = item["species"].asVector<string>(2);
            vector<double> h(2), s(2), vh(2), vs(2);
            if (item.hasKey("excess-enthalpy")) {
                h = item.convertVector("excess-enthalpy", "J/kmol", 2);
            }
            if (item.hasKey("excess-entropy")) {
                s = item.convertVector("excess-entropy", "J/kmol/K", 2);
            }
            if (item.hasKey("excess-volume-enthalpy")) {
                vh = item.convertVector("excess-volume-enthalpy", "m^3/kmol", 2);
            }
            if (item.hasKey("excess-volume-entropy")) {
                vs = item.convertVector("excess-volume-entropy", "m^3/kmol/K", 2);
            }
            addBinaryInteraction(species[0], species[1],
                h[0], h[1], s[0], s[1], vh[0], vh[1], vs[0], vs[1]);
        }
    }
    GibbsExcessVPSSTP::initThermo();
}

void MargulesVPSSTP::getParameters(AnyMap& phaseNode) const
{
    GibbsExcessVPSSTP::getParameters(phaseNode);
    vector<AnyMap> interactions;
    for (size_t n = 0; n < m_pSpecies_A_ij.size(); n++) {
        AnyMap interaction;
        interaction["species"] = vector<string>{
            speciesName(m_pSpecies_A_ij[n]), speciesName(m_pSpecies_B_ij[n])};
        if (m_HE_b_ij[n] != 0 || m_HE_c_ij[n] != 0) {
            interaction["excess-enthalpy"].setQuantity(
                {m_HE_b_ij[n], m_HE_c_ij[n]}, "J/kmol");
        }
        if (m_SE_b_ij[n] != 0 || m_SE_c_ij[n] != 0) {
            interaction["excess-entropy"].setQuantity(
                {m_SE_b_ij[n], m_SE_c_ij[n]}, "J/kmol/K");
        }
        if (m_VHE_b_ij[n] != 0 || m_VHE_c_ij[n] != 0) {
            interaction["excess-volume-enthalpy"].setQuantity(
                {m_VHE_b_ij[n], m_VHE_c_ij[n]}, "m^3/kmol");
        }
        if (m_VSE_b_ij[n] != 0 || m_VSE_c_ij[n] != 0) {
            interaction["excess-volume-entropy"].setQuantity(
                {m_VSE_b_ij[n], m_VSE_c_ij[n]}, "m^3/kmol/K");
        }
        interactions.push_back(std::move(interaction));
    }
    phaseNode["interactions"] = std::move(interactions);
}

void MargulesVPSSTP::initLengths()
{
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
}

void MargulesVPSSTP::addBinaryInteraction(const string& speciesA,
    const string& speciesB, double h0, double h1, double s0, double s1,
    double vh0, double vh1, double vs0, double vs1)
{
    size_t kA = speciesIndex(speciesA);
    size_t kB = speciesIndex(speciesB);
    // The interaction is silently ignored if either species is not defined in
    // the current phase.
    if (kA == npos || kB == npos) {
        return;
    }
    m_pSpecies_A_ij.push_back(kA);
    m_pSpecies_B_ij.push_back(kB);

    m_HE_b_ij.push_back(h0);
    m_HE_c_ij.push_back(h1);
    m_SE_b_ij.push_back(s0);
    m_SE_c_ij.push_back(s1);
    m_VHE_b_ij.push_back(vh0);
    m_VHE_c_ij.push_back(vh1);
    m_VSE_b_ij.push_back(vs0);
    m_VSE_c_ij.push_back(vs1);
    numBinaryInteractions_++;
}


void MargulesVPSSTP::s_update_lnActCoeff() const
{
    double T = temperature();
    lnActCoeff_Scaled_.assign(m_kk, 0.0);
    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT();
        double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT();
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        const double XAXB = XA * XB;
        const double g0g1XB = (g0 + g1 * XB);
        const double all = -1.0 * XAXB * g0g1XB - XAXB * XB * g1;
        for (size_t iK = 0; iK < m_kk; iK++) {
            lnActCoeff_Scaled_[iK] += all;
        }
        lnActCoeff_Scaled_[iA] += XB * g0g1XB;
        lnActCoeff_Scaled_[iB] += XA * g0g1XB + XAXB * g1;
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dT() const
{
    double invT = 1.0 / temperature();
    double invRTT = 1.0 / GasConstant*invT*invT;
    dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
    d2lnActCoeffdT2_Scaled_.assign(m_kk, 0.0);
    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        double g0 = -m_HE_b_ij[i] * invRTT;
        double g1 = -m_HE_c_ij[i] * invRTT;
        const double XAXB = XA * XB;
        const double g0g1XB = (g0 + g1 * XB);
        const double all = -1.0 * XAXB * g0g1XB - XAXB * XB * g1;
        const double mult = 2.0 * invT;
        const double dT2all = mult * all;
        for (size_t iK = 0; iK < m_kk; iK++) {
            dlnActCoeffdT_Scaled_[iK] += all;
            d2lnActCoeffdT2_Scaled_[iK] -= dT2all;
        }
        dlnActCoeffdT_Scaled_[iA] += XB * g0g1XB;
        dlnActCoeffdT_Scaled_[iB] += XA * g0g1XB + XAXB * g1;
        d2lnActCoeffdT2_Scaled_[iA] -= mult * XB * g0g1XB;
        d2lnActCoeffdT2_Scaled_[iB] -= mult * (XA * g0g1XB + XAXB * g1);
    }
}

void MargulesVPSSTP::getdlnActCoeffdT(double* dlnActCoeffdT) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
}

void MargulesVPSSTP::getd2lnActCoeffdT2(double* d2lnActCoeffdT2) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        d2lnActCoeffdT2[k] = d2lnActCoeffdT2_Scaled_[k];
    }
}

void MargulesVPSSTP::getdlnActCoeffds(const double dTds, const double* const dXds,
                                       double* dlnActCoeffds) const
{
    double T = temperature();
    s_update_dlnActCoeff_dT();
    for (size_t iK = 0; iK < m_kk; iK++) {
        dlnActCoeffds[iK] = 0.0;
    }

    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        double dXA = dXds[iA];
        double dXB = dXds[iB];
        double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT();
        double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT();
        const double g02g1XB = g0 + 2*g1*XB;
        const double g2XAdXB = 2*g1*XA*dXB;
        const double all = (-XB * dXA - XA *dXB) * g02g1XB - XB *g2XAdXB;
        for (size_t iK = 0; iK < m_kk; iK++) {
            dlnActCoeffds[iK] += all + dlnActCoeffdT_Scaled_[iK]*dTds;
        }
        dlnActCoeffds[iA] += dXB * g02g1XB;
        dlnActCoeffds[iB] += dXA * g02g1XB + g2XAdXB;
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dlnN_diag() const
{
    double T = temperature();
    dlnActCoeffdlnN_diag_.assign(m_kk, 0.0);

    for (size_t iK = 0; iK < m_kk; iK++) {
        double XK = moleFractions_[iK];

        for (size_t i = 0; i < numBinaryInteractions_; i++) {
            size_t iA = m_pSpecies_A_ij[i];
            size_t iB = m_pSpecies_B_ij[i];
            size_t delAK = 0;
            size_t delBK = 0;

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            double XA = moleFractions_[iA];
            double XB = moleFractions_[iB];

            double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT();
            double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT();

            dlnActCoeffdlnN_diag_[iK] += 2*(delBK-XB)*(g0*(delAK-XA)+g1*(2*(delAK-XA)*XB+XA*(delBK-XB)));
        }
        dlnActCoeffdlnN_diag_[iK] = XK*dlnActCoeffdlnN_diag_[iK];
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dlnN() const
{
    double T = temperature();
    dlnActCoeffdlnN_.zero();

    // Loop over the activity coefficient gamma_k
    for (size_t iK = 0; iK < m_kk; iK++) {
        for (size_t iM = 0; iM < m_kk; iM++) {
            double XM = moleFractions_[iM];
            for (size_t i = 0; i < numBinaryInteractions_; i++) {
                size_t iA = m_pSpecies_A_ij[i];
                size_t iB = m_pSpecies_B_ij[i];
                double delAK = 0.0;
                double delBK = 0.0;
                double delAM = 0.0;
                double delBM = 0.0;
                if (iA==iK) {
                    delAK = 1.0;
                } else if (iB==iK) {
                    delBK = 1.0;
                }
                if (iA==iM) {
                    delAM = 1.0;
                } else if (iB==iM) {
                    delBM = 1.0;
                }

                double XA = moleFractions_[iA];
                double XB = moleFractions_[iB];
                double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT();
                double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT();
                dlnActCoeffdlnN_(iK,iM) += g0*((delAM-XA)*(delBK-XB)+(delAK-XA)*(delBM-XB));
                dlnActCoeffdlnN_(iK,iM) += 2*g1*((delAM-XA)*(delBK-XB)*XB+(delAK-XA)*(delBM-XB)*XB+(delBM-XB)*(delBK-XB)*XA);
            }
            dlnActCoeffdlnN_(iK,iM) = XM*dlnActCoeffdlnN_(iK,iM);
        }
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dlnX_diag() const
{
    double T = temperature();
    dlnActCoeffdlnX_diag_.assign(m_kk, 0.0);

    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];

        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];

        double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT();
        double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT();

        dlnActCoeffdlnX_diag_[iA] += XA*XB*(2*g1*-2*g0-6*g1*XB);
        dlnActCoeffdlnX_diag_[iB] += XA*XB*(2*g1*-2*g0-6*g1*XB);
    }
}

void MargulesVPSSTP::getdlnActCoeffdlnN_diag(double* dlnActCoeffdlnN_diag) const
{
    s_update_dlnActCoeff_dlnN_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnN_diag[k] = dlnActCoeffdlnN_diag_[k];
    }
}

void MargulesVPSSTP::getdlnActCoeffdlnX_diag(double* dlnActCoeffdlnX_diag) const
{
    s_update_dlnActCoeff_dlnX_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}

void MargulesVPSSTP::getdlnActCoeffdlnN(const size_t ld, double* dlnActCoeffdlnN)
{
    s_update_dlnActCoeff_dlnN();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}

}
