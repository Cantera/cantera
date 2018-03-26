/**
 *  @file PhaseCombo_Interaction.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PhaseCombo_Interaction.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{
PhaseCombo_Interaction::PhaseCombo_Interaction() :
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    warn_deprecated("Class PhaseCombo_Interaction", "To be removed after Cantera 2.4");
}

PhaseCombo_Interaction::PhaseCombo_Interaction(const std::string& inputFile,
        const std::string& id_) :
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    warn_deprecated("Class PhaseCombo_Interaction", "To be removed after Cantera 2.4");
    initThermoFile(inputFile, id_);
}

PhaseCombo_Interaction::PhaseCombo_Interaction(XML_Node& phaseRoot,
        const std::string& id_) :
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    warn_deprecated("Class PhaseCombo_Interaction", "To be removed after Cantera 2.4");
    importPhase(phaseRoot, this);
}

// - Activities, Standard States, Activity Concentrations -----------

void PhaseCombo_Interaction::getActivityCoefficients(doublereal* ac) const
{
    // Update the activity coefficients
    s_update_lnActCoeff();

    // take the exp of the internally stored coefficients.
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = exp(lnActCoeff_Scaled_[k]);
    }
}

// ------------ Partial Molar Properties of the Solution ------------

void PhaseCombo_Interaction::getChemPotentials(doublereal* mu) const
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

doublereal PhaseCombo_Interaction::enthalpy_mole() const
{
    double h = 0;
    vector_fp hbar(m_kk);
    getPartialMolarEnthalpies(&hbar[0]);
    for (size_t i = 0; i < m_kk; i++) {
        h += moleFractions_[i]*hbar[i];
    }
    return h;
}

doublereal PhaseCombo_Interaction::entropy_mole() const
{
    double s = 0;
    vector_fp sbar(m_kk);
    getPartialMolarEntropies(&sbar[0]);
    for (size_t i = 0; i < m_kk; i++) {
        s += moleFractions_[i]*sbar[i];
    }
    return s;
}

doublereal PhaseCombo_Interaction::cp_mole() const
{
    double cp = 0;
    vector_fp cpbar(m_kk);
    getPartialMolarCp(&cpbar[0]);
    for (size_t i = 0; i < m_kk; i++) {
        cp += moleFractions_[i]*cpbar[i];
    }
    return cp;
}

doublereal PhaseCombo_Interaction::cv_mole() const
{
    return cp_mole() - GasConstant;
}

void PhaseCombo_Interaction::getPartialMolarEnthalpies(doublereal* hbar) const
{
    // Get the nondimensional standard state enthalpies
    getEnthalpy_RT(hbar);
    // dimensionalize it.
    double T = temperature();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= GasConstant * T;
    }

    // Update the activity coefficients, This also update the internally stored
    // molalities.
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= GasConstant * T * T * dlnActCoeffdT_Scaled_[k];
    }
}

void PhaseCombo_Interaction::getPartialMolarCp(doublereal* cpbar) const
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

void PhaseCombo_Interaction::getPartialMolarEntropies(doublereal* sbar) const
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
        sbar[k] += - lnActCoeff_Scaled_[k] - log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }

    // dimensionalize it.
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void PhaseCombo_Interaction::getPartialMolarVolumes(doublereal* vbar) const
{
    double T = temperature();

    // Get the standard state values in m^3 kmol-1
    getStandardVolumes(vbar);

    for (size_t iK = 0; iK < m_kk; iK++) {
        int delAK = 0;
        int delBK = 0;
        for (size_t i = 0; i < numBinaryInteractions_; i++) {
            size_t iA = m_pSpecies_A_ij[i];
            size_t iB = m_pSpecies_B_ij[i];

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            double XA = moleFractions_[iA];
            double XB = moleFractions_[iB];
            double g0 = (m_VHE_b_ij[i] - T * m_VSE_b_ij[i]);
            double g1 = (m_VHE_c_ij[i] - T * m_VSE_c_ij[i]);
            vbar[iK] += XA*XB*(g0+g1*XB)+((delAK-XA)*XB+XA*(delBK-XB))*(g0+g1*XB)+XA*XB*(delBK-XB)*g1;
        }
    }
}

void PhaseCombo_Interaction::initThermo()
{
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}

void PhaseCombo_Interaction::initLengths()
{
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
}

void PhaseCombo_Interaction::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    if ((int) id.size() > 0 && phaseNode.id() != id) {
        throw CanteraError("PhaseCombo_Interaction::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    // Check on the thermo field. Must have:
    // <thermo model="PhaseCombo_Interaction" />
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("PhaseCombo_Interaction::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    if (!caseInsensitiveEquals(thermoNode["model"], "phasecombo_interaction")) {
        throw CanteraError("PhaseCombo_Interaction::initThermoXML",
                           "model name isn't PhaseCombo_Interaction: " + thermoNode["model"]);
    }

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        if (!caseInsensitiveEquals(acNode["model"], "margules")) {
            throw CanteraError("PhaseCombo_Interaction::initThermoXML",
                               "Unknown activity coefficient model: " + acNode["model"]);
        }
        for (size_t i = 0; i < acNode.nChildren(); i++) {
            XML_Node& xmlACChild = acNode.child(i);

            // Process a binary salt field, or any of the other XML fields that
            // make up the Pitzer Database. Entries will be ignored if any of
            // the species in the entry isn't in the solution.
            if (caseInsensitiveEquals(xmlACChild.name(), "binaryneutralspeciesparameters")) {
                readXMLBinarySpecies(xmlACChild);
            }
        }
    }

    // Go down the chain
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id);
}

void PhaseCombo_Interaction::s_update_lnActCoeff() const
{
    doublereal T = temperature();
    lnActCoeff_Scaled_.assign(m_kk, 0.0);

    for (size_t iK = 0; iK < m_kk; iK++) {
        // We never sample the end of the mole fraction domains
        double xx = std::max(moleFractions_[iK], SmallNumber);

        // First wipe out the ideal solution mixing term
        lnActCoeff_Scaled_[iK] = - log(xx);

        // Then add in the Margules interaction terms. that's it!
        for (size_t i = 0; i < numBinaryInteractions_; i++) {
            size_t iA = m_pSpecies_A_ij[i];
            size_t iB = m_pSpecies_B_ij[i];
            int delAK = 0;
            int delBK = 0;
            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }
            double XA = moleFractions_[iA];
            double XB = moleFractions_[iB];
            double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / (GasConstant*T);
            double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / (GasConstant*T);
            lnActCoeff_Scaled_[iK] += (delAK * XB + XA * delBK - XA * XB) * (g0 + g1 * XB) + XA * XB * (delBK - XB) * g1;
        }
    }
}

void PhaseCombo_Interaction::s_update_dlnActCoeff_dT() const
{
    doublereal T = temperature();
    dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
    d2lnActCoeffdT2_Scaled_.assign(m_kk, 0.0);
    for (size_t iK = 0; iK < m_kk; iK++) {
        for (size_t i = 0; i < numBinaryInteractions_; i++) {
            size_t iA = m_pSpecies_A_ij[i];
            size_t iB = m_pSpecies_B_ij[i];
            int delAK = 0;
            int delBK = 0;
            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }
            double XA = moleFractions_[iA];
            double XB = moleFractions_[iB];
            double g0 = -m_HE_b_ij[i] / (GasConstant*T*T);
            double g1 = -m_HE_c_ij[i] / (GasConstant*T*T);
            double temp = (delAK * XB + XA * delBK - XA * XB) * (g0 + g1 * XB) + XA * XB * (delBK - XB) * g1;
            dlnActCoeffdT_Scaled_[iK] += temp;
            d2lnActCoeffdT2_Scaled_[iK] -= 2.0 * temp / T;
        }
    }
}

void PhaseCombo_Interaction::getdlnActCoeffdT(doublereal* dlnActCoeffdT) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
}

void PhaseCombo_Interaction::getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        d2lnActCoeffdT2[k] = d2lnActCoeffdT2_Scaled_[k];
    }
}

void PhaseCombo_Interaction::getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
        doublereal* dlnActCoeffds) const
{
    doublereal T = temperature();
    s_update_dlnActCoeff_dT();

    for (size_t iK = 0; iK < m_kk; iK++) {
        // We never sample the end of the mole fraction domains
        double xx = std::max(moleFractions_[iK], SmallNumber);

        // First wipe out the ideal solution mixing term
        if (xx > SmallNumber) {
            dlnActCoeffds[iK] += - 1.0 / xx;
        }

        for (size_t i = 0; i < numBinaryInteractions_; i++) {
            size_t iA = m_pSpecies_A_ij[i];
            size_t iB = m_pSpecies_B_ij[i];
            int delAK = 0;
            int delBK = 0;

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            double XA = moleFractions_[iA];
            double XB = moleFractions_[iB];
            double dXA = dXds[iA];
            double dXB = dXds[iB];
            double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / (GasConstant*T);
            double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / (GasConstant*T);
            dlnActCoeffds[iK] += ((delBK-XB)*dXA + (delAK-XA)*dXB)*(g0+2*g1*XB) + (delBK-XB)*2*g1*XA*dXB
                                 + dlnActCoeffdT_Scaled_[iK]*dTds;
        }
    }
}

void PhaseCombo_Interaction::s_update_dlnActCoeff_dlnN_diag() const
{
    doublereal T = temperature();
    dlnActCoeffdlnN_diag_.assign(m_kk, 0.0);

    for (size_t iK = 0; iK < m_kk; iK++) {
        double XK = moleFractions_[iK];
        // We never sample the end of the mole fraction domains
        double xx = std::max(moleFractions_[iK], SmallNumber);

        // First wipe out the ideal solution mixing term
        if (xx > SmallNumber) {
            dlnActCoeffdlnN_diag_[iK] = - 1.0 + xx;
        }

        for (size_t i = 0; i < numBinaryInteractions_; i++) {
            size_t iA = m_pSpecies_A_ij[i];
            size_t iB = m_pSpecies_B_ij[i];
            int delAK = 0;
            int delBK = 0;

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            double XA = moleFractions_[iA];
            double XB = moleFractions_[iB];
            double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / (GasConstant*T);
            double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / (GasConstant*T);
            dlnActCoeffdlnN_diag_[iK] += 2*(delBK-XB)*(g0*(delAK-XA)+g1*(2*(delAK-XA)*XB+XA*(delBK-XB)));
        }
        dlnActCoeffdlnN_diag_[iK] = XK*dlnActCoeffdlnN_diag_[iK];
    }
}

void PhaseCombo_Interaction::s_update_dlnActCoeff_dlnN() const
{
    double T = temperature();
    dlnActCoeffdlnN_.zero();

    // Loop over the activity coefficient gamma_k
    for (size_t iK = 0; iK < m_kk; iK++) {
        // We never sample the end of the mole fraction domains
        double xx = std::max(moleFractions_[iK], SmallNumber);
        for (size_t iM = 0; iM < m_kk; iM++) {
            double XM = moleFractions_[iM];
            if (xx > SmallNumber) {
                double delKM = 0.0;
                if (iK == iM) {
                    delKM = 1.0;
                }
                // this gets multiplied by XM at the bottom
                dlnActCoeffdlnN_(iK,iM) += - delKM/XM + 1.0;
            }

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
                double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / (GasConstant*T);
                double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / (GasConstant*T);
                dlnActCoeffdlnN_(iK,iM) += g0*((delAM-XA)*(delBK-XB)+(delAK-XA)*(delBM-XB));
                dlnActCoeffdlnN_(iK,iM) += 2*g1*((delAM-XA)*(delBK-XB)*XB+(delAK-XA)*(delBM-XB)*XB+(delBM-XB)*(delBK-XB)*XA);
            }
            dlnActCoeffdlnN_(iK,iM) = XM * dlnActCoeffdlnN_(iK,iM);
        }
    }
}

void PhaseCombo_Interaction::s_update_dlnActCoeff_dlnX_diag() const
{
    doublereal T = temperature();
    dlnActCoeffdlnX_diag_.assign(m_kk, 0.0);
    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];

        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];

        double g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / (GasConstant * T);
        double g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / (GasConstant * T);

        dlnActCoeffdlnX_diag_[iA] += XA*XB*(2*g1*-2*g0-6*g1*XB);
        dlnActCoeffdlnX_diag_[iB] += XA*XB*(2*g1*-2*g0-6*g1*XB);
    }
    throw CanteraError("PhaseCombo_Interaction::s_update_dlnActCoeff_dlnX_diag", "unimplemented");
}

void PhaseCombo_Interaction::getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const
{
    s_update_dlnActCoeff_dlnN_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnN_diag[k] = dlnActCoeffdlnN_diag_[k];
    }
}

void PhaseCombo_Interaction::getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const
{
    s_update_dlnActCoeff_dlnX_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}

void PhaseCombo_Interaction::getdlnActCoeffdlnN(const size_t ld, doublereal* dlnActCoeffdlnN)
{
    s_update_dlnActCoeff_dlnN();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}

void PhaseCombo_Interaction::resizeNumInteractions(const size_t num)
{
    numBinaryInteractions_ = num;
    m_HE_b_ij.resize(num, 0.0);
    m_HE_c_ij.resize(num, 0.0);
    m_HE_d_ij.resize(num, 0.0);
    m_SE_b_ij.resize(num, 0.0);
    m_SE_c_ij.resize(num, 0.0);
    m_SE_d_ij.resize(num, 0.0);
    m_VHE_b_ij.resize(num, 0.0);
    m_VHE_c_ij.resize(num, 0.0);
    m_VHE_d_ij.resize(num, 0.0);
    m_VSE_b_ij.resize(num, 0.0);
    m_VSE_c_ij.resize(num, 0.0);
    m_VSE_d_ij.resize(num, 0.0);
    m_pSpecies_A_ij.resize(num, npos);
    m_pSpecies_B_ij.resize(num, npos);
}

void PhaseCombo_Interaction::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    string xname = xmLBinarySpecies.name();
    if (xname != "binaryNeutralSpeciesParameters") {
        throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies",
                           "Incorrect name for processing this routine: " + xname);
    }
    vector_fp vParams;
    string iName = xmLBinarySpecies.attrib("speciesA");
    if (iName == "") {
        throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies", "no speciesA attrib");
    }
    string jName = xmLBinarySpecies.attrib("speciesB");
    if (jName == "") {
        throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies", "no speciesB attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    string ispName = speciesName(iSpecies);
    if (charge(iSpecies) != 0) {
        throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies", "speciesA charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    string jspName = speciesName(jSpecies);
    if (charge(jSpecies) != 0) {
        throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies", "speciesB charge problem");
    }

    resizeNumInteractions(numBinaryInteractions_ + 1);
    size_t iSpot = numBinaryInteractions_ - 1;
    m_pSpecies_A_ij[iSpot] = iSpecies;
    m_pSpecies_B_ij[iSpot] = jSpecies;

    for (size_t iChild = 0; iChild < xmLBinarySpecies.nChildren(); iChild++) {
        XML_Node& xmlChild = xmLBinarySpecies.child(iChild);
        string nodeName = toLowerCopy(xmlChild.name());

        // Process the binary species interaction child elements
        if (nodeName == "excessenthalpy") {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, true, "toSI", "excessEnthalpy");
            if (vParams.size() != 2) {
                throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies::excessEnthalpy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_HE_b_ij[iSpot] = vParams[0];
            m_HE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessentropy") {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, true, "toSI", "excessEntropy");
            if (vParams.size() != 2) {
                throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies::excessEntropy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_SE_b_ij[iSpot] = vParams[0];
            m_SE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessvolume_enthalpy") {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, true, "toSI", "excessVolume_Enthalpy");
            if (vParams.size() != 2) {
                throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies::excessVolume_Enthalpy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_VHE_b_ij[iSpot] = vParams[0];
            m_VHE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessvolume_entropy") {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, true, "toSI", "excessVolume_Entropy");
            if (vParams.size() != 2) {
                throw CanteraError("PhaseCombo_Interaction::readXMLBinarySpecies::excessVolume_Entropy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found");
            }
            m_VSE_b_ij[iSpot] = vParams[0];
            m_VSE_c_ij[iSpot] = vParams[1];
        }
    }
}

}
