/**
 *  @file RedlichKisterVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess Gibbs free energy formulations related to RedlichKister
 *   expansions (see \ref thermoprops
 *    and class \link Cantera::RedlichKisterVPSSTP RedlichKisterVPSSTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{
RedlichKisterVPSSTP::RedlichKisterVPSSTP() :
    numBinaryInteractions_(0),
    formRedlichKister_(0),
    formTempModel_(0)
{
}

RedlichKisterVPSSTP::RedlichKisterVPSSTP(const std::string& inputFile,
        const std::string& id_) :
    numBinaryInteractions_(0),
    formRedlichKister_(0),
    formTempModel_(0)
{
    initThermoFile(inputFile, id_);
}

RedlichKisterVPSSTP::RedlichKisterVPSSTP(XML_Node& phaseRoot,
        const std::string& id_) :
    numBinaryInteractions_(0),
    formRedlichKister_(0),
    formTempModel_(0)
{
    importPhase(phaseRoot, this);
}

// - Activities, Standard States, Activity Concentrations -----------

void RedlichKisterVPSSTP::getLnActivityCoefficients(doublereal* lnac) const
{
    // Update the activity coefficients
    s_update_lnActCoeff();

    for (size_t k = 0; k < m_kk; k++) {
        lnac[k] = lnActCoeff_Scaled_[k];
    }
}

// ------------ Partial Molar Properties of the Solution ------------

void RedlichKisterVPSSTP::getChemPotentials(doublereal* mu) const
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

doublereal RedlichKisterVPSSTP::enthalpy_mole() const
{
    double h = 0;
    vector_fp hbar(m_kk);
    getPartialMolarEnthalpies(&hbar[0]);
    for (size_t i = 0; i < m_kk; i++) {
        h += moleFractions_[i]*hbar[i];
    }
    return h;
}

doublereal RedlichKisterVPSSTP::entropy_mole() const
{
    double s = 0;
    vector_fp sbar(m_kk);
    getPartialMolarEntropies(&sbar[0]);
    for (size_t i = 0; i < m_kk; i++) {
        s += moleFractions_[i]*sbar[i];
    }
    return s;
}

doublereal RedlichKisterVPSSTP::cp_mole() const
{
    double cp = 0;
    vector_fp cpbar(m_kk);
    getPartialMolarCp(&cpbar[0]);
    for (size_t i = 0; i < m_kk; i++) {
        cp += moleFractions_[i]*cpbar[i];
    }
    return cp;
}

doublereal RedlichKisterVPSSTP::cv_mole() const
{
    return cp_mole() - GasConstant;
}

void RedlichKisterVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
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

void RedlichKisterVPSSTP::getPartialMolarCp(doublereal* cpbar) const
{
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

void RedlichKisterVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
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

void RedlichKisterVPSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    // Get the standard state values in m^3 kmol-1
    getStandardVolumes(vbar);
    for (size_t iK = 0; iK < m_kk; iK++) {
        vbar[iK] += 0.0;
    }
}

void RedlichKisterVPSSTP::initThermo()
{
    if (m_input.hasKey("interactions")) {
        for (const auto& item : m_input["interactions"].asVector<AnyMap>()) {
            auto& species = item["species"].asVector<string>(2);
            vector_fp h_excess = item.convertVector("excess-enthalpy", "J/kmol");
            vector_fp s_excess = item.convertVector("excess-entropy", "J/kmol/K");
            addBinaryInteraction(species[0], species[1],
                                 h_excess.data(), h_excess.size(),
                                 s_excess.data(), s_excess.size());
        }
    }
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}

void RedlichKisterVPSSTP::initLengths()
{
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
}

void RedlichKisterVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if ((int) id_.size() > 0 && phaseNode.id() != id_) {
        throw CanteraError("RedlichKisterVPSSTP::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    // Check on the thermo field. Must have:
    // <thermo model="Redlich-Kister" />
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("RedlichKisterVPSSTP::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");
    if (!caseInsensitiveEquals(thermoNode["model"], "redlich-kister")) {
        throw CanteraError("RedlichKisterVPSSTP::initThermoXML",
                           "Unknown thermo model: " + thermoNode["model"]
                           + " - This object only knows \"Redlich-Kister\" ");
    }

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        if (!caseInsensitiveEquals(acNode["model"], "redlich-kister")) {
            throw CanteraError("RedlichKisterVPSSTP::initThermoXML",
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
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id_);
}

void RedlichKisterVPSSTP::s_update_lnActCoeff() const
{
    doublereal T = temperature();
    lnActCoeff_Scaled_.assign(m_kk, 0.0);

    // Scaling: I moved the division of RT higher so that we are always dealing
    // with G/RT dimensionless terms within the routine. There is a severe
    // problem with roundoff error in these calculations. The dimensionless
    // terms help.
    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        doublereal deltaX = XA - XB;
        size_t N = m_N_ij[i];
        vector_fp& he_vec = m_HE_m_ij[i];
        vector_fp& se_vec = m_SE_m_ij[i];
        doublereal poly = 1.0;
        doublereal polyMm1 = 1.0;
        doublereal sum = 0.0;
        doublereal sumMm1 = 0.0;
        doublereal sum2 = 0.0;
        for (size_t m = 0; m < N; m++) {
            doublereal A_ge = (he_vec[m] - T * se_vec[m]) / (GasConstant * T);
            sum += A_ge * poly;
            sum2 += A_ge * (m + 1) * poly;
            poly *= deltaX;
            if (m >= 1) {
                sumMm1 += (A_ge * polyMm1 * m);
                polyMm1 *= deltaX;
            }
        }
        doublereal oneMXA = 1.0 - XA;
        doublereal oneMXB = 1.0 - XB;
        for (size_t k = 0; k < m_kk; k++) {
            if (iA == k) {
                lnActCoeff_Scaled_[k] += (oneMXA * XB * sum) + (XA * XB * sumMm1 * (oneMXA + XB));
            } else if (iB == k) {
                lnActCoeff_Scaled_[k] += (oneMXB * XA * sum) + (XA * XB * sumMm1 * (-oneMXB - XA));
            } else {
                lnActCoeff_Scaled_[k] += -(XA * XB * sum2);
            }
        }
        // Debug against formula in literature
    }
}

void RedlichKisterVPSSTP::s_update_dlnActCoeff_dT() const
{
    dlnActCoeffdT_Scaled_.assign(m_kk, 0.0);
    d2lnActCoeffdT2_Scaled_.assign(m_kk, 0.0);

    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        doublereal deltaX = XA - XB;
        size_t N = m_N_ij[i];
        doublereal poly = 1.0;
        doublereal sum = 0.0;
        vector_fp& se_vec = m_SE_m_ij[i];
        doublereal sumMm1 = 0.0;
        doublereal polyMm1 = 1.0;
        doublereal sum2 = 0.0;
        for (size_t m = 0; m < N; m++) {
            doublereal A_ge = - se_vec[m];
            sum += A_ge * poly;
            sum2 += A_ge * (m + 1) * poly;
            poly *= deltaX;
            if (m >= 1) {
                sumMm1 += (A_ge * polyMm1 * m);
                polyMm1 *= deltaX;
            }
        }
        doublereal oneMXA = 1.0 - XA;
        doublereal oneMXB = 1.0 - XB;
        for (size_t k = 0; k < m_kk; k++) {
            if (iA == k) {
                dlnActCoeffdT_Scaled_[k] += (oneMXA * XB * sum) + (XA * XB * sumMm1 * (oneMXA + XB));
            } else if (iB == k) {
                dlnActCoeffdT_Scaled_[k] += (oneMXB * XA * sum) + (XA * XB * sumMm1 * (-oneMXB - XA));
            } else {
                dlnActCoeffdT_Scaled_[k] += -(XA * XB * sum2);
            }
        }
    }
}

void RedlichKisterVPSSTP::getdlnActCoeffdT(doublereal* dlnActCoeffdT) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
}

void RedlichKisterVPSSTP::getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        d2lnActCoeffdT2[k] = d2lnActCoeffdT2_Scaled_[k];
    }
}

void RedlichKisterVPSSTP::s_update_dlnActCoeff_dlnX_diag() const
{
    double T = temperature();
    dlnActCoeffdlnX_diag_.assign(m_kk, 0.0);

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
      size_t iA =  m_pSpecies_A_ij[i];
      size_t iB =  m_pSpecies_B_ij[i];
      double XA = moleFractions_[iA];
      double XB = moleFractions_[iB];
      double deltaX = XA - XB;
      size_t N = m_N_ij[i];
      double poly = 1.0;
      double sum = 0.0;
      vector_fp& he_vec = m_HE_m_ij[i];
      vector_fp& se_vec = m_SE_m_ij[i];
      double sumMm1 = 0.0;
      double polyMm1 = 1.0;
      double polyMm2 = 1.0;
      double sumMm2 = 0.0;
      for (size_t m = 0; m < N; m++) {
          double A_ge = (he_vec[m] -  T * se_vec[m]) / (GasConstant * T);;
          sum += A_ge * poly;
          poly *= deltaX;
          if (m >= 1) {
              sumMm1  += (A_ge * polyMm1 * m);
              polyMm1 *= deltaX;
          }
          if (m >= 2) {
              sumMm2 += (A_ge * polyMm2 * m * (m - 1.0));
              polyMm2 *= deltaX;
          }
      }

      for (size_t k = 0; k < m_kk; k++) {
          if (iA == k) {
              dlnActCoeffdlnX_diag_[k] +=
                  XA * (- (1-XA+XB) * sum + 2*(1.0 - XA) * XB * sumMm1
                        + sumMm1 * (XB * (1 - 2*XA + XB) - XA * (1 - XA + 2*XB))
                        + 2 * XA * XB * sumMm2 * (1.0 - XA + XB));
          } else  if (iB == k) {
              dlnActCoeffdlnX_diag_[k] +=
                  XB * (- (1-XB+XA) * sum - 2*(1.0 - XB) * XA * sumMm1
                        + sumMm1 * (XA * (2*XB - XA - 1) - XB * (-2*XA + XB - 1))
                        - 2 * XA * XB * sumMm2 * (-XA - 1 + XB));
          }
      }
    }
}

void RedlichKisterVPSSTP::s_update_dlnActCoeff_dX_() const
{
    doublereal T = temperature();
    dlnActCoeff_dX_.zero();

    for (size_t i = 0; i < numBinaryInteractions_; i++) {
        size_t iA = m_pSpecies_A_ij[i];
        size_t iB = m_pSpecies_B_ij[i];
        double XA = moleFractions_[iA];
        double XB = moleFractions_[iB];
        doublereal deltaX = XA - XB;
        size_t N = m_N_ij[i];
        doublereal poly = 1.0;
        doublereal sum = 0.0;
        vector_fp& he_vec = m_HE_m_ij[i];
        vector_fp& se_vec = m_SE_m_ij[i];
        doublereal sumMm1 = 0.0;
        doublereal polyMm1 = 1.0;
        doublereal polyMm2 = 1.0;
        doublereal sum2 = 0.0;
        doublereal sum2Mm1 = 0.0;
        doublereal sumMm2 = 0.0;
        for (size_t m = 0; m < N; m++) {
            doublereal A_ge = he_vec[m] - T * se_vec[m];
            sum += A_ge * poly;
            sum2 += A_ge * (m + 1) * poly;
            poly *= deltaX;
            if (m >= 1) {
                sumMm1 += (A_ge * polyMm1 * m);
                sum2Mm1 += (A_ge * polyMm1 * m * (1.0 + m));
                polyMm1 *= deltaX;
            }
            if (m >= 2) {
                sumMm2 += (A_ge * polyMm2 * m * (m - 1.0));
                polyMm2 *= deltaX;
            }
        }

        for (size_t k = 0; k < m_kk; k++) {
            if (iA == k) {
                dlnActCoeff_dX_(k, iA) += (- XB * sum + (1.0 - XA) * XB * sumMm1
                                           + XB * sumMm1 * (1.0 - 2.0 * XA + XB)
                                           + XA * XB * sumMm2 * (1.0 - XA + XB));

                dlnActCoeff_dX_(k, iB) += ((1.0 - XA) * sum - (1.0 - XA) * XB * sumMm1
                                           + XA * sumMm1 * (1.0 + 2.0 * XB - XA)
                                           - XA * XB * sumMm2 * (1.0 - XA + XB));
            } else if (iB == k) {
                dlnActCoeff_dX_(k, iA) += ((1.0 - XB) * sum + (1.0 - XA) * XB * sumMm1
                                           + XB * sumMm1 * (1.0 - 2.0 * XA + XB)
                                           + XA * XB * sumMm2 * (1.0 - XA + XB));

                dlnActCoeff_dX_(k, iB) += (- XA * sum - (1.0 - XB) * XA * sumMm1
                                           + XA * sumMm1 * (XB - XA - (1.0 - XB))
                                           - XA * XB * sumMm2 * (-XA - (1.0 - XB)));
            } else {
                dlnActCoeff_dX_(k, iA) += (- XB * sum2 - XA * XB * sum2Mm1);
                dlnActCoeff_dX_(k, iB) += (- XA * sum2 + XA * XB * sum2Mm1);
            }
        }
    }
}

void RedlichKisterVPSSTP::getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
        doublereal* dlnActCoeffds) const
{
    s_update_dlnActCoeff_dT();
    s_update_dlnActCoeff_dX_();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffds[k] = dlnActCoeffdT_Scaled_[k] * dTds;
        for (size_t j = 0; j < m_kk; j++) {
            dlnActCoeffds[k] += dlnActCoeff_dX_(k, j) * dXds[j];
        }
    }
}

void RedlichKisterVPSSTP::getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const
{
    s_update_dlnActCoeff_dX_();
    for (size_t j = 0; j < m_kk; j++) {
        dlnActCoeffdlnN_diag[j] = dlnActCoeff_dX_(j, j);
        for (size_t k = 0; k < m_kk; k++) {
            dlnActCoeffdlnN_diag[k] -= dlnActCoeff_dX_(j, k) * moleFractions_[k];
        }
    }
}

void RedlichKisterVPSSTP::getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const
{
    s_update_dlnActCoeff_dlnX_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}

void RedlichKisterVPSSTP::getdlnActCoeffdlnN(const size_t ld, doublereal* dlnActCoeffdlnN)
{
    s_update_dlnActCoeff_dX_();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}

void RedlichKisterVPSSTP::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    std::string xname = xmLBinarySpecies.name();
    if (xname != "binaryNeutralSpeciesParameters") {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies",
                           "Incorrect name for processing this routine: " + xname);
    }
    vector_fp hParams, sParams;
    std::string iName = xmLBinarySpecies.attrib("speciesA");
    if (iName == "") {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies", "no speciesA attrib");
    }
    std::string jName = xmLBinarySpecies.attrib("speciesB");
    if (jName == "") {
        throw CanteraError("RedlichKisterVPSSTP::readXMLBinarySpecies", "no speciesB attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species. This means that the interaction doesn't occur for
    // the current implementation of the phase.
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }

    // Ok we have found a valid interaction
    for (size_t iChild = 0; iChild < xmLBinarySpecies.nChildren(); iChild++) {
        XML_Node& xmlChild = xmLBinarySpecies.child(iChild);
        string nodeName = toLowerCopy(xmlChild.name());

        // Process the binary species interaction child elements
        if (nodeName == "excessenthalpy") {
            getFloatArray(xmlChild, hParams, true, "toSI", "excessEnthalpy");
        } else if (nodeName == "excessentropy") {
            getFloatArray(xmlChild, sParams, true, "toSI", "excessEntropy");
        }
    }
    addBinaryInteraction(iName, jName, hParams.data(), hParams.size(),
                         sParams.data(), sParams.size());
}

void RedlichKisterVPSSTP::addBinaryInteraction(
    const std::string& speciesA, const std::string& speciesB,
    const double* excess_enthalpy, size_t n_enthalpy,
    const double* excess_entropy, size_t n_entropy)
{
    size_t kA = speciesIndex(speciesA);
    size_t kB = speciesIndex(speciesB);
    if (kA == npos) {
        throw CanteraError("RedlichKisterVPSSTP::addBinaryInteraction",
            "Species '{}' not present in phase", speciesA);
    } else if (kB == npos) {
        throw CanteraError("RedlichKisterVPSSTP::addBinaryInteraction",
            "Species '{}' not present in phase", speciesB);
    }
    if (charge(kA) != 0) {
        throw CanteraError("RedlichKisterVPSSTP::addBinaryInteraction",
            "Species '{}' should be neutral", speciesA);
    } else if (charge(kB) != 0) {
        throw CanteraError("RedlichKisterVPSSTP::addBinaryInteraction",
            "Species '{}' should be neutral", speciesB);
    }

    m_pSpecies_A_ij.push_back(kA);
    m_pSpecies_B_ij.push_back(kB);
    m_HE_m_ij.emplace_back(excess_enthalpy, excess_enthalpy + n_enthalpy);
    m_SE_m_ij.emplace_back(excess_entropy, excess_entropy + n_entropy);
    size_t N = max(n_enthalpy, n_entropy);
    m_HE_m_ij.back().resize(N, 0.0);
    m_SE_m_ij.back().resize(N, 0.0);
    m_N_ij.push_back(N);
    dlnActCoeff_dX_.resize(N, N, 0.0);
    numBinaryInteractions_++;
}

}
