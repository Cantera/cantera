/**
 *  @file HMWSoln_input.cpp
 *    Definitions for the HMWSoln ThermoPhase object, which models concentrated
 *    electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::HMWSoln HMWSoln \endlink) .
 *
 * This file contains definitions for reading in the interaction terms
 * in the formulation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/electrolytes.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

#include <fstream>

using namespace std;

namespace Cantera
{
int HMWSoln::interp_est(const std::string& estString)
{
    if (ba::iequals(estString, "solvent")) {
        return cEST_solvent;
    } else if (ba::iequals(estString, "chargedspecies")) {
        return cEST_chargedSpecies;
    } else if (ba::iequals(estString, "weakacidassociated")) {
        return cEST_weakAcidAssociated;
    } else if (ba::iequals(estString, "strongacidassociated")) {
        return cEST_strongAcidAssociated;
    } else if (ba::iequals(estString, "polarneutral")) {
        return cEST_polarNeutral;
    } else if (ba::iequals(estString, "nonpolarneutral")) {
        return cEST_nonpolarNeutral;
    }
    int retn, rval;
    if ((retn = sscanf(estString.c_str(), "%d", &rval)) != 1) {
        return -1;
    }
    return rval;
}

void HMWSoln::readXMLBinarySalt(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    if (xname != "binarySaltParameters") {
        throw CanteraError("HMWSoln::readXMLBinarySalt",
                           "Incorrect name for processing this routine: " + xname);
    }

    vector_fp vParams;
    string iName = BinSalt.attrib("cation");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "no cation attrib");
    }
    string jName = BinSalt.attrib("anion");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "no anion attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    string ispName = speciesName(iSpecies);
    if (charge(iSpecies) <= 0) {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "cation charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    string jspName = speciesName(jSpecies);
    if (charge(jSpecies) >= 0) {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "anion charge problem");
    }

    size_t n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    for (size_t iChild = 0; iChild < BinSalt.nChildren(); iChild++) {
        XML_Node& xmlChild = BinSalt.child(iChild);
        string nodeName = xmlChild.name();

        // Process the binary salt child elements
        if (ba::iequals(nodeName, "beta0")) {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, false, "", "beta0");
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta0 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Beta0MX_ij[counter] = vParams[0];
                m_Beta0MX_ij_coeff(0,counter) = m_Beta0MX_ij[counter];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta0 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Beta0MX_ij_coeff(0,counter) = vParams[0];
                m_Beta0MX_ij_coeff(1,counter) = vParams[1];
                m_Beta0MX_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta0 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                for (size_t i = 0; i < nParamsFound; i++) {
                    m_Beta0MX_ij_coeff(i, counter) = vParams[i];
                }
                m_Beta0MX_ij[counter] = vParams[0];
            }
        }
        if (ba::iequals(nodeName, "beta1")) {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, false, "", "beta1");
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta1 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Beta1MX_ij[counter] = vParams[0];
                m_Beta1MX_ij_coeff(0,counter) = m_Beta1MX_ij[counter];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta1 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Beta1MX_ij_coeff(0,counter) = vParams[0];
                m_Beta1MX_ij_coeff(1,counter) = vParams[1];
                m_Beta1MX_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta1 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                for (size_t i = 0; i < nParamsFound; i++) {
                    m_Beta1MX_ij_coeff(i, counter) = vParams[i];
                }
                m_Beta1MX_ij[counter] = vParams[0];
            }
        }
        if (ba::iequals(nodeName, "beta2")) {
            getFloatArray(xmlChild, vParams, false, "", "beta2");
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta2 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Beta2MX_ij[counter] = vParams[0];
                m_Beta2MX_ij_coeff(0,counter) = m_Beta2MX_ij[counter];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta2 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Beta2MX_ij_coeff(0,counter) = vParams[0];
                m_Beta2MX_ij_coeff(1,counter) = vParams[1];
                m_Beta2MX_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::beta2 for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                for (size_t i = 0; i < nParamsFound; i++) {
                    m_Beta2MX_ij_coeff(i, counter) = vParams[i];
                }
                m_Beta2MX_ij[counter] = vParams[0];
            }
        }
        if (ba::iequals(nodeName, "cphi")) {
            // Get the string containing all of the values
            getFloatArray(xmlChild, vParams, false, "", "Cphi");
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::Cphi for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_CphiMX_ij[counter] = vParams[0];
                m_CphiMX_ij_coeff(0,counter) = m_CphiMX_ij[counter];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::Cphi for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_CphiMX_ij_coeff(0,counter) = vParams[0];
                m_CphiMX_ij_coeff(1,counter) = vParams[1];
                m_CphiMX_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLBinarySalt::Cphi for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                for (size_t i = 0; i < nParamsFound; i++) {
                    m_CphiMX_ij_coeff(i, counter) = vParams[i];
                }
                m_CphiMX_ij[counter] = vParams[0];
            }
        }

        if (ba::iequals(nodeName, "alpha1")) {
            m_Alpha1MX_ij[counter] = fpValueCheck(xmlChild.value());
        }

        if (ba::iequals(nodeName, "alpha2")) {
            m_Alpha2MX_ij[counter] = fpValueCheck(xmlChild.value());
        }
    }
}

void HMWSoln::readXMLThetaAnion(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    vector_fp vParams;
    if (xname != "thetaAnion") {
        throw CanteraError("HMWSoln::readXMLThetaAnion",
                           "Incorrect name for processing this routine: " + xname);
    }
    string ispName = BinSalt.attrib("anion1");
    if (ispName == "") {
        throw CanteraError("HMWSoln::readXMLThetaAnion", "no anion1 attrib");
    }
    string jspName = BinSalt.attrib("anion2");
    if (jspName == "") {
        throw CanteraError("HMWSoln::readXMLThetaAnion", "no anion2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(ispName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) >= 0) {
        throw CanteraError("HMWSoln::readXMLThetaAnion", "anion1 charge problem");
    }
    size_t jSpecies = speciesIndex(jspName);
    if (jSpecies == npos) {
        return;
    }
    if (charge(jSpecies) >= 0) {
        throw CanteraError("HMWSoln::readXMLThetaAnion", "anion2 charge problem");
    }

    size_t n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "theta")) {
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLThetaAnion::Theta for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Theta_ij_coeff(0,counter) = vParams[0];
                m_Theta_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLThetaAnion::Theta for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Theta_ij_coeff(0,counter) = vParams[0];
                m_Theta_ij_coeff(1,counter) = vParams[1];
                m_Theta_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLThetaAnion::Theta for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Theta_ij_coeff(j, counter) = vParams[j];
                }
                m_Theta_ij[counter] = vParams[0];
            }
        }
    }
}

void HMWSoln::readXMLThetaCation(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    vector_fp vParams;
    if (xname != "thetaCation") {
        throw CanteraError("HMWSoln::readXMLThetaCation",
                           "Incorrect name for processing this routine: " + xname);
    }
    string ispName = BinSalt.attrib("cation1");
    if (ispName == "") {
        throw CanteraError("HMWSoln::readXMLThetaCation", "no cation1 attrib");
    }
    string jspName = BinSalt.attrib("cation2");
    if (jspName == "") {
        throw CanteraError("HMWSoln::readXMLThetaCation", "no cation2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(ispName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) <= 0) {
        throw CanteraError("HMWSoln::readXMLThetaCation", "cation1 charge problem");
    }
    size_t jSpecies = speciesIndex(jspName);
    if (jSpecies == npos) {
        return;
    }
    if (charge(jSpecies) <= 0) {
        throw CanteraError("HMWSoln::readXMLThetaCation", "cation2 charge problem");
    }

    size_t n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "theta")) {
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLThetaCation::Theta for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Theta_ij_coeff(0,counter) = vParams[0];
                m_Theta_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLThetaCation::Theta for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                m_Theta_ij_coeff(0,counter) = vParams[0];
                m_Theta_ij_coeff(1,counter) = vParams[1];
                m_Theta_ij[counter] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLThetaCation::Theta for " + ispName
                                       + "::" + jspName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Theta_ij_coeff(j, counter) = vParams[j];
                }
                m_Theta_ij[counter] = vParams[0];
            }
        }
    }
}

void HMWSoln::readXMLPsiCommonCation(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    if (xname != "psiCommonCation") {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation",
                           "Incorrect name for processing this routine: " + xname);
    }
    vector_fp vParams;
    string kName = BinSalt.attrib("cation");
    if (kName == "") {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation", "no cation attrib");
    }
    string iName = BinSalt.attrib("anion1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation", "no anion1 attrib");
    }
    string jName = BinSalt.attrib("anion2");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation", "no anion2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t kSpecies = speciesIndex(kName);
    if (kSpecies == npos) {
        return;
    }
    if (charge(kSpecies) <= 0) {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation",
                           "cation charge problem");
    }
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) >= 0) {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation",
                           "anion1 charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    if (charge(jSpecies) >= 0) {
        throw CanteraError("HMWSoln::readXMLPsiCommonCation",
                           "anion2 charge problem");
    }

    size_t n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "theta")) {
            double old = m_Theta_ij[counter];
            m_Theta_ij[counter] = fpValueCheck(xmlChild.value());
            if (old != 0.0 && old != m_Theta_ij[counter]) {
                throw CanteraError("HMWSoln::readXMLPsiCommonCation",
                                   "conflicting values");
            }
        }
        if (ba::iequals(xmlChild.name(), "psi")) {
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            n = iSpecies * m_kk *m_kk + jSpecies * m_kk + kSpecies;

            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLPsiCommonCation::Psi for "
                                       + kName + "::" + iName + "::" + jName,
                                       "wrong number of params found");
                }
                m_Psi_ijk_coeff(0,n) = vParams[0];
                m_Psi_ijk[n] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLPsiCation::Psi for "
                                       + kName + "::" + iName + "::" + jName,
                                       "wrong number of params found");
                }
                m_Psi_ijk_coeff(0,n) = vParams[0];
                m_Psi_ijk_coeff(1,n) = vParams[1];
                m_Psi_ijk[n] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLPsiCation::Psi for "
                                       + kName + "::" + iName + "::" + jName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Psi_ijk_coeff(j, n) = vParams[j];
                }
                m_Psi_ijk[n] = vParams[0];
            }

            // fill in the duplicate entries
            n = iSpecies * m_kk *m_kk + kSpecies * m_kk + jSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = jSpecies * m_kk *m_kk + iSpecies * m_kk + kSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = jSpecies * m_kk *m_kk + kSpecies * m_kk + iSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = kSpecies * m_kk *m_kk + jSpecies * m_kk + iSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = kSpecies * m_kk *m_kk + iSpecies * m_kk + jSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];
        }
    }
}

void HMWSoln::readXMLPsiCommonAnion(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    if (xname != "psiCommonAnion") {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion",
                           "Incorrect name for processing this routine: " + xname);
    }
    vector_fp vParams;
    string kName = BinSalt.attrib("anion");
    if (kName == "") {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "no anion attrib");
    }
    string iName = BinSalt.attrib("cation1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "no cation1 attrib");
    }
    string jName = BinSalt.attrib("cation2");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "no cation2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t kSpecies = speciesIndex(kName);
    if (kSpecies == npos) {
        return;
    }
    if (charge(kSpecies) >= 0) {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion", "anion charge problem");
    }
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) <= 0) {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion",
                           "cation1 charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    if (charge(jSpecies) <= 0) {
        throw CanteraError("HMWSoln::readXMLPsiCommonAnion",
                           "cation2 charge problem");
    }

    size_t n = iSpecies * m_kk + jSpecies;
    int counter = m_CounterIJ[n];
    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "theta")) {
            double old = m_Theta_ij[counter];
            m_Theta_ij[counter] = fpValueCheck(xmlChild.value());
            if (old != 0.0 && old != m_Theta_ij[counter]) {
                throw CanteraError("HMWSoln::readXMLPsiCommonAnion",
                                   "conflicting values");
            }
        }
        if (ba::iequals(xmlChild.name(), "psi")) {
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            n = iSpecies * m_kk *m_kk + jSpecies * m_kk + kSpecies;

            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLPsiCommonAnion::Psi for "
                                       + kName + "::" + iName + "::" + jName,
                                       "wrong number of params found");
                }
                m_Psi_ijk_coeff(0,n) = vParams[0];
                m_Psi_ijk[n] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLPsiAnion::Psi for "
                                       + kName + "::" + iName + "::" + jName,
                                       "wrong number of params found");
                }
                m_Psi_ijk_coeff(0,n) = vParams[0];
                m_Psi_ijk_coeff(1,n) = vParams[1];
                m_Psi_ijk[n] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLPsiAnion::Psi for "
                                       + kName + "::" + iName + "::" + jName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Psi_ijk_coeff(j, n) = vParams[j];
                }
                m_Psi_ijk[n] = vParams[0];
            }

            // fill in the duplicate entries
            n = iSpecies * m_kk *m_kk + kSpecies * m_kk + jSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = jSpecies * m_kk *m_kk + iSpecies * m_kk + kSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = jSpecies * m_kk *m_kk + kSpecies * m_kk + iSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = kSpecies * m_kk *m_kk + jSpecies * m_kk + iSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];

            n = kSpecies * m_kk *m_kk + iSpecies * m_kk + jSpecies;
            for (size_t j = 0; j < nParamsFound; j++) {
                m_Psi_ijk_coeff(j, n) = vParams[j];
            }
            m_Psi_ijk[n] = vParams[0];
        }
    }
}

void HMWSoln::readXMLLambdaNeutral(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    vector_fp vParams;
    if (xname != "lambdaNeutral") {
        throw CanteraError("HMWSoln::readXMLLanbdaNeutral",
                           "Incorrect name for processing this routine: " + xname);
    }
    string stemp;
    string iName = BinSalt.attrib("species1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no species1 attrib");
    }
    string jName = BinSalt.attrib("species2");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no species2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) != 0) {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral",
                           "neutral charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }

    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "lambda")) {
            size_t nCount = iSpecies*m_kk + jSpecies;
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLLambdaNeutral::Lambda for " + iName
                                       + "::" + jName,
                                       "wrong number of params found");
                }
                m_Lambda_nj_coeff(0,nCount) = vParams[0];
                m_Lambda_nj(iSpecies,jSpecies) = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLLambdaNeutral::Lambda for " + iName
                                       + "::" + jName,
                                       "wrong number of params found");
                }
                m_Lambda_nj_coeff(0,nCount) = vParams[0];
                m_Lambda_nj_coeff(1,nCount) = vParams[1];
                m_Lambda_nj(iSpecies, jSpecies) = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLLambdaNeutral::Lambda for " + iName
                                       + "::" + jName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Lambda_nj_coeff(j,nCount) = vParams[j];
                }
                m_Lambda_nj(iSpecies, jSpecies) = vParams[0];
            }
        }
    }
}

void HMWSoln::readXMLMunnnNeutral(XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    vector_fp vParams;
    if (xname != "MunnnNeutral") {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral",
                           "Incorrect name for processing this routine: " + xname);
    }
    string iName = BinSalt.attrib("species1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral", "no species1 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) != 0) {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral",
                           "neutral charge problem");
    }

    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "Munnn")) {
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLMunnnNeutral::Munnn for " + iName,
                                       "wrong number of params found");
                }
                m_Mu_nnn_coeff(0,iSpecies) = vParams[0];
                m_Mu_nnn[iSpecies] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLMunnnNeutral::Munnn for " + iName,
                                       "wrong number of params found");
                }
                m_Mu_nnn_coeff(0, iSpecies) = vParams[0];
                m_Mu_nnn_coeff(1, iSpecies) = vParams[1];
                m_Mu_nnn[iSpecies] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLMunnnNeutral::Munnn for " + iName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Mu_nnn_coeff(j, iSpecies) = vParams[j];
                }
                m_Mu_nnn[iSpecies] = vParams[0];
            }
        }
    }
}

void HMWSoln::readXMLZetaCation(const XML_Node& BinSalt)
{
    string xname = BinSalt.name();
    if (xname != "zetaCation") {
        throw CanteraError("HMWSoln::readXMLZetaCation",
                           "Incorrect name for processing this routine: " + xname);
    }
    vector_fp vParams;

    string iName = BinSalt.attrib("neutral");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no neutral attrib");
    }

    string jName = BinSalt.attrib("cation1");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no cation1 attrib");
    }

    string kName = BinSalt.attrib("anion1");
    if (kName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no anion1 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    if (charge(iSpecies) != 0.0) {
        throw CanteraError("HMWSoln::readXMLZetaCation", "neutral charge problem");
    }

    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    if (charge(jSpecies) <= 0.0) {
        throw CanteraError("HMWSoln::readXLZetaCation", "cation1 charge problem");
    }

    size_t kSpecies = speciesIndex(kName);
    if (kSpecies == npos) {
        return;
    }
    if (charge(kSpecies) >= 0.0) {
        throw CanteraError("HMWSoln::readXMLZetaCation", "anion1 charge problem");
    }

    for (size_t i = 0; i < BinSalt.nChildren(); i++) {
        XML_Node& xmlChild = BinSalt.child(i);
        if (ba::iequals(xmlChild.name(), "zeta")) {
            getFloatArray(xmlChild, vParams, false, "", xmlChild.name());
            size_t nParamsFound = vParams.size();
            size_t n = iSpecies * m_kk *m_kk + jSpecies * m_kk + kSpecies;

            if (m_formPitzerTemp == PITZER_TEMP_CONSTANT) {
                if (nParamsFound != 1) {
                    throw CanteraError("HMWSoln::readXMLZetaCation::Zeta for "
                                       + iName + "::" + jName + "::" + kName,
                                       "wrong number of params found");
                }
                m_Psi_ijk_coeff(0,n) = vParams[0];
                m_Psi_ijk[n] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR) {
                if (nParamsFound != 2) {
                    throw CanteraError("HMWSoln::readXMLZetaCation::Zeta for "
                                       + iName + "::" + jName + "::" + kName,
                                       "wrong number of params found");
                }
                m_Psi_ijk_coeff(0,n) = vParams[0];
                m_Psi_ijk_coeff(1,n) = vParams[1];
                m_Psi_ijk[n] = vParams[0];
            } else if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
                if (nParamsFound == 1) {
                    vParams.resize(5, 0.0);
                    nParamsFound = 5;
                } else if (nParamsFound != 5) {
                    throw CanteraError("HMWSoln::readXMLZetaCation::Zeta for "
                                       + iName + "::" + jName + "::" + kName,
                                       "wrong number of params found");
                }
                for (size_t j = 0; j < nParamsFound; j++) {
                    m_Psi_ijk_coeff(j, n) = vParams[j];
                }
                m_Psi_ijk[n] = vParams[0];
            }
            // There are no duplicate entries
        }
    }
}

void HMWSoln::readXMLCroppingCoefficients(const XML_Node& acNode)
{
    if (acNode.hasChild("croppingCoefficients")) {
        XML_Node& cropNode = acNode.child("croppingCoefficients");
        if (cropNode.hasChild("ln_gamma_k_min")) {
            XML_Node& gkminNode = cropNode.child("ln_gamma_k_min");
            getOptionalFloat(gkminNode, "pureSolventValue", CROP_ln_gamma_k_min);
        }
        if (cropNode.hasChild("ln_gamma_k_max")) {
            XML_Node& gkmaxNode = cropNode.child("ln_gamma_k_max");
            getOptionalFloat(gkmaxNode, "pureSolventValue", CROP_ln_gamma_k_max);
        }

        if (cropNode.hasChild("ln_gamma_o_min")) {
            XML_Node& gominNode = cropNode.child("ln_gamma_o_min");
            getOptionalFloat(gominNode, "pureSolventValue", CROP_ln_gamma_o_min);
        }

        if (cropNode.hasChild("ln_gamma_o_max")) {
            XML_Node& gomaxNode = cropNode.child("ln_gamma_o_max");
            getOptionalFloat(gomaxNode, "pureSolventValue", CROP_ln_gamma_o_max);
        }
    }
}

void HMWSoln::initThermo()
{
    MolalityVPSSTP::initThermo();
    for (int i = 0; i < 17; i++) {
        elambda[i] = 0.0;
        elambda1[i] = 0.0;
    }
    initLengths();
}

void HMWSoln::constructPhaseFile(std::string inputFile, std::string id_)
{
    warn_deprecated("HMWSoln::constructPhaseFile",
        "Use initThermoFile instead. To be removed after Cantera 2.3.");

    initThermoFile(inputFile, id_);
}

void HMWSoln::constructPhaseXML(XML_Node& phaseNode, std::string id_)
{
    warn_deprecated("HMWSoln::constructPhaseXML",
        "Use importPhase instead. To be removed after Cantera 2.3.");
    importPhase(phaseNode, this);
}

void HMWSoln::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id_) {
            throw CanteraError("HMWSoln::initThermoXML",
                               "phasenode and Id are incompatible");
        }
    }

    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("HMWSoln::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    // Possibly change the form of the standard concentrations
    if (thermoNode.hasChild("standardConc")) {
        XML_Node& scNode = thermoNode.child("standardConc");
        m_formGC = 2;
        string formString = scNode.attrib("model");
        if (formString != "") {
            if (ba::iequals(formString, "unity")) {
                m_formGC = 0;
                throw CanteraError("HMWSoln::initThermoXML",
                                   "standardConc = unity not done");
            } else if (ba::iequals(formString, "molar_volume")) {
                m_formGC = 1;
                throw CanteraError("HMWSoln::initThermoXML",
                                   "standardConc = molar_volume not done");
            } else if (ba::iequals(formString, "solvent_volume")) {
                m_formGC = 2;
            } else {
                throw CanteraError("HMWSoln::initThermoXML",
                                   "Unknown standardConc model: " + formString);
            }
        }
    }

    // Determine the form of the Pitzer model, We will use this information to
    // size arrays below.
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& scNode = thermoNode.child("activityCoefficients");
        string formString = scNode.attrib("model");
        if (formString != "") {
            if (ba::iequals(formString, "pitzer") || ba::iequals(formString, "default")) {
                m_formPitzer = PITZERFORM_BASE;
            } else if (ba::iequals(formString, "base")) {
                m_formPitzer = PITZERFORM_BASE;
            } else {
                throw CanteraError("HMWSoln::initThermoXML",
                                   "Unknown Pitzer ActivityCoeff model: "
                                   + formString);
            }
        }

        // Determine the form of the temperature dependence of the Pitzer
        // activity coefficient model.
        formString = scNode.attrib("TempModel");
        if (formString != "") {
            if (ba::iequals(formString, "constant") || ba::iequals(formString, "default")) {
                m_formPitzerTemp = PITZER_TEMP_CONSTANT;
            } else if (ba::iequals(formString, "linear")) {
                m_formPitzerTemp = PITZER_TEMP_LINEAR;
            } else if (ba::iequals(formString, "complex") || ba::iequals(formString, "complex1")) {
                m_formPitzerTemp = PITZER_TEMP_COMPLEX1;
            } else {
                throw CanteraError("HMWSoln::initThermoXML",
                                   "Unknown Pitzer ActivityCoeff Temp model: "
                                   + formString);
            }
        }

        // Determine the reference temperature of the Pitzer activity
        // coefficient model's temperature dependence formulation: defaults to
        // 25C
        formString = scNode.attrib("TempReference");
        if (formString != "") {
            m_TempPitzerRef = fpValueCheck(formString);
        } else {
            m_TempPitzerRef = 273.15 + 25;
        }
    }

    // Get the Name of the Solvent:
    //      <solvent> solventName </solvent>
    string solventName = "";
    if (thermoNode.hasChild("solvent")) {
        XML_Node& scNode = thermoNode.child("solvent");
        vector<string> nameSolventa;
        getStringArray(scNode, nameSolventa);
        if (nameSolventa.size() != 1) {
            throw CanteraError("HMWSoln::initThermoXML",
                               "badly formed solvent XML node");
        }
        solventName = nameSolventa[0];
    }

    // Initialize all of the lengths of arrays in the object
    // now that we know what species are in the phase.
    initLengths();

    // Reconcile the solvent name and index.
    for (size_t k = 0; k < m_kk; k++) {
        string sname = speciesName(k);
        if (solventName == sname) {
            setSolvent(k);
            if (k != 0) {
                throw CanteraError("HMWSoln::initThermoXML",
                                   "Solvent must be species 0 atm");
            }
            m_indexSolvent = k;
            break;
        }
    }
    if (m_indexSolvent == npos) {
        std::cout << "HMWSoln::initThermo: Solvent Name not found"
                  << std::endl;
        throw CanteraError("HMWSoln::initThermoXML",
                           "Solvent name not found");
    }
    if (m_indexSolvent != 0) {
        throw CanteraError("HMWSoln::initThermoXML",
                           "Solvent " + solventName +
                           " should be first species");
    }

    // Now go get the specification of the standard states for species in the
    // solution. This includes the molar volumes data blocks for incompressible
    // species.
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB =
        get_XML_NameID("speciesData", speciesList["datasrc"],
                       &phaseNode.root());
    const vector<string>&sss = speciesNames();

    for (size_t k = 0; k < m_kk; k++) {
        XML_Node* s = speciesDB->findByAttr("name", sss[k]);
        if (!s) {
            throw CanteraError("HMWSoln::initThermoXML",
                               "Species Data Base " + sss[k] + " not found");
        }
        XML_Node* ss = s->findByName("standardState");
        if (!ss) {
            throw CanteraError("HMWSoln::initThermoXML",
                               "Species " + sss[k] +
                               " standardState XML block not found");
        }
        string modelString = ba::to_lower_copy(ss->attrib("model"));
        if (modelString == "") {
            throw CanteraError("HMWSoln::initThermoXML",
                               "Species " + sss[k] +
                               " standardState XML block model attribute not found");
        }
        if (k == 0) {
            if (modelString == "wateriapws" || modelString == "real_water" ||
                    modelString == "waterpdss") {

                // Store a local pointer to the water standard state model.
                // We've hardcoded it to a PDSS_Water model, so this is ok.
                m_waterSS = dynamic_cast<PDSS_Water*>(providePDSS(0));
                if (!m_waterSS) {
                    throw CanteraError("HMWSoln::initThermoXML",
                                       "Dynamic cast to PDSS_Water failed");
                }

                // Fill in the molar volume of water (m3/kmol) at standard
                // conditions to fill in the m_speciesSize entry with something
                // reasonable.
                m_waterSS->setState_TP(300., OneAtm);
                double dens = m_waterSS->density();
                double mw = m_waterSS->molecularWeight();
                m_speciesSize[0] = mw / dens;
            } else {
                m_waterSS = providePDSS(0);
                m_waterSS->setState_TP(300., OneAtm);
                double dens = m_waterSS->density();
                double mw = m_waterSS->molecularWeight();
                m_speciesSize[0] = mw / dens;
            }
        } else {
            if (modelString != "constant_incompressible" && modelString != "hkft") {
                throw CanteraError("HMWSoln::initThermoXML",
                                   "Solute SS Model \"" + modelString +
                                   "\" is not known");
            }
            if (modelString ==  "constant_incompressible") {
                m_speciesSize[k] = getFloat(*ss, "molarVolume", "toSI");
            }
            // HKM Note, have to fill up m_speciesSize[] for HKFT species
        }
    }

    // Initialize the water property calculator. It will share the internal eos
    // water calculator.
    m_waterProps.reset(new WaterProps(&dynamic_cast<PDSS_Water&>(*m_waterSS)));

    // Fill in parameters for the calculation of the stoichiometric Ionic
    // Strength. The default is that stoich charge is the same as the regular
    // charge.
    for (size_t k = 0; k < m_kk; k++) {
        m_speciesCharge_Stoich[k] = charge(k);
    }

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;

        // Look for parameters for A_Debye
        if (acNode.hasChild("A_Debye")) {
            XML_Node& ADebye = acNode.child("A_Debye");
            m_form_A_Debye = A_DEBYE_CONST;
            string stemp = "model";
            if (ADebye.hasAttrib("model")) {
                if (ba::iequals(ADebye.attrib("model"), "water")) {
                    m_form_A_Debye = A_DEBYE_WATER;
                }
            }
            if (m_form_A_Debye == A_DEBYE_CONST) {
                m_A_Debye = getFloat(acNode, "A_Debye");
            }
        }

        // Look for Parameters for the Maximum Ionic Strength
        if (acNode.hasChild("maxIonicStrength")) {
            m_maxIionicStrength = getFloat(acNode, "maxIonicStrength");
        }

        // Look for parameters for the Ionic radius
        if (acNode.hasChild("ionicRadius")) {
            XML_Node& irNode = acNode.child("ionicRadius");
            double Afactor = 1.0;
            if (irNode.hasAttrib("units")) {
                string Aunits = irNode.attrib("units");
                Afactor = toSI(Aunits);
            }

            if (irNode.hasAttrib("default")) {
                string ads = irNode.attrib("default");
                double ad = fpValue(ads);
                for (size_t k = 0; k < m_kk; k++) {
                    m_Aionic[k] = ad * Afactor;
                }
            }
        }

        // First look at the species database. Look for the subelement
        // "stoichIsMods" in each of the species SS databases.
        std::vector<const XML_Node*> xspecies = speciesData();
        for (size_t k = 0; k < m_kk; k++) {
            size_t jmap = npos;
            string kname = speciesName(k);
            for (size_t j = 0; j < xspecies.size(); j++) {
                const XML_Node& sp = *xspecies[j];
                string jname = sp["name"];
                if (jname == kname) {
                    jmap = j;
                    break;
                }
            }
            if (jmap != npos) {
                const XML_Node& sp = *xspecies[jmap];
                getOptionalFloat(sp, "stoichIsMods", m_speciesCharge_Stoich[k]);
            }
        }

        // Now look at the activity coefficient database
        if (acNodePtr && acNodePtr->hasChild("stoichIsMods")) {
            XML_Node& sIsNode = acNodePtr->child("stoichIsMods");
            map<string, string> msIs;
            getMap(sIsNode, msIs);
            for (const auto& b : msIs) {
                size_t kk = speciesIndex(b.first);
                if (kk != npos) {
                    double val = fpValue(b.second);
                    m_speciesCharge_Stoich[kk] = val;
                }
            }
        }

        // Loop through the children getting multiple instances of parameters
        if (acNodePtr) {
            for (size_t i = 0; i < acNodePtr->nChildren(); i++) {
                XML_Node& xmlACChild = acNodePtr->child(i);
                string nodeName = xmlACChild.name();

                // Process a binary salt field, or any of the other XML fields
                // that make up the Pitzer Database. Entries will be ignored
                // if any of the species in the entry isn't in the solution.
                if (ba::iequals(nodeName, "binarysaltparameters")) {
                    readXMLBinarySalt(xmlACChild);
                } else if (ba::iequals(nodeName, "thetaanion")) {
                    readXMLThetaAnion(xmlACChild);
                } else if (ba::iequals(nodeName, "thetacation")) {
                    readXMLThetaCation(xmlACChild);
                } else if (ba::iequals(nodeName, "psicommonanion")) {
                    readXMLPsiCommonAnion(xmlACChild);
                } else if (ba::iequals(nodeName, "psicommoncation")) {
                    readXMLPsiCommonCation(xmlACChild);
                } else if (ba::iequals(nodeName, "lambdaneutral")) {
                    readXMLLambdaNeutral(xmlACChild);
                } else if (ba::iequals(nodeName, "zetacation")) {
                    readXMLZetaCation(xmlACChild);
                }
            }
        }

        // Go look up the optional Cropping parameters
        readXMLCroppingCoefficients(acNode);
    }

    // Fill in the vector specifying the electrolyte species type
    //
    // First fill in default values. Everything is either a charge species, a
    // nonpolar neutral, or the solvent.
    for (size_t k = 0; k < m_kk; k++) {
        if (fabs(charge(k)) > 0.0001) {
            m_electrolyteSpeciesType[k] = cEST_chargedSpecies;
            if (fabs(m_speciesCharge_Stoich[k] - charge(k)) > 0.0001) {
                m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
            }
        } else if (fabs(m_speciesCharge_Stoich[k]) > 0.0001) {
            m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
        } else {
            m_electrolyteSpeciesType[k] = cEST_nonpolarNeutral;
        }
    }
    m_electrolyteSpeciesType[m_indexSolvent] = cEST_solvent;

    // First look at the species database. Look for the subelement
    // "stoichIsMods" in each of the species SS databases.
    std::vector<const XML_Node*> xspecies = speciesData();
    for (size_t k = 0; k < m_kk; k++) {
        const XML_Node* spPtr = xspecies[k];
        if (spPtr && spPtr->hasChild("electrolyteSpeciesType")) {
            string est = getChildValue(*spPtr, "electrolyteSpeciesType");
            if ((m_electrolyteSpeciesType[k] = interp_est(est)) == -1) {
                throw CanteraError("HMWSoln::initThermoXML",
                                   "Bad electrolyte type: " + est);
            }
        }
    }

    // Then look at the phase thermo specification
    if (acNodePtr && acNodePtr->hasChild("electrolyteSpeciesType")) {
        XML_Node& ESTNode = acNodePtr->child("electrolyteSpeciesType");
        map<string, string> msEST;
        getMap(ESTNode, msEST);
        for (const auto& b : msEST) {
            size_t kk = speciesIndex(b.first);
            if (kk != npos) {
                string est = b.second;
                if ((m_electrolyteSpeciesType[kk] = interp_est(est))  == -1) {
                    throw CanteraError("HMWSoln::initThermoXML",
                                       "Bad electrolyte type: " + est);
                }
            }
        }
    }

    IMS_typeCutoff_ = 2;
    if (IMS_typeCutoff_ == 2) {
        calcIMSCutoffParams_();
    }
    calcMCCutoffParams_();
    setMoleFSolventMin(1.0E-5);

    MolalityVPSSTP::initThermoXML(phaseNode, id_);

    // Lastly calculate the charge balance and then add stuff until the charges
    // compensate
    vector_fp mf(m_kk, 0.0);
    getMoleFractions(mf.data());
    bool notDone = true;

    while (notDone) {
        double sum = 0.0;
        size_t kMaxC = npos;
        double MaxC = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            sum += mf[k] * charge(k);
            if (fabs(mf[k] * charge(k)) > MaxC) {
                kMaxC = k;
            }
        }
        size_t kHp = speciesIndex("H+");
        size_t kOHm = speciesIndex("OH-");

        if (fabs(sum) > 1.0E-30) {
            if (kHp != npos) {
                if (mf[kHp] > sum * 1.1) {
                    mf[kHp] -= sum;
                    mf[0] += sum;
                    notDone = false;
                } else {
                    if (sum > 0.0) {
                        mf[kHp] *= 0.5;
                        mf[0] += mf[kHp];
                        sum -= mf[kHp];
                    }
                }
            }
            if (notDone) {
                if (kOHm != npos) {
                    if (mf[kOHm] > -sum * 1.1) {
                        mf[kOHm] += sum;
                        mf[0] -= sum;
                        notDone = false;
                    } else {
                        if (sum < 0.0) {
                            mf[kOHm] *= 0.5;
                            mf[0] += mf[kOHm];
                            sum += mf[kOHm];
                        }
                    }
                }
                if (notDone && kMaxC != npos) {
                    if (mf[kMaxC] > (1.1 * sum / charge(kMaxC))) {
                        mf[kMaxC] -= sum / charge(kMaxC);
                        mf[0] += sum / charge(kMaxC);
                    } else {
                        mf[kMaxC] *= 0.5;
                        mf[0] += mf[kMaxC];
                        notDone = true;
                    }
                }
            }
            setMoleFractions(mf.data());
        } else {
            notDone = false;
        }
    }
}

void HMWSoln::calcIMSCutoffParams_()
{
    IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_);
    IMS_efCut_ = 0.0;
    bool converged = false;
    double oldV = 0.0;
    for (int its = 0; its < 100 && !converged; its++) {
        oldV = IMS_efCut_;
        IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_) -IMS_efCut_;
        IMS_bfCut_ = IMS_afCut_ / IMS_cCut_ + IMS_slopefCut_ - 1.0;
        IMS_dfCut_ = ((- IMS_afCut_/IMS_cCut_ + IMS_bfCut_ - IMS_bfCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_afCut_ + IMS_X_o_cutoff_*(IMS_bfCut_ + IMS_dfCut_ *IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_efCut_ = - eterm * tmp;
        if (fabs(IMS_efCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("HMWSoln::calcIMSCutoffParams_()",
                           " failed to converge on the f polynomial");
    }
    converged = false;
    double f_0 = IMS_afCut_ + IMS_efCut_;
    double f_prime_0 = 1.0 - IMS_afCut_ / IMS_cCut_ + IMS_bfCut_;
    IMS_egCut_ = 0.0;
    for (int its = 0; its < 100 && !converged; its++) {
        oldV = IMS_egCut_;
        double lng_0 = -log(IMS_gamma_o_min_) - f_prime_0 / f_0;
        IMS_agCut_ = exp(lng_0) - IMS_egCut_;
        IMS_bgCut_ = IMS_agCut_ / IMS_cCut_ + IMS_slopegCut_ - 1.0;
        IMS_dgCut_ = ((- IMS_agCut_/IMS_cCut_ + IMS_bgCut_ - IMS_bgCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_agCut_ + IMS_X_o_cutoff_*(IMS_bgCut_ + IMS_dgCut_ *IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_egCut_ = - eterm * tmp;
        if (fabs(IMS_egCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("HMWSoln::calcIMSCutoffParams_()",
                           " failed to converge on the g polynomial");
    }
}

void HMWSoln::calcMCCutoffParams_()
{
    MC_X_o_min_ = 0.35;
    MC_X_o_cutoff_ = 0.6;
    MC_slopepCut_ = 0.02;
    MC_cpCut_ = 0.25;

    // Initial starting values
    MC_apCut_ = MC_X_o_min_;
    MC_epCut_ = 0.0;
    bool converged = false;
    double oldV = 0.0;
    double damp = 0.5;
    for (int its = 0; its < 500 && !converged; its++) {
        oldV = MC_epCut_;
        MC_apCut_ = damp *(MC_X_o_min_ - MC_epCut_) + (1-damp) * MC_apCut_;
        double MC_bpCutNew = MC_apCut_ / MC_cpCut_ + MC_slopepCut_ - 1.0;
        MC_bpCut_ = damp * MC_bpCutNew + (1-damp) * MC_bpCut_;
        double MC_dpCutNew = ((- MC_apCut_/MC_cpCut_ + MC_bpCut_ - MC_bpCut_ * MC_X_o_cutoff_/MC_cpCut_)
                              /
                              (MC_X_o_cutoff_ * MC_X_o_cutoff_/MC_cpCut_ - 2.0 * MC_X_o_cutoff_));
        MC_dpCut_ = damp * MC_dpCutNew + (1-damp) * MC_dpCut_;
        double tmp = MC_apCut_ + MC_X_o_cutoff_*(MC_bpCut_ + MC_dpCut_ * MC_X_o_cutoff_);
        double eterm = std::exp(- MC_X_o_cutoff_ / MC_cpCut_);
        MC_epCut_ = - eterm * tmp;
        double diff = MC_epCut_ - oldV;
        if (fabs(diff) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError("HMWSoln::calcMCCutoffParams_()",
                           " failed to converge on the p polynomial");
    }
}

}
