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

static void check_nParams(const std::string& method, size_t nParams,
                          size_t m_formPitzerTemp)
{
    if (m_formPitzerTemp == PITZER_TEMP_CONSTANT && nParams != 1) {
        throw CanteraError(method, "'constant' temperature model requires one"
            " coefficient for each of parameter, but {} were given", nParams);
    } else if (m_formPitzerTemp == PITZER_TEMP_LINEAR && nParams != 2) {
        throw CanteraError(method, "'linear' temperature model requires two"
            " coefficients for each parameter, but {} were given", nParams);
    }
    if (m_formPitzerTemp == PITZER_TEMP_COMPLEX1 && nParams != 5) {
        throw CanteraError(method, "'complex' temperature model requires five"
            " coefficients for each parameter, but {} were given", nParams);
    }
}

void HMWSoln::readXMLBinarySalt(XML_Node& BinSalt)
{
    if (BinSalt.name() != "binarySaltParameters") {
        throw CanteraError("HMWSoln::readXMLBinarySalt",
                           "Incorrect name for processing this routine: " + BinSalt.name());
    }

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
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos) {
        return;
    }

    vector_fp beta0, beta1, beta2, Cphi;
    getFloatArray(BinSalt, beta0, false, "", "beta0");
    getFloatArray(BinSalt, beta1, false, "", "beta1");
    getFloatArray(BinSalt, beta2, false, "", "beta2");
    getFloatArray(BinSalt, Cphi, false, "", "Cphi");
    if (beta0.size() != beta1.size() || beta0.size() != beta2.size() ||
        beta0.size() != Cphi.size()) {
        throw CanteraError("HMWSoln::readXMLBinarySalt", "Inconsistent"
            " array sizes ({}, {}, {}, {})", beta0.size(), beta1.size(),
            beta2.size(), Cphi.size());
    }
    if (beta0.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        beta0.resize(5, 0.0);
        beta1.resize(5, 0.0);
        beta2.resize(5, 0.0);
        Cphi.resize(5, 0.0);
    }
    double alpha1 = getFloat(BinSalt, "Alpha1");
    double alpha2 = 0.0;
    getOptionalFloat(BinSalt, "Alpha2", alpha2);
    setBinarySalt(iName, jName, beta0.size(), beta0.data(), beta1.data(),
        beta2.data(), Cphi.data(), alpha1, alpha2);
}

void HMWSoln::setBinarySalt(const std::string& sp1, const std::string& sp2,
    size_t nParams, double* beta0, double* beta1, double* beta2,
    double* Cphi, double alpha1, double alpha2)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setBinarySalt", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setBinarySalt", "Species '{}' not found", sp2);
    }
    if (charge(k1) < 0 && charge(k2) > 0) {
        std::swap(k1, k2);
    } else if (charge(k1) * charge(k2) >= 0) {
        throw CanteraError("HMWSoln::setBinarySalt", "Species '{}' and '{}' "
            "do not have opposite charges ({}, {})", sp1, sp2,
            charge(k1), charge(k2));
    }
    check_nParams("HMWSoln::setBinarySalt", nParams, m_formPitzerTemp);

    size_t c = m_CounterIJ[k1 * m_kk + k2];
    m_Beta0MX_ij[c] = beta0[0];
    m_Beta1MX_ij[c] = beta1[0];
    m_Beta2MX_ij[c] = beta2[0];
    m_CphiMX_ij[c] = Cphi[0];
    for (size_t n = 0; n < nParams; n++) {
        m_Beta0MX_ij_coeff(n, c) = beta0[n];
        m_Beta1MX_ij_coeff(n, c) = beta1[n];
        m_Beta2MX_ij_coeff(n, c) = beta2[n];
        m_CphiMX_ij_coeff(n, c) = Cphi[n];
    }
    m_Alpha1MX_ij[c] = alpha1;
    m_Alpha2MX_ij[c] = alpha2;
}


void HMWSoln::readXMLTheta(XML_Node& node)
{
    string ispName, jspName;
    if (node.name() == "thetaAnion") {
        ispName = node.attrib("anion1");
        if (ispName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no anion1 attrib");
        }
        jspName = node.attrib("anion2");
        if (jspName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no anion2 attrib");
        }
    } else if (node.name() == "thetaCation") {
        ispName = node.attrib("cation1");
        if (ispName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no cation1 attrib");
        }
        jspName = node.attrib("cation2");
        if (jspName == "") {
            throw CanteraError("HMWSoln::readXMLTheta", "no cation2 attrib");
        }
    } else {
        throw CanteraError("HMWSoln::readXMLTheta",
                "Incorrect name for processing this routine: " + node.name());
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(ispName) == npos || speciesIndex(jspName) == npos) {
        return;
    }

    vector_fp theta;
    getFloatArray(node, theta, false, "", "theta");
    if (theta.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        theta.resize(5, 0.0);
    }
    setTheta(ispName, jspName, theta.size(), theta.data());
}

void HMWSoln::setTheta(const std::string& sp1, const std::string& sp2,
        size_t nParams, double* theta)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setTheta", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setTheta", "Species '{}' not found", sp2);
    }
    if (charge(k1) * charge(k2) <= 0) {
        throw CanteraError("HMWSoln::setTheta", "Species '{}' and '{}' "
            "should both have the same (non-zero) charge ({}, {})", sp1, sp2,
            charge(k1), charge(k2));
    }
    check_nParams("HMWSoln::setTheta", nParams, m_formPitzerTemp);
    size_t c = m_CounterIJ[k1 * m_kk + k2];
    m_Theta_ij[c] = theta[0];
    for (size_t n = 0; n < nParams; n++) {
        m_Theta_ij_coeff(n, c) = theta[n];
    }
}

void HMWSoln::readXMLPsi(XML_Node& node)
{
    string iName, jName, kName;
    if (node.name() == "psiCommonCation") {
        kName = node.attrib("cation");
        if (kName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no cation attrib");
        }
        iName = node.attrib("anion1");
        if (iName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no anion1 attrib");
        }
        jName = node.attrib("anion2");
        if (jName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no anion2 attrib");
        }
    } else if (node.name() == "psiCommonAnion") {
        kName = node.attrib("anion");
        if (kName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no anion attrib");
        }
        iName = node.attrib("cation1");
        if (iName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no cation1 attrib");
        }
        jName = node.attrib("cation2");
        if (jName == "") {
            throw CanteraError("HMWSoln::readXMLPsi", "no cation2 attrib");
        }
    } else {
        throw CanteraError("HMWSoln::readXMLPsi",
                   "Incorrect name for processing this routine: " + node.name());
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos ||
        speciesIndex(kName) == npos) {
        return;
    }

    vector_fp psi;
    getFloatArray(node, psi, false, "", "psi");
    if (psi.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        psi.resize(5, 0.0);
    }
    setPsi(iName, jName, kName, psi.size(), psi.data());
}

void HMWSoln::setPsi(const std::string& sp1, const std::string& sp2,
        const std::string& sp3, size_t nParams, double* psi)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    size_t k3 = speciesIndex(sp3);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setPsi", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setPsi", "Species '{}' not found", sp2);
    } else if (k3 == npos) {
        throw CanteraError("HMWSoln::setPsi", "Species '{}' not found", sp3);
    }

    if (!charge(k1) || !charge(k2) || !charge(k3) ||
        std::abs(sign(charge(k1) + sign(charge(k2)) + sign(charge(k3)))) != 1) {
        throw CanteraError("HMWSoln::setPsi", "All species must be ions and"
            " must include at least one cation and one anion, but given species"
            " (charges) were: {} ({}), {} ({}), and {} ({}).",
            sp1, charge(k1), sp2, charge(k2), sp3, charge(k3));
    }
    check_nParams("HMWSoln::setPsi", nParams, m_formPitzerTemp);
    auto cc = {k1*m_kk*m_kk + k2*m_kk + k3,
               k1*m_kk*m_kk + k3*m_kk + k2,
               k2*m_kk*m_kk + k1*m_kk + k3,
               k2*m_kk*m_kk + k3*m_kk + k1,
               k3*m_kk*m_kk + k2*m_kk + k1,
               k3*m_kk*m_kk + k1*m_kk + k2};
    for (auto c : cc) {
        for (size_t n = 0; n < nParams; n++) {
            m_Psi_ijk_coeff(n, c) = psi[n];
        }
        m_Psi_ijk[c] = psi[0];
    }
}

void HMWSoln::readXMLLambdaNeutral(XML_Node& node)
{
    vector_fp vParams;
    if (node.name() != "lambdaNeutral") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral",
                           "Incorrect name for processing this routine: " + node.name());
    }
    string iName = node.attrib("species1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no species1 attrib");
    }
    string jName = node.attrib("species2");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLLambdaNeutral", "no species2 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos) {
        return;
    }

    vector_fp lambda;
    getFloatArray(node, lambda, false, "", "lambda");
    if (lambda.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        lambda.resize(5, 0.0);
    }
    setLambda(iName, jName, lambda.size(), lambda.data());
}

void HMWSoln::setLambda(const std::string& sp1, const std::string& sp2,
        size_t nParams, double* lambda)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setLambda", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setLambda", "Species '{}' not found", sp2);
    }

    if (charge(k1) != 0 && charge(k2) != 0) {
        throw CanteraError("HMWSoln::setLambda", "Expected at least one neutral"
            " species, but given species (charges) were: {} ({}) and {} ({}).",
            sp1, charge(k1), sp2, charge(k2));
    }
    if (charge(k1) != 0) {
        std::swap(k1, k2);
    }
    check_nParams("HMWSoln::setLambda", nParams, m_formPitzerTemp);
    size_t c = k1*m_kk + k2;
    for (size_t n = 0; n < nParams; n++) {
        m_Lambda_nj_coeff(n, c) = lambda[n];
    }
    m_Lambda_nj(k1, k2) = lambda[0];
}

void HMWSoln::readXMLMunnnNeutral(XML_Node& node)
{
    if (node.name() != "MunnnNeutral") {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral",
                           "Incorrect name for processing this routine: " + node.name());
    }
    string iName = node.attrib("species1");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLMunnnNeutral", "no species1 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos) {
        return;
    }

    vector_fp munnn;
    getFloatArray(node, munnn, false, "", "munnn");
    if (munnn.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        munnn.resize(5, 0.0);
    }
    setMunnn(iName, munnn.size(), munnn.data());
}

void HMWSoln::setMunnn(const std::string& sp, size_t nParams, double* munnn)
{
    size_t k = speciesIndex(sp);
    if (k == npos) {
        throw CanteraError("HMWSoln::setMunnn", "Species '{}' not found", sp);
    }

    if (charge(k) != 0) {
        throw CanteraError("HMWSoln::setMunnn", "Expected a neutral species,"
                " got {} ({}).", sp, charge(k));
    }
    check_nParams("HMWSoln::setMunnn", nParams, m_formPitzerTemp);
    for (size_t n = 0; n < nParams; n++) {
        m_Mu_nnn_coeff(n, k) = munnn[n];
    }
    m_Mu_nnn[k] = munnn[0];

}

void HMWSoln::readXMLZetaCation(const XML_Node& node)
{
    if (node.name() != "zetaCation") {
        throw CanteraError("HMWSoln::readXMLZetaCation",
                           "Incorrect name for processing this routine: " + node.name());
    }

    string iName = node.attrib("neutral");
    if (iName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no neutral attrib");
    }
    string jName = node.attrib("cation1");
    if (jName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no cation1 attrib");
    }
    string kName = node.attrib("anion1");
    if (kName == "") {
        throw CanteraError("HMWSoln::readXMLZetaCation", "no anion1 attrib");
    }

    // Find the index of the species in the current phase. It's not an error to
    // not find the species
    if (speciesIndex(iName) == npos || speciesIndex(jName) == npos ||
        speciesIndex(kName) == npos) {
        return;
    }

    vector_fp zeta;
    getFloatArray(node, zeta, false, "", "zeta");
    if (zeta.size() == 1 && m_formPitzerTemp == PITZER_TEMP_COMPLEX1) {
        zeta.resize(5, 0.0);
    }
    setZeta(iName, jName, kName, zeta.size(), zeta.data());
}

void HMWSoln::setZeta(const std::string& sp1, const std::string& sp2,
        const std::string& sp3, size_t nParams, double* psi)
{
    size_t k1 = speciesIndex(sp1);
    size_t k2 = speciesIndex(sp2);
    size_t k3 = speciesIndex(sp3);
    if (k1 == npos) {
        throw CanteraError("HMWSoln::setZeta", "Species '{}' not found", sp1);
    } else if (k2 == npos) {
        throw CanteraError("HMWSoln::setZeta", "Species '{}' not found", sp2);
    } else if (k3 == npos) {
        throw CanteraError("HMWSoln::setZeta", "Species '{}' not found", sp3);
    }

    if (charge(k1)*charge(k2)*charge(k3) != 0 ||
        sign(charge(k1)) + sign(charge(k2)) + sign(charge(k3)) != 0) {
        throw CanteraError("HMWSoln::setZeta", "Requires one neutral species, "
            "one cation, and one anion, but given species (charges) were: "
            "{} ({}), {} ({}), and {} ({}).",
            sp1, charge(k1), sp2, charge(k2), sp3, charge(k3));
    }

    //! Make k1 the neutral species
    if (charge(k2) == 0) {
        std::swap(k1, k2);
    } else if (charge(k3) == 0) {
        std::swap(k1, k3);
    }

    // Make k2 the cation
    if (charge(k3) > 0) {
        std::swap(k2, k3);
    }

    check_nParams("HMWSoln::setZeta", nParams, m_formPitzerTemp);
    // In contrast to setPsi, there are no duplicate entries
    size_t c = k1 * m_kk *m_kk + k2 * m_kk + k3;
    for (size_t n = 0; n < nParams; n++) {
        m_Psi_ijk_coeff(n, c) = psi[n];
    }
    m_Psi_ijk[c] = psi[0];
}

void HMWSoln::setPitzerTempModel(const std::string& model)
{
    if (ba::iequals(model, "constant") || ba::iequals(model, "default")) {
        m_formPitzerTemp = PITZER_TEMP_CONSTANT;
    } else if (ba::iequals(model, "linear")) {
        m_formPitzerTemp = PITZER_TEMP_LINEAR;
    } else if (ba::iequals(model, "complex") || ba::iequals(model, "complex1")) {
        m_formPitzerTemp = PITZER_TEMP_COMPLEX1;
    } else {
        throw CanteraError("HMWSoln::setPitzerTempModel",
                           "Unknown Pitzer ActivityCoeff Temp model: {}", model);
    }
}

void HMWSoln::setA_Debye(double A)
{
    if (A < 0) {
        m_form_A_Debye = A_DEBYE_WATER;
    } else {
        m_form_A_Debye = A_DEBYE_CONST;
        m_A_Debye = A;
    }
}

void HMWSoln::setCroppingCoefficients(double ln_gamma_k_min,
    double ln_gamma_k_max, double ln_gamma_o_min, double ln_gamma_o_max)
{
        CROP_ln_gamma_k_min = ln_gamma_k_min;
        CROP_ln_gamma_k_max = ln_gamma_k_max;
        CROP_ln_gamma_o_min = ln_gamma_o_min;
        CROP_ln_gamma_o_max = ln_gamma_o_max;
}

void HMWSoln::initThermo()
{
    MolalityVPSSTP::initThermo();
    initLengths();

    for (int i = 0; i < 17; i++) {
        elambda[i] = 0.0;
        elambda1[i] = 0.0;
    }
    for (size_t k = 0; k < nSpecies(); k++) {
        m_speciesSize[k] = providePDSS(k)->molarVolume();
    }

    // Store a local pointer to the water standard state model.
    m_waterSS = providePDSS(0);

    // Initialize the water property calculator. It will share the internal eos
    // water calculator.
    m_waterProps.reset(new WaterProps(dynamic_cast<PDSS_Water*>(m_waterSS)));

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

    calcIMSCutoffParams_();
    calcMCCutoffParams_();
    setMoleFSolventMin(1.0E-5);
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

    // Determine the form of the Pitzer model, We will use this information to
    // size arrays below.
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& scNode = thermoNode.child("activityCoefficients");

        // Determine the form of the temperature dependence of the Pitzer
        // activity coefficient model.
        string formString = scNode.attrib("TempModel");
        if (formString != "") {
            setPitzerTempModel(formString);
        }

        // Determine the reference temperature of the Pitzer activity
        // coefficient model's temperature dependence formulation: defaults to
        // 25C
        formString = scNode.attrib("TempReference");
        if (formString != "") {
            setPitzerRefTemperature(fpValueCheck(formString));
        }
    }

    // Initialize all of the lengths of arrays in the object
    // now that we know what species are in the phase.
    initLengths();

    // Go get all of the coefficients and factors in the activityCoefficients
    // XML block
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");

        // Look for parameters for A_Debye
        if (acNode.hasChild("A_Debye")) {
            XML_Node& ADebye = acNode.child("A_Debye");
            if (ba::iequals(ADebye["model"], "water")) {
                setA_Debye(-1);
            } else {
                setA_Debye(getFloat(acNode, "A_Debye"));
            }
        }

        // Look for Parameters for the Maximum Ionic Strength
        if (acNode.hasChild("maxIonicStrength")) {
            setMaxIonicStrength(getFloat(acNode, "maxIonicStrength"));
        }

        for (const auto& xmlACChild : acNode.children()) {
            string nodeName = xmlACChild->name();

            // Process any of the XML fields that make up the Pitzer Database.
            // Entries will be ignored if any of the species in the entry aren't
            // in the solution.
            if (ba::iequals(nodeName, "binarysaltparameters")) {
                readXMLBinarySalt(*xmlACChild);
            } else if (ba::iequals(nodeName, "thetaanion")) {
                readXMLTheta(*xmlACChild);
            } else if (ba::iequals(nodeName, "thetacation")) {
                readXMLTheta(*xmlACChild);
            } else if (ba::iequals(nodeName, "psicommonanion")) {
                readXMLPsi(*xmlACChild);
            } else if (ba::iequals(nodeName, "psicommoncation")) {
                readXMLPsi(*xmlACChild);
            } else if (ba::iequals(nodeName, "lambdaneutral")) {
                readXMLLambdaNeutral(*xmlACChild);
            } else if (ba::iequals(nodeName, "zetacation")) {
                readXMLZetaCation(*xmlACChild);
            }
        }

        // Go look up the optional Cropping parameters
        if (acNode.hasChild("croppingCoefficients")) {
            XML_Node& cropNode = acNode.child("croppingCoefficients");
            setCroppingCoefficients(
                getFloat(cropNode.child("ln_gamma_k_min"), "pureSolventValue"),
                getFloat(cropNode.child("ln_gamma_k_max"), "pureSolventValue"),
                getFloat(cropNode.child("ln_gamma_o_min"), "pureSolventValue"),
                getFloat(cropNode.child("ln_gamma_o_max"), "pureSolventValue"));
        }
    }

    MolalityVPSSTP::initThermoXML(phaseNode, id_);

}

void HMWSoln::calcIMSCutoffParams_()
{
    double IMS_gamma_o_min_ = 1.0E-5; // value at the zero solvent point
    double IMS_gamma_k_min_ = 10.0; // minimum at the zero solvent point
    double IMS_slopefCut_ = 0.6; // slope of the f function at the zero solvent point

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
    double MC_X_o_min_ = 0.35; // value at the zero solvent point
    MC_X_o_cutoff_ = 0.6;
    double MC_slopepCut_ = 0.02; // slope of the p function at the zero solvent point
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
