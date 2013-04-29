/**
 *  @file MargulesVPSSTP.cpp
 *   Definitions for ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations related to Margules
 *   expansions (see \ref thermoprops
 *    and class \link Cantera::MargulesVPSSTP MargulesVPSSTP\endlink).
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

#include <iomanip>
#include <fstream>

using namespace std;

namespace Cantera
{
MargulesVPSSTP::MargulesVPSSTP() :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
}

MargulesVPSSTP::MargulesVPSSTP(const std::string& inputFile, const std::string& id_) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    initThermoFile(inputFile, id_);
}

MargulesVPSSTP::MargulesVPSSTP(XML_Node& phaseRoot, const std::string& id_) :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{
    importPhase(*findXMLPhase(&phaseRoot, id_), this);
}

MargulesVPSSTP::MargulesVPSSTP(const MargulesVPSSTP& b) :
    GibbsExcessVPSSTP()
{
    MargulesVPSSTP::operator=(b);
}

MargulesVPSSTP& MargulesVPSSTP::
operator=(const MargulesVPSSTP& b)
{
    if (&b == this) {
        return *this;
    }

    GibbsExcessVPSSTP::operator=(b);

    numBinaryInteractions_      = b.numBinaryInteractions_ ;
    m_HE_b_ij                   = b.m_HE_b_ij;
    m_HE_c_ij                   = b.m_HE_c_ij;
    m_HE_d_ij                   = b.m_HE_d_ij;
    m_SE_b_ij                   = b.m_SE_b_ij;
    m_SE_c_ij                   = b.m_SE_c_ij;
    m_SE_d_ij                   = b.m_SE_d_ij;
    m_VHE_b_ij                  = b.m_VHE_b_ij;
    m_VHE_c_ij                  = b.m_VHE_c_ij;
    m_VHE_d_ij                  = b.m_VHE_d_ij;
    m_VSE_b_ij                  = b.m_VSE_b_ij;
    m_VSE_c_ij                  = b.m_VSE_c_ij;
    m_VSE_d_ij                  = b.m_VSE_d_ij;
    m_pSpecies_A_ij             = b.m_pSpecies_A_ij;
    m_pSpecies_B_ij             = b.m_pSpecies_B_ij;
    formMargules_               = b.formMargules_;
    formTempModel_              = b.formTempModel_;

    return *this;
}

ThermoPhase*
MargulesVPSSTP::duplMyselfAsThermoPhase() const
{
    return new MargulesVPSSTP(*this);
}

MargulesVPSSTP::MargulesVPSSTP(int testProb)  :
    GibbsExcessVPSSTP(),
    numBinaryInteractions_(0),
    formMargules_(0),
    formTempModel_(0)
{

    initThermoFile("LiKCl_liquid.xml", "");


    numBinaryInteractions_ = 1;

    m_HE_b_ij.resize(1);
    m_HE_c_ij.resize(1);
    m_HE_d_ij.resize(1);

    m_SE_b_ij.resize(1);
    m_SE_c_ij.resize(1);
    m_SE_d_ij.resize(1);

    m_VHE_b_ij.resize(1);
    m_VHE_c_ij.resize(1);
    m_VHE_d_ij.resize(1);

    m_VSE_b_ij.resize(1);
    m_VSE_c_ij.resize(1);
    m_VSE_d_ij.resize(1);

    m_pSpecies_A_ij.resize(1);
    m_pSpecies_B_ij.resize(1);



    m_HE_b_ij[0] = -17570E3;
    m_HE_c_ij[0] = -377.0E3;
    m_HE_d_ij[0] = 0.0;

    m_SE_b_ij[0] = -7.627E3;
    m_SE_c_ij[0] =  4.958E3;
    m_SE_d_ij[0] =  0.0;


    size_t iLiCl = speciesIndex("LiCl(L)");
    if (iLiCl == npos) {
        throw CanteraError("MargulesVPSSTP test1 constructor",
                           "Unable to find LiCl(L)");
    }
    m_pSpecies_B_ij[0] = iLiCl;


    size_t iKCl = speciesIndex("KCl(L)");
    if (iKCl == npos) {
        throw CanteraError("MargulesVPSSTP test1 constructor",
                           "Unable to find KCl(L)");
    }
    m_pSpecies_A_ij[0] = iKCl;
}

/*
 *  -------------- Utilities -------------------------------
 */

int MargulesVPSSTP::eosType() const
{
    return 0;
}

/*
 * - Activities, Standard States, Activity Concentrations -----------
 */

void MargulesVPSSTP::getLnActivityCoefficients(doublereal* lnac) const
{
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();

    /*
     * take the exp of the internally stored coefficients.
     */
    for (size_t k = 0; k < m_kk; k++) {
        lnac[k] = lnActCoeff_Scaled_[k];
    }
}

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */

void MargulesVPSSTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}

void MargulesVPSSTP::getChemPotentials(doublereal* mu) const
{
    doublereal xx;
    /*
     * First get the standard chemical potentials in
     * molar form.
     *  -> this requires updates of standard state as a function
     *     of T and P
     */
    getStandardChemPotentials(mu);
    /*
     * Update the activity coefficients
     */
    s_update_lnActCoeff();
    doublereal RT = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(moleFractions_[k], SmallNumber);
        mu[k] += RT * (log(xx) + lnActCoeff_Scaled_[k]);
    }
}

doublereal MargulesVPSSTP::enthalpy_mole() const
{
    size_t kk = nSpecies();
    double h = 0;
    vector_fp hbar(kk);
    getPartialMolarEnthalpies(&hbar[0]);
    for (size_t i = 0; i < kk; i++) {
        h += moleFractions_[i]*hbar[i];
    }
    return h;
}

doublereal MargulesVPSSTP::entropy_mole() const
{
    size_t kk = nSpecies();
    double s = 0;
    vector_fp sbar(kk);
    getPartialMolarEntropies(&sbar[0]);
    for (size_t i = 0; i < kk; i++) {
        s += moleFractions_[i]*sbar[i];
    }
    return s;
}

doublereal MargulesVPSSTP::cp_mole() const
{
    size_t kk = nSpecies();
    double cp = 0;
    vector_fp cpbar(kk);
    getPartialMolarCp(&cpbar[0]);
    for (size_t i = 0; i < kk; i++) {
        cp += moleFractions_[i]*cpbar[i];
    }
    return cp;
}

doublereal MargulesVPSSTP::cv_mole() const
{
    return cp_mole() - GasConstant;
}

void MargulesVPSSTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * dimensionalize it.
     */
    double T = temperature();
    double RT = GasConstant * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();
    double RTT = RT * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] -= RTT * dlnActCoeffdT_Scaled_[k];
    }
}

void MargulesVPSSTP::getPartialMolarCp(doublereal* cpbar) const
{
    /*
     * Get the nondimensional standard state entropies
     */
    getCp_R(cpbar);
    double T = temperature();
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] -= 2 * T * dlnActCoeffdT_Scaled_[k] + T * T * d2lnActCoeffdT2_Scaled_[k];
    }
    /*
     * dimensionalize it.
     */
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void MargulesVPSSTP::getPartialMolarEntropies(doublereal* sbar) const
{
    double xx;
    /*
     * Get the nondimensional standard state entropies
     */
    getEntropy_R(sbar);
    double T = temperature();
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnActCoeff();
    s_update_dlnActCoeff_dT();

    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(moleFractions_[k], SmallNumber);
        sbar[k] += - lnActCoeff_Scaled_[k] -log(xx) - T * dlnActCoeffdT_Scaled_[k];
    }
    /*
     * dimensionalize it.
     */
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void MargulesVPSSTP::getPartialMolarVolumes(doublereal* vbar) const
{

    size_t iA, iB, delAK, delBK;
    double XA, XB, g0 , g1;
    double T = temperature();

    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);

    for (size_t iK = 0; iK < m_kk; iK++) {
        delAK = 0;
        delBK = 0;
    }
    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        iA =  m_pSpecies_A_ij[i];
        iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        g0 = (m_VHE_b_ij[i] - T * m_VSE_b_ij[i]);
        g1 = (m_VHE_c_ij[i] - T * m_VSE_c_ij[i]);
        const doublereal temp1 = g0 + g1 * XB;
        const doublereal all = -1.0*XA*XB*temp1 - XA*XB*XB*g1;

        for (size_t iK = 0; iK < m_kk; iK++) {
            vbar[iK] += all;
            // vbar[iK] += XA*XB*temp1+((delAK-XA)*XB+XA*(delBK-XB))*temp1+XB*XA*(delBK-XB)*g1;
        }
        vbar[iA] += XB * temp1;
        vbar[iB] += XA * temp1 + XA*XB*g1;
    }
}

doublereal MargulesVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("MargulesVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}

void MargulesVPSSTP::initThermo()
{
    initLengths();
    GibbsExcessVPSSTP::initThermo();
}

void  MargulesVPSSTP::initLengths()
{
    m_kk = nSpecies();
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
}

void MargulesVPSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    string stemp;
    string subname = "MargulesVPSSTP::initThermoXML";
    if ((int) id_.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id_) {
            throw CanteraError(subname, "phasenode and Id are incompatible");
        }
    }

    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError(subname,
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Make sure that the thermo model is Margules
     */
    stemp = thermoNode.attrib("model");
    string formString = lowercase(stemp);
    if (formString != "margules") {
        throw CanteraError(subname,
                           "model name isn't Margules: " + formString);

    }

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;
        string mStringa = acNode.attrib("model");
        string mString = lowercase(mStringa);
        if (mString != "margules") {
            throw CanteraError(subname.c_str(),
                               "Unknown activity coefficient model: " + mStringa);
        }
        for (size_t i = 0; i < acNodePtr->nChildren(); i++) {
            XML_Node& xmlACChild = acNodePtr->child(i);
            stemp = xmlACChild.name();
            string nodeName = lowercase(stemp);
            /*
             * Process a binary salt field, or any of the other XML fields
             * that make up the Pitzer Database. Entries will be ignored
             * if any of the species in the entry isn't in the solution.
             */
            if (nodeName == "binaryneutralspeciesparameters") {
                readXMLBinarySpecies(xmlACChild);

            }
        }
    }

    /*
     * Go down the chain
     */
    GibbsExcessVPSSTP::initThermoXML(phaseNode, id_);


}

void MargulesVPSSTP::s_update_lnActCoeff() const
{
    size_t iA, iB, iK;
    double XA, XB, g0 , g1;
    double T = temperature();
    double invRT = 1.0 / (GasConstant*T);
    lnActCoeff_Scaled_.resize(m_kk);
    for (iK = 0; iK < m_kk; iK++) {
        lnActCoeff_Scaled_[iK] = 0.0;
    }
    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        iA =  m_pSpecies_A_ij[i];
        iB =  m_pSpecies_B_ij[i];
        g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) * invRT;
        g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) * invRT;
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        const doublereal XAXB = XA * XB;
        const doublereal g0g1XB = (g0 + g1 * XB);
        const doublereal all = -1.0 * XAXB * g0g1XB - XAXB * XB * g1;
        for (iK = 0; iK < m_kk; iK++) {
            lnActCoeff_Scaled_[iK] += all;
        }
        lnActCoeff_Scaled_[iA] += XB * g0g1XB;
        lnActCoeff_Scaled_[iB] += XA * g0g1XB + XAXB * g1;
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dT() const
{
    size_t iA, iB, iK;
    doublereal XA, XB, g0, g1;
    doublereal invT = 1.0 / temperature();
    doublereal invRTT = 1.0 / (GasConstant)*invT*invT;
    dlnActCoeffdT_Scaled_.resize(m_kk);
    d2lnActCoeffdT2_Scaled_.resize(m_kk);
    for (iK = 0; iK < m_kk; iK++) {
        dlnActCoeffdT_Scaled_[iK] = 0.0;
        d2lnActCoeffdT2_Scaled_[iK] = 0.0;
    }
    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        iA =  m_pSpecies_A_ij[i];
        iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        g0 = -m_HE_b_ij[i] * invRTT;
        g1 = -m_HE_c_ij[i] * invRTT;
        const doublereal XAXB = XA * XB;
        const doublereal g0g1XB = (g0 + g1 * XB);
        const doublereal all = -1.0 * XAXB * g0g1XB - XAXB * XB * g1;
        const doublereal mult = 2.0 * invT;
        const doublereal dT2all = mult * all;
        for (iK = 0; iK < m_kk; iK++) {
            // double temp = (delAK * XB + XA * delBK - XA * XB) * (g0 + g1 * XB) + XA * XB * (delBK - XB) * g1;
            dlnActCoeffdT_Scaled_[iK] += all;
            d2lnActCoeffdT2_Scaled_[iK] -= dT2all;
        }
        dlnActCoeffdT_Scaled_[iA] += XB * g0g1XB;
        dlnActCoeffdT_Scaled_[iB] += XA * g0g1XB + XAXB * g1;
        d2lnActCoeffdT2_Scaled_[iA] -= mult * XB * g0g1XB;
        d2lnActCoeffdT2_Scaled_[iB] -= mult * (XA * g0g1XB + XAXB * g1);
    }
}

void MargulesVPSSTP::getdlnActCoeffdT(doublereal* dlnActCoeffdT) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdT[k] = dlnActCoeffdT_Scaled_[k];
    }
}

void MargulesVPSSTP::getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const
{
    s_update_dlnActCoeff_dT();
    for (size_t k = 0; k < m_kk; k++) {
        d2lnActCoeffdT2[k] = d2lnActCoeffdT2_Scaled_[k];
    }
}

void  MargulesVPSSTP::getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
                                       doublereal* dlnActCoeffds) const
{
    size_t iA, iB, iK;
    double XA, XB, g0 , g1, dXA, dXB;
    double T = temperature();
    double RT = GasConstant*T;

    //fvo_zero_dbl_1(dlnActCoeff, m_kk);
    s_update_dlnActCoeff_dT();

    for (iK = 0; iK < m_kk; iK++) {
        dlnActCoeffds[iK] = 0.0;
    }

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        iA =  m_pSpecies_A_ij[i];
        iB =  m_pSpecies_B_ij[i];
        XA = moleFractions_[iA];
        XB = moleFractions_[iB];
        dXA = dXds[iA];
        dXB = dXds[iB];
        g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
        g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;
        const doublereal g02g1XB = g0 + 2*g1*XB;
        const doublereal g2XAdXB = 2*g1*XA*dXB;
        const doublereal all = (-XB * dXA - XA *dXB) * g02g1XB - XB *g2XAdXB;
        for (iK = 0; iK < m_kk; iK++) {
            // dlnActCoeffds[iK] += ((delBK-XB)*dXA + (delAK-XA)*dXB)*(g0+2*g1*XB) + (delBK-XB)*2*g1*XA*dXB
            //                      + dlnActCoeffdT_Scaled_[iK]*dTds;
            dlnActCoeffds[iK] += all + dlnActCoeffdT_Scaled_[iK]*dTds;
        }
        dlnActCoeffds[iA] += dXB * g02g1XB;
        dlnActCoeffds[iB] += dXA * g02g1XB + g2XAdXB;
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dlnN_diag() const
{
    size_t iA, iB, iK, delAK, delBK;
    double XA, XB, XK, g0 , g1;
    double T = temperature();
    double RT = GasConstant*T;

    dlnActCoeffdlnN_diag_.assign(m_kk, 0.0);

    for (iK = 0; iK < m_kk; iK++) {

        XK = moleFractions_[iK];

        for (size_t i = 0; i <  numBinaryInteractions_; i++) {

            iA =  m_pSpecies_A_ij[i];
            iB =  m_pSpecies_B_ij[i];

            delAK = 0;
            delBK = 0;

            if (iA==iK) {
                delAK = 1;
            } else if (iB==iK) {
                delBK = 1;
            }

            XA = moleFractions_[iA];
            XB = moleFractions_[iB];

            g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
            g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

            dlnActCoeffdlnN_diag_[iK] += 2*(delBK-XB)*(g0*(delAK-XA)+g1*(2*(delAK-XA)*XB+XA*(delBK-XB)));

            // double gfac = g0 + g1 * XB;
            // double gggg = (delBK - XB) * g1;


            // dlnActCoeffdlnN_diag_[iK] += gfac * delAK * ( - XB + delBK);

            // dlnActCoeffdlnN_diag_[iK] += gfac * delBK * ( - XA + delAK);

            // dlnActCoeffdlnN_diag_[iK] += gfac * (2.0 * XA * XB - delAK * XB - XA * delBK);

            // dlnActCoeffdlnN_diag_[iK] += (delAK * XB + XA * delBK - XA * XB) * g1 * (-XB + delBK);

            // dlnActCoeffdlnN_diag_[iK] += gggg * ( - 2.0 * XA * XB + delAK * XB + XA * delBK);

            // dlnActCoeffdlnN_diag_[iK] += - g1 * XA * XB * (- XB + delBK);
        }
        dlnActCoeffdlnN_diag_[iK] = XK*dlnActCoeffdlnN_diag_[iK];//-XK;
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dlnN() const
{
    size_t iA, iB;
    doublereal delAK, delBK;
    double XA, XB, g0, g1,XM;
    double T = temperature();
    double RT = GasConstant*T;

    doublereal delAM, delBM;

    dlnActCoeffdlnN_.zero();

    /*
     *  Loop over the activity coefficient gamma_k
     */
    for (size_t iK = 0; iK < m_kk; iK++) {
        for (size_t iM = 0; iM < m_kk; iM++) {
            XM = moleFractions_[iM];
            for (size_t i = 0; i <  numBinaryInteractions_; i++) {

                iA =  m_pSpecies_A_ij[i];
                iB =  m_pSpecies_B_ij[i];

                delAK = 0.0;
                delBK = 0.0;
                delAM = 0.0;
                delBM = 0.0;
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

                XA = moleFractions_[iA];
                XB = moleFractions_[iB];

                g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
                g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

                dlnActCoeffdlnN_(iK,iM) += g0*((delAM-XA)*(delBK-XB)+(delAK-XA)*(delBM-XB));
                dlnActCoeffdlnN_(iK,iM) += 2*g1*((delAM-XA)*(delBK-XB)*XB+(delAK-XA)*(delBM-XB)*XB+(delBM-XB)*(delBK-XB)*XA);

                // double gfac = g0 + g1 * XB;
                // double gggg = (delBK - XB) * g1;


                // dlnActCoeffdlnN_(iK, iM) += gfac * delAK * ( - XB + delBM);

                // dlnActCoeffdlnN_(iK, iM) += gfac * delBK * ( - XA + delAM);

                // dlnActCoeffdlnN_(iK, iM) += gfac * (2.0 * XA * XB - delAM * XB - XA * delBM);

                // dlnActCoeffdlnN_(iK, iM) += (delAK * XB + XA * delBK - XA * XB) * g1 * (-XB + delBM);

                // dlnActCoeffdlnN_(iK, iM) += gggg * ( - 2.0 * XA * XB + delAM * XB + XA * delBM);

                // dlnActCoeffdlnN_(iK, iM) += - g1 * XA * XB * (- XB + delBM);
            }
            dlnActCoeffdlnN_(iK,iM) = XM*dlnActCoeffdlnN_(iK,iM);
        }
    }
}

void MargulesVPSSTP::s_update_dlnActCoeff_dlnX_diag() const
{
    doublereal T = temperature();
    dlnActCoeffdlnX_diag_.assign(m_kk, 0.0);
    doublereal RT = GasConstant * T;

    for (size_t i = 0; i <  numBinaryInteractions_; i++) {
        size_t iA =  m_pSpecies_A_ij[i];
        size_t iB =  m_pSpecies_B_ij[i];

        doublereal XA = moleFractions_[iA];
        doublereal XB = moleFractions_[iB];

        doublereal g0 = (m_HE_b_ij[i] - T * m_SE_b_ij[i]) / RT;
        doublereal g1 = (m_HE_c_ij[i] - T * m_SE_c_ij[i]) / RT;

        dlnActCoeffdlnX_diag_[iA] += XA*XB*(2*g1*-2*g0-6*g1*XB);
        dlnActCoeffdlnX_diag_[iB] += XA*XB*(2*g1*-2*g0-6*g1*XB);
    }
}

void MargulesVPSSTP::getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const
{
    s_update_dlnActCoeff_dlnN_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnN_diag[k] = dlnActCoeffdlnN_diag_[k];
    }
}

void MargulesVPSSTP::getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const
{
    s_update_dlnActCoeff_dlnX_diag();
    for (size_t k = 0; k < m_kk; k++) {
        dlnActCoeffdlnX_diag[k] = dlnActCoeffdlnX_diag_[k];
    }
}

void MargulesVPSSTP::getdlnActCoeffdlnN(const size_t ld, doublereal* dlnActCoeffdlnN)
{
    s_update_dlnActCoeff_dlnN();
    double* data =  & dlnActCoeffdlnN_(0,0);
    for (size_t k = 0; k < m_kk; k++) {
        for (size_t m = 0; m < m_kk; m++) {
            dlnActCoeffdlnN[ld * k + m] = data[m_kk * k + m];
        }
    }
}

void MargulesVPSSTP::resizeNumInteractions(const size_t num)
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

void MargulesVPSSTP::readXMLBinarySpecies(XML_Node& xmLBinarySpecies)
{
    string xname = xmLBinarySpecies.name();
    if (xname != "binaryNeutralSpeciesParameters") {
        throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies",
                           "Incorrect name for processing this routine: " + xname);
    }
    string stemp;
    size_t nParamsFound;
    vector_fp vParams;
    string iName = xmLBinarySpecies.attrib("speciesA");
    if (iName == "") {
        throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "no speciesA attrib");
    }
    string jName = xmLBinarySpecies.attrib("speciesB");
    if (jName == "") {
        throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "no speciesB attrib");
    }
    /*
     * Find the index of the species in the current phase. It's not
     * an error to not find the species
     */
    size_t iSpecies = speciesIndex(iName);
    if (iSpecies == npos) {
        return;
    }
    string ispName = speciesName(iSpecies);
    if (charge(iSpecies) != 0) {
        throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "speciesA charge problem");
    }
    size_t jSpecies = speciesIndex(jName);
    if (jSpecies == npos) {
        return;
    }
    string jspName = speciesName(jSpecies);
    if (charge(jSpecies) != 0) {
        throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies", "speciesB charge problem");
    }

    resizeNumInteractions(numBinaryInteractions_ + 1);
    size_t iSpot = numBinaryInteractions_ - 1;
    m_pSpecies_A_ij[iSpot] = iSpecies;
    m_pSpecies_B_ij[iSpot] = jSpecies;

    for (size_t iChild = 0; iChild < xmLBinarySpecies.nChildren(); iChild++) {
        XML_Node& xmlChild = xmLBinarySpecies.child(iChild);
        stemp = xmlChild.name();
        string nodeName = lowercase(stemp);
        /*
         * Process the binary species interaction child elements
         */
        if (nodeName == "excessenthalpy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessEnthalpy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies::excessEnthalpy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found. Need 2");
            }
            m_HE_b_ij[iSpot] = vParams[0];
            m_HE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessentropy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessEntropy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies::excessEntropy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found. Need 2");
            }
            m_SE_b_ij[iSpot] = vParams[0];
            m_SE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessvolume_enthalpy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessVolume_Enthalpy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies::excessVolume_Enthalpy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found. Need 2");
            }
            m_VHE_b_ij[iSpot] = vParams[0];
            m_VHE_c_ij[iSpot] = vParams[1];
        }

        if (nodeName == "excessvolume_entropy") {
            /*
             * Get the string containing all of the values
             */
            ctml::getFloatArray(xmlChild, vParams, true, "toSI", "excessVolume_Entropy");
            nParamsFound = vParams.size();

            if (nParamsFound != 2) {
                throw CanteraError("MargulesVPSSTP::readXMLBinarySpecies::excessVolume_Entropy for " + ispName
                                   + "::" + jspName,
                                   "wrong number of params found. Need 2");
            }
            m_VSE_b_ij[iSpot] = vParams[0];
            m_VSE_c_ij[iSpot] = vParams[1];
        }
    }
}

}
