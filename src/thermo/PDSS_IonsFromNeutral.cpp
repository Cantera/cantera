/**
 * @file PDSS_IonsFromNeutral.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

PDSS_IonsFromNeutral::PDSS_IonsFromNeutral()
    : neutralMoleculePhase_(0)
    , numMult_(0)
    , add2RTln2_(true)
    , specialSpecies_(0)
{
}

void PDSS_IonsFromNeutral::setParent(VPStandardStateTP* phase, size_t k)
{
    neutralMoleculePhase_ = dynamic_cast<IonsFromNeutralVPSSTP&>(*phase).neutralMoleculePhase_;
}

void PDSS_IonsFromNeutral::setParametersFromXML(const XML_Node& speciesNode)
{
    PDSS::setParametersFromXML(speciesNode);
    const XML_Node* tn = speciesNode.findByName("thermo");
    if (!tn) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSXML",
                           "no thermo Node for species " + speciesNode.name());
    }
    if (!ba::iequals(tn->attrib("model"), "ionfromneutral")) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSXML",
                           "thermo model for species isn't IonsFromNeutral: "
                           + speciesNode.name());
    }
    const XML_Node* nsm = tn->findByName("neutralSpeciesMultipliers");
    if (!nsm) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSXML",
                           "no Thermo::neutralSpeciesMultipliers Node for species " + speciesNode.name());
    }

    neutralSpeciesMultipliers_ = parseCompString(nsm->value());
    numMult_ = neutralSpeciesMultipliers_.size();

    specialSpecies_ = 0;
    const XML_Node* ss = tn->findByName("specialSpecies");
    if (ss) {
        specialSpecies_ = 1;
    }
    const XML_Node* sss = tn->findByName("secondSpecialSpecies");
    if (sss) {
        specialSpecies_ = 2;
    }
    add2RTln2_ = true;
    if (specialSpecies_ == 1) {
        add2RTln2_ = false;
    }
}

void PDSS_IonsFromNeutral::initThermo()
{
    PDSS::initThermo();
    m_p0 = neutralMoleculePhase_->refPressure();
    m_minTemp = neutralMoleculePhase_->minTemp();
    m_maxTemp = neutralMoleculePhase_->maxTemp();
    tmpNM.resize(neutralMoleculePhase_->nSpecies());
    for (auto multiplier : neutralSpeciesMultipliers_) {
        idNeutralMoleculeVec.push_back( neutralMoleculePhase_->speciesIndex(multiplier.first));
        factorVec.push_back(multiplier.second);
    }
}

doublereal PDSS_IonsFromNeutral::enthalpy_RT() const
{
    neutralMoleculePhase_->getEnthalpy_RT(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::intEnergy_mole() const
{
    return (m_h0_RT - 1.0) * GasConstant * m_temp;
}

doublereal PDSS_IonsFromNeutral::entropy_R() const
{
    neutralMoleculePhase_->getEntropy_R(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val -= 2.0 * log(2.0);
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::gibbs_RT() const
{
    neutralMoleculePhase_->getGibbs_RT(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val += 2.0 * log(2.0);
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::cp_R() const
{
    neutralMoleculePhase_->getCp_R(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::molarVolume() const
{
    neutralMoleculePhase_->getStandardVolumes(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::density() const
{
    return (m_pres * m_mw / (GasConstant * m_temp));
}

doublereal PDSS_IonsFromNeutral::gibbs_RT_ref() const
{
    neutralMoleculePhase_->getGibbs_RT_ref(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val += 2.0 * log(2.0);
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::enthalpy_RT_ref() const
{
    neutralMoleculePhase_->getEnthalpy_RT_ref(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::entropy_R_ref() const
{
    neutralMoleculePhase_->getEntropy_R_ref(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val -= 2.0 * log(2.0);
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::cp_R_ref() const
{
    neutralMoleculePhase_->getCp_R_ref(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::molarVolume_ref() const
{
    neutralMoleculePhase_->getStandardVolumes_ref(tmpNM.data());
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

void PDSS_IonsFromNeutral::setState_TP(doublereal temp, doublereal pres)
{
    m_pres = pres;
    m_temp = temp;
}

void PDSS_IonsFromNeutral::setState_TR(doublereal temp, doublereal rho)
{
}

}
