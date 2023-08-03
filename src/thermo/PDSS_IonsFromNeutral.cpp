/**
 * @file PDSS_IonsFromNeutral.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{

PDSS_IonsFromNeutral::PDSS_IonsFromNeutral()
{
    warn_deprecated("class PDSS_IonsFromNeutral", "To be removed after Cantera 3.0");
}

void PDSS_IonsFromNeutral::setParent(VPStandardStateTP* phase, size_t k)
{
    neutralMoleculePhase_ = dynamic_cast<IonsFromNeutralVPSSTP&>(*phase).getNeutralMoleculePhase();
}

void PDSS_IonsFromNeutral::setNeutralSpeciesMultiplier(const string& species, double mult)
{
    neutralSpeciesMultipliers_[species] = mult;
    numMult_++;
}

void PDSS_IonsFromNeutral::setSpecialSpecies(bool special) {
    add2RTln2_ = !special;
}

void PDSS_IonsFromNeutral::getParameters(AnyMap& eosNode) const
{
    PDSS::getParameters(eosNode);
    eosNode["model"] = "ions-from-neutral-molecule";
    if (!add2RTln2_) {
        eosNode["special-species"] = true;
    }
    if (!neutralSpeciesMultipliers_.empty()) {
        eosNode["multipliers"] = neutralSpeciesMultipliers_;
    }
}

void PDSS_IonsFromNeutral::initThermo()
{
    PDSS::initThermo();
    if (m_input.getBool("special-species", false)) {
        setSpecialSpecies();
    }
    if (m_input.hasKey("multipliers")) {
        for (const auto& [species, multiplier] : m_input["multipliers"].asMap<double>()) {
            setNeutralSpeciesMultiplier(species, multiplier);
        }
    }

    m_p0 = neutralMoleculePhase_->refPressure();
    m_minTemp = neutralMoleculePhase_->minTemp();
    m_maxTemp = neutralMoleculePhase_->maxTemp();
    tmpNM.resize(neutralMoleculePhase_->nSpecies());
    for (auto [species, multiplier] : neutralSpeciesMultipliers_) {
        idNeutralMoleculeVec.push_back(neutralMoleculePhase_->speciesIndex(species));
        factorVec.push_back(multiplier);
    }
}

double PDSS_IonsFromNeutral::enthalpy_RT() const
{
    neutralMoleculePhase_->getEnthalpy_RT(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

double PDSS_IonsFromNeutral::intEnergy_mole() const
{
    return (m_h0_RT - 1.0) * GasConstant * m_temp;
}

double PDSS_IonsFromNeutral::entropy_R() const
{
    neutralMoleculePhase_->getEntropy_R(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val -= 2.0 * log(2.0);
    }
    return val;
}

double PDSS_IonsFromNeutral::gibbs_RT() const
{
    neutralMoleculePhase_->getGibbs_RT(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val += 2.0 * log(2.0);
    }
    return val;
}

double PDSS_IonsFromNeutral::cp_R() const
{
    neutralMoleculePhase_->getCp_R(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

double PDSS_IonsFromNeutral::molarVolume() const
{
    neutralMoleculePhase_->getStandardVolumes(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

double PDSS_IonsFromNeutral::density() const
{
    return (m_pres * m_mw / (GasConstant * m_temp));
}

double PDSS_IonsFromNeutral::gibbs_RT_ref() const
{
    neutralMoleculePhase_->getGibbs_RT_ref(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val += 2.0 * log(2.0);
    }
    return val;
}

double PDSS_IonsFromNeutral::enthalpy_RT_ref() const
{
    neutralMoleculePhase_->getEnthalpy_RT_ref(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

double PDSS_IonsFromNeutral::entropy_R_ref() const
{
    neutralMoleculePhase_->getEntropy_R_ref(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val -= 2.0 * log(2.0);
    }
    return val;
}

double PDSS_IonsFromNeutral::cp_R_ref() const
{
    neutralMoleculePhase_->getCp_R_ref(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

double PDSS_IonsFromNeutral::molarVolume_ref() const
{
    neutralMoleculePhase_->getStandardVolumes_ref(tmpNM.data());
    double val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

void PDSS_IonsFromNeutral::setState_TP(double temp, double pres)
{
    neutralMoleculePhase_->setState_TP(temp, pres);
    m_pres = pres;
    m_temp = temp;
}

}
