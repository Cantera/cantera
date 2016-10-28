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
PDSS_IonsFromNeutral::PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex) :
    PDSS(tp, spindex),
    neutralMoleculePhase_(0),
    numMult_(0),
    add2RTln2_(true),
    specialSpecies_(0)
{
    m_pdssType = cPDSS_IONSFROMNEUTRAL;
}

PDSS_IonsFromNeutral::PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex,
        const std::string& inputFile, const std::string& id) :
    PDSS(tp, spindex),
    neutralMoleculePhase_(0),
    numMult_(0),
    add2RTln2_(true),
    specialSpecies_(0)
{
    warn_deprecated("PDSS_IonsFromNeutral constructor from XML input file",
                    "To be removed after Cantera 2.3.");
    m_pdssType = cPDSS_IONSFROMNEUTRAL;
    constructPDSSFile(tp, spindex, inputFile, id);
}

PDSS_IonsFromNeutral::PDSS_IonsFromNeutral(VPStandardStateTP* tp, size_t spindex, const XML_Node& speciesNode,
        const XML_Node& phaseRoot, bool spInstalled) :
    PDSS(tp, spindex),
    neutralMoleculePhase_(0),
    numMult_(0),
    add2RTln2_(true),
    specialSpecies_(0)
{
    if (!spInstalled) {
        throw CanteraError("PDSS_IonsFromNeutral", "sp installing not done yet");
    }
    m_pdssType = cPDSS_IONSFROMNEUTRAL;
    std::string id = "";
    constructPDSSXML(tp, spindex, speciesNode, phaseRoot, id);
}

PDSS_IonsFromNeutral::PDSS_IonsFromNeutral(const PDSS_IonsFromNeutral& b) :
    PDSS(b)
{
    // Use the assignment operator to do the brunt of the work for the copy
    // constructor.
    *this = b;
}

PDSS_IonsFromNeutral& PDSS_IonsFromNeutral::operator=(const PDSS_IonsFromNeutral& b)
{
    if (&b == this) {
        return *this;
    }

    PDSS::operator=(b);

    // The shallow pointer copy in the next step will be insufficient in most
    // cases. However, its functionally the best we can do for this assignment
    // operator. We fix up the pointer in the initAllPtrs() function.
    neutralMoleculePhase_ = b.neutralMoleculePhase_;

    numMult_ = b.numMult_;
    idNeutralMoleculeVec = b.idNeutralMoleculeVec;
    factorVec = b.factorVec;
    add2RTln2_ = b.add2RTln2_;
    tmpNM = b.tmpNM;
    specialSpecies_ = b.specialSpecies_;

    return *this;
}

PDSS* PDSS_IonsFromNeutral::duplMyselfAsPDSS() const
{
    return new PDSS_IonsFromNeutral(*this);
}

void PDSS_IonsFromNeutral::initAllPtrs(VPStandardStateTP* tp, VPSSMgr* vpssmgr_ptr,
                                       MultiSpeciesThermo* spthermo)
{
    PDSS::initAllPtrs(tp, vpssmgr_ptr, spthermo);

    IonsFromNeutralVPSSTP* ionPhase = dynamic_cast<IonsFromNeutralVPSSTP*>(tp);
    if (!ionPhase) {
        throw CanteraError("PDSS_IonsFromNeutral::initAllPts", "Dynamic cast failed");
    }
    neutralMoleculePhase_ = ionPhase->neutralMoleculePhase_;
}

void PDSS_IonsFromNeutral::constructPDSSXML(VPStandardStateTP* tp, size_t spindex,
        const XML_Node& speciesNode,
        const XML_Node& phaseNode, const std::string& id)
{
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

    IonsFromNeutralVPSSTP* ionPhase = dynamic_cast<IonsFromNeutralVPSSTP*>(tp);
    if (!ionPhase) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSXML", "Dynamic cast failed");
    }
    neutralMoleculePhase_ = ionPhase->neutralMoleculePhase_;

    std::vector<std::string> key;
    std::vector<std::string> val;
    numMult_ = getPairs(*nsm, key, val);
    idNeutralMoleculeVec.resize(numMult_);
    factorVec.resize(numMult_);
    tmpNM.resize(neutralMoleculePhase_->nSpecies());
    for (size_t i = 0; i < numMult_; i++) {
        idNeutralMoleculeVec[i] = neutralMoleculePhase_->speciesIndex(key[i]);
        factorVec[i] = fpValueCheck(val[i]);
    }
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

void PDSS_IonsFromNeutral::constructPDSSFile(VPStandardStateTP* tp, size_t spindex,
        const std::string& inputFile, const std::string& id)
{
    warn_deprecated("PDSS_IonsFromNeutral::constructPDSSFile",
                    "To be removed after Cantera 2.3.");
    if (inputFile.size() == 0) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSFile",
                           "input file is null");
    }

    // The phase object automatically constructs an XML object. Use this object
    // to store information.
    XML_Node fxml;
    fxml.build(findInputFile(inputFile));
    XML_Node* fxml_phase = findXMLPhase(&fxml, id);
    if (!fxml_phase) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSFile",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }

    XML_Node& speciesList = fxml_phase->child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &fxml_phase->root());
    const XML_Node* s = speciesDB->findByAttr("name", tp->speciesName(spindex));
    constructPDSSXML(tp, spindex, *s, *fxml_phase, id);
}

void PDSS_IonsFromNeutral::initThermo()
{
    PDSS::initThermo();
    m_p0 = m_tp->speciesThermo().refPressure(m_spindex);
    m_minTemp = m_spthermo->minTemp(m_spindex);
    m_maxTemp = m_spthermo->maxTemp(m_spindex);
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
    return (m_h0_RT_ptr[m_spindex] - 1.0) * GasConstant * m_temp;
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

doublereal PDSS_IonsFromNeutral::temperature() const
{
    // Obtain the temperature from the owning VPStandardStateTP object if you
    // can.
    m_temp = m_vpssmgr_ptr->temperature();
    return m_temp;
}

void PDSS_IonsFromNeutral::setTemperature(doublereal temp)
{
    m_temp = temp;
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
