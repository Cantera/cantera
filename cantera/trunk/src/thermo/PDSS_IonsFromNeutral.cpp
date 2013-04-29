/**
 * @file PDSS_IonsFromNeutral.cpp
 * Implementation of a pressure dependent standard state
 * virtual function.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ct_defs.h"
#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"

#include <fstream>

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
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

PDSS_IonsFromNeutral& PDSS_IonsFromNeutral::operator=(const PDSS_IonsFromNeutral& b)
{
    if (&b == this) {
        return *this;
    }

    PDSS::operator=(b);

    /*
     *  The shallow pointer copy in the next step will be insufficient in most cases. However, its
     *  functionally the best we can do for this assignment operator. We fix up the pointer in the
     *  initAllPtrs() function.
     */
    neutralMoleculePhase_   = b.neutralMoleculePhase_;

    numMult_                = b.numMult_;
    idNeutralMoleculeVec    = b.idNeutralMoleculeVec;
    factorVec               = b.factorVec;
    add2RTln2_              = b.add2RTln2_;
    tmpNM                   = b.tmpNM;
    specialSpecies_         = b.specialSpecies_;

    return *this;
}

PDSS* PDSS_IonsFromNeutral::duplMyselfAsPDSS() const
{
    return new PDSS_IonsFromNeutral(*this);
}

void PDSS_IonsFromNeutral::initAllPtrs(VPStandardStateTP* tp, VPSSMgr* vpssmgr_ptr,
                                       SpeciesThermo* spthermo)
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
    std::string model = lowercase((*tn)["model"]);
    if (model != "ionfromneutral") {
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

    numMult_ = ctml::getPairs(*nsm,  key, val);
    idNeutralMoleculeVec.resize(numMult_);
    factorVec.resize(numMult_);
    tmpNM.resize(neutralMoleculePhase_->nSpecies());

    for (size_t i = 0; i < numMult_; i++) {
        idNeutralMoleculeVec[i] = neutralMoleculePhase_->speciesIndex(key[i]);
        factorVec[i] =  fpValueCheck(val[i]);
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
    if (inputFile.size() == 0) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSFile",
                           "input file is null");
    }
    std::string path = findInputFile(inputFile);
    ifstream fin(path.c_str());
    if (!fin) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSFile","could not open "
                           +path+" for reading.");
    }
    /*
     * The phase object automatically constructs an XML object.
     * Use this object to store information.
     */

    XML_Node* fxml = new XML_Node();
    fxml->build(fin);
    XML_Node* fxml_phase = findXMLPhase(fxml, id);
    if (!fxml_phase) {
        throw CanteraError("PDSS_IonsFromNeutral::constructPDSSFile",
                           "ERROR: Can not find phase named " +
                           id + " in file named " + inputFile);
    }

    XML_Node& speciesList = fxml_phase->child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &(fxml_phase->root()));
    const vector<string>&sss = tp->speciesNames();

    const XML_Node* s =  speciesDB->findByAttr("name", sss[spindex]);

    constructPDSSXML(tp, spindex, *s, *fxml_phase, id);
    delete fxml;
}

void PDSS_IonsFromNeutral::initThermoXML(const XML_Node& phaseNode, const std::string& id)
{
    PDSS::initThermoXML(phaseNode, id);
}

void PDSS_IonsFromNeutral::initThermo()
{
    PDSS::initThermo();
    SpeciesThermo& sp = m_tp->speciesThermo();
    m_p0 = sp.refPressure(m_spindex);
    m_minTemp = m_spthermo->minTemp(m_spindex);
    m_maxTemp = m_spthermo->maxTemp(m_spindex);
}

doublereal
PDSS_IonsFromNeutral::enthalpy_mole() const
{
    doublereal val = enthalpy_RT();
    doublereal RT = GasConstant * m_temp;
    return val * RT;
}

doublereal
PDSS_IonsFromNeutral::enthalpy_RT() const
{
    neutralMoleculePhase_->getEnthalpy_RT(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut = idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal
PDSS_IonsFromNeutral::intEnergy_mole() const
{
    doublereal val = m_h0_RT_ptr[m_spindex] - 1.0;
    doublereal RT = GasConstant * m_temp;
    return val * RT;
}

doublereal
PDSS_IonsFromNeutral::entropy_mole() const
{
    doublereal val = entropy_R();
    return val * GasConstant;
}

doublereal
PDSS_IonsFromNeutral::entropy_R() const
{
    neutralMoleculePhase_->getEntropy_R(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val -= 2.0 * log(2.0);
    }
    return val;
}

doublereal
PDSS_IonsFromNeutral::gibbs_mole() const
{
    doublereal val = gibbs_RT();
    doublereal RT = GasConstant * m_temp;
    return val * RT;
}

doublereal
PDSS_IonsFromNeutral::gibbs_RT() const
{
    neutralMoleculePhase_->getGibbs_RT(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val += 2.0 * log(2.0);
    }
    return val;
}

doublereal
PDSS_IonsFromNeutral::cp_mole() const
{
    doublereal val = cp_R();
    return val * GasConstant;
}

doublereal
PDSS_IonsFromNeutral::cp_R() const
{
    neutralMoleculePhase_->getCp_R(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal
PDSS_IonsFromNeutral::molarVolume() const
{
    neutralMoleculePhase_->getStandardVolumes(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal
PDSS_IonsFromNeutral::density() const
{
    return (m_pres * m_mw / (GasConstant * m_temp));
}

doublereal
PDSS_IonsFromNeutral::cv_mole() const
{
    throw CanteraError("PDSS_IonsFromNeutral::cv_mole()", "unimplemented");
    return 0.0;
}

doublereal
PDSS_IonsFromNeutral::gibbs_RT_ref() const
{
    neutralMoleculePhase_->getGibbs_RT_ref(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val += 2.0 * log(2.0);
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::enthalpy_RT_ref() const
{
    neutralMoleculePhase_->getEnthalpy_RT_ref(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::entropy_R_ref() const
{
    neutralMoleculePhase_->getEntropy_R_ref(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    if (add2RTln2_) {
        val -= 2.0 * log(2.0);
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::cp_R_ref() const
{
    neutralMoleculePhase_->getCp_R_ref(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal PDSS_IonsFromNeutral::molarVolume_ref() const
{
    neutralMoleculePhase_->getStandardVolumes_ref(DATA_PTR(tmpNM));
    doublereal val = 0.0;
    for (size_t i = 0; i < numMult_; i++) {
        size_t jNeut =  idNeutralMoleculeVec[i];
        val += factorVec[i] * tmpNM[jNeut];
    }
    return val;
}

doublereal  PDSS_IonsFromNeutral::pressure() const
{
    return m_pres;
}

void PDSS_IonsFromNeutral::setPressure(doublereal p)
{
    m_pres = p;
}

doublereal PDSS_IonsFromNeutral::critTemperature() const
{
    throw CanteraError("PDSS_IonsFromNeutral::critTemperature()", "unimplemented");
    return 0.0;
}

doublereal PDSS_IonsFromNeutral::critPressure() const
{
    throw CanteraError("PDSS_IonsFromNeutral::critPressure()", "unimplemented");
    return 0.0;
}

doublereal PDSS_IonsFromNeutral::critDensity() const
{
    throw CanteraError("PDSS_IonsFromNeutral::critDensity()", "unimplemented");
    return 0.0;
}

doublereal PDSS_IonsFromNeutral::temperature() const
{
    /*
     * Obtain the temperature from the owning VPStandardStateTP object if you can.
     */
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

void  PDSS_IonsFromNeutral::setState_TR(doublereal temp, doublereal rho)
{
}

doublereal PDSS_IonsFromNeutral::satPressure(doublereal t)
{
    throw CanteraError("PDSS_IonsFromNeutral::satPressure()", "unimplemented");
    /*NOTREACHED*/
    return 0.0;
}

}
