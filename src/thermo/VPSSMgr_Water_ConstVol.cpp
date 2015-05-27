/**
 *  @file VPSSMgr_Water_ConstVol.cpp
 * Definition file for a derived class that handles the calculation
 * of standard state thermo properties for pure water and
 * a set of species which have a constant molar volume pressure
 * dependence.
 * (see \ref thermoprops and class
 * \link Cantera::VPSSMgr_Water_ConstVol VPSSMgr_Water_ConstVol\endlink).
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/VPSSMgr_Water_ConstVol.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{
VPSSMgr_Water_ConstVol::VPSSMgr_Water_ConstVol(VPStandardStateTP* vp_ptr,
        SpeciesThermo* spth) :
    VPSSMgr(vp_ptr, spth),
    m_waterSS(0)
{
    m_useTmpRefStateStorage      = true;
    m_useTmpStandardStateStorage = true;
}

VPSSMgr_Water_ConstVol::VPSSMgr_Water_ConstVol(const VPSSMgr_Water_ConstVol& right) :
    VPSSMgr(right.m_vptp_ptr, right.m_spthermo)
{
    m_useTmpRefStateStorage = true;
    m_useTmpStandardStateStorage = true;
    *this = right;
}

VPSSMgr_Water_ConstVol&
VPSSMgr_Water_ConstVol::operator=(const VPSSMgr_Water_ConstVol& b)
{
    if (&b == this) {
        return *this;
    }
    VPSSMgr::operator=(b);
    return *this;
}

VPSSMgr*
VPSSMgr_Water_ConstVol::duplMyselfAsVPSSMgr() const
{
    return new VPSSMgr_Water_ConstVol(*this);
}

void
VPSSMgr_Water_ConstVol::initAllPtrs(VPStandardStateTP* vp_ptr,
                                    SpeciesThermo* sp_ptr)
{
    VPSSMgr::initAllPtrs(vp_ptr, sp_ptr);
    m_waterSS = dynamic_cast<PDSS_Water*>(m_vptp_ptr->providePDSS(0));
    if (!m_waterSS) {
        throw CanteraError("VPSSMgr_Water_ConstVol::initAllPtrs",
                           "bad dynamic cast");
    }
}

void
VPSSMgr_Water_ConstVol::getEnthalpy_RT_ref(doublereal* hrt) const
{
    // Everything should be OK except for the water SS
    m_p0 = m_waterSS->pref_safe(m_tlast);
    if (m_p0 != m_plast) {
        m_waterSS->setState_TP(m_tlast, m_p0);
        m_h0_RT[0] = (m_waterSS->enthalpy_mole()) / (GasConstant * m_tlast);
        m_waterSS->setState_TP(m_tlast, m_plast);
    } else {
        m_h0_RT[0] = m_hss_RT[0];
    }
    copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
}

void
VPSSMgr_Water_ConstVol::getGibbs_RT_ref(doublereal* grt) const
{
    // Everything should be OK except for the water SS
    m_p0 = m_waterSS->pref_safe(m_tlast);
    if (m_p0 != m_plast) {
        m_waterSS->setState_TP(m_tlast, m_p0);
        m_g0_RT[0] = (m_waterSS->gibbs_mole()) / (GasConstant * m_tlast);
        m_waterSS->setState_TP(m_tlast, m_plast);
    } else {
        m_g0_RT[0] = m_gss_RT[0];
    }
    copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
}

void
VPSSMgr_Water_ConstVol::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= GasConstant * m_tlast;
    }
}

void
VPSSMgr_Water_ConstVol::getEntropy_R_ref(doublereal* sr) const
{
    // Everything should be OK except for the water SS
    m_p0 = m_waterSS->pref_safe(m_tlast);
    if (m_p0 != m_plast) {
        m_waterSS->setState_TP(m_tlast, m_p0);
        m_s0_R[0] = (m_waterSS->entropy_mole()) / GasConstant;
        m_waterSS->setState_TP(m_tlast, m_plast);
    } else {
        m_s0_R[0] = m_sss_R[0];
    }
    copy(m_s0_R.begin(), m_s0_R.end(), sr);
}

void
VPSSMgr_Water_ConstVol::getCp_R_ref(doublereal* cpr) const
{
    // Everything should be OK except for the water SS
    m_p0 = m_waterSS->pref_safe(m_tlast);
    if (m_p0 != m_plast) {
        m_waterSS->setState_TP(m_tlast, m_p0);
        m_cp0_R[0] = (m_waterSS->cp_mole()) / GasConstant;
        m_waterSS->setState_TP(m_tlast, m_plast);
    } else {
        m_cp0_R[0] = m_cpss_R[0];
    }
    copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
}

void
VPSSMgr_Water_ConstVol::getStandardVolumes_ref(doublereal* vol) const
{
    // Everything should be OK except for the water SS
    m_p0 = m_waterSS->pref_safe(m_tlast);
    if (m_p0 != m_plast) {
        m_waterSS->setState_TP(m_tlast, m_p0);
        m_V0[0] = m_vptp_ptr->molecularWeight(0) / m_waterSS->density();
        m_waterSS->setState_TP(m_tlast, m_plast);
    } else {
        m_V0[0] = m_Vss[0];
    }
    copy(m_V0.begin(), m_V0.end(), vol);
}

void VPSSMgr_Water_ConstVol::_updateRefStateThermo() const
{
    m_p0 = m_waterSS->pref_safe(m_tlast);
    m_spthermo->update(m_tlast, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
    for (size_t k = 0; k < m_kk; k++) {
        m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        m_vptp_ptr->providePDSS(k)->setTemperature(m_tlast);
    }
    m_waterSS->setState_TP(m_tlast, m_p0);
    m_h0_RT[0] = (m_waterSS->enthalpy_mole()) / (GasConstant * m_tlast);
    m_s0_R[0]  = (m_waterSS->entropy_mole()) / GasConstant;
    m_cp0_R[0] = (m_waterSS->cp_mole()) / GasConstant;
    m_g0_RT[0] = (m_hss_RT[0] - m_sss_R[0]);
    m_V0[0]    =  m_vptp_ptr->molecularWeight(0) / (m_waterSS->density());
    m_waterSS->setState_TP(m_tlast, m_plast);
}

void VPSSMgr_Water_ConstVol::_updateStandardStateThermo()
{
    doublereal del_pRT = (m_plast - OneAtm) / (GasConstant * m_tlast);

    for (size_t k = 1; k < m_kk; k++) {
        m_hss_RT[k]  = m_h0_RT[k] + del_pRT * m_Vss[k];
        m_cpss_R[k]  = m_cp0_R[k];
        m_sss_R[k]   = m_s0_R[k];
        m_gss_RT[k]  = m_hss_RT[k] - m_sss_R[k];
        // m_Vss[k] constant
        PDSS* kPDSS = m_vptp_ptr->providePDSS(k);
        kPDSS->setState_TP(m_tlast, m_plast);
    }
    // Do the water
    m_waterSS->setState_TP(m_tlast, m_plast);
    m_hss_RT[0] = (m_waterSS->enthalpy_mole()) / (GasConstant * m_tlast);
    m_sss_R[0]  = (m_waterSS->entropy_mole()) / GasConstant;
    m_cpss_R[0] = (m_waterSS->cp_mole())      / GasConstant;
    m_gss_RT[0] = (m_hss_RT[0] - m_sss_R[0]);
    m_Vss[0]    = (m_vptp_ptr->molecularWeight(0) / m_waterSS->density());
}

void VPSSMgr_Water_ConstVol::initThermo()
{
    VPSSMgr::initThermo();
}

void
VPSSMgr_Water_ConstVol::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    VPSSMgr::initThermoXML(phaseNode, id);

    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &phaseNode.root());

    if (!m_waterSS) {
        throw CanteraError("VPSSMgr_Water_ConstVol::initThermoXML",
                           "bad dynamic cast");
    }

    m_waterSS->setState_TP(300., OneAtm);
    m_Vss[0] = (m_waterSS->density())      / m_vptp_ptr->molecularWeight(0);

    for (size_t k = 1; k < m_kk; k++) {
        const XML_Node* s = speciesDB->findByAttr("name", m_vptp_ptr->speciesName(k));
        if (!s) {
            throw CanteraError("VPSSMgr_Water_ConstVol::initThermoXML",
                               "no species Node for species " + m_vptp_ptr->speciesName(k));
        }
        const XML_Node* ss = s->findByName("standardState");
        if (!ss) {
            throw CanteraError("VPSSMgr_Water_ConstVol::initThermoXML",
                               "no standardState Node for species " + s->attrib("name"));
        }
        if (ss->attrib("model") != "constant_incompressible") {
            throw CanteraError("VPSSMgr_Water_ConstVol::initThermoXML",
                               "standardState model for species isn't "
                               "constant_incompressible: " + s->attrib("name"));
        }
        m_Vss[k] = getFloat(*ss, "molarVolume", "toSI");
    }
}

PDSS*
VPSSMgr_Water_ConstVol::createInstallPDSS(size_t k, const XML_Node& speciesNode,
        const XML_Node* const phaseNode_ptr)
{
    PDSS* kPDSS = 0;
    // Will have to do something for water
    // -> make sure it's species 0
    // -> make sure it's designated as a real water EOS
    if (k == 0) {
        string xn = speciesNode["name"];
        if (xn != "H2O(L)") {
            throw CanteraError("VPSSMgr_Water_ConstVol::installSpecies",
                               "h2o wrong name: " + xn);
        }
        const XML_Node* ss = speciesNode.findByName("standardState");
        std::string model = ss->attrib("model");
        if (model != "waterIAPWS" && model != "waterPDSS") {
            throw CanteraError("VPSSMgr_Water_ConstVol::installSpecies",
                               "wrong SS mode: " + model);
        }
        delete m_waterSS;
        m_waterSS = new PDSS_Water(m_vptp_ptr, 0);
        GeneralSpeciesThermo* genSpthermo = dynamic_cast<GeneralSpeciesThermo*>(m_spthermo);
        if (!genSpthermo) {
            throw CanteraError("VPSSMgr_Water_ConstVol::installSpecies",
                               "failed dynamic cast");
        }
        genSpthermo->installPDSShandler(k, m_waterSS, this);
        kPDSS = m_waterSS;
    } else {

        VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);

        const XML_Node* ss = speciesNode.findByName("standardState");
        if (!ss) {
            throw CanteraError("VPSSMgr_Water_ConstVol::installSpecies",
                               "no standardState Node for species " + speciesNode.name());
        }
        if (ss->attrib("model") != "constant_incompressible") {
            throw CanteraError("VPSSMgr_Water_ConstVol::initThermoXML",
                               "standardState model for species isn't "
                               "constant_incompressible: " + speciesNode.name());
        }
        if (m_Vss.size() < k+1) {
            m_Vss.resize(k+1, 0.0);
        }
        m_Vss[k] = getFloat(*ss, "molarVolume", "toSI");

        // instantiate a new kPDSS object
        kPDSS = new PDSS_ConstVol(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
    }
    return kPDSS;
}

PDSS_enumType VPSSMgr_Water_ConstVol::reportPDSSType(int k) const
{
    return cPDSS_UNDEF;
}

VPSSMgr_enumType VPSSMgr_Water_ConstVol::reportVPSSMgrType() const
{
    return cVPSSMGR_WATER_CONSTVOL;
}
}
