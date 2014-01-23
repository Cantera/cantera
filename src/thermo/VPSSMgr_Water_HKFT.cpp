/**
 *  @file VPSSMgr_Water_HKFT.cpp
 * Definition file for a derived class that handles the calculation
 * of standard state thermo properties for pure water and
 *  a set of species which obey the HKFT standard state
 * dependence
 * (see \ref thermoprops and class
 * \link Cantera::VPSSMgr_Water_HKFT VPSSMgr_Water_HKFT\endlink).
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/VPSSMgr_Water_HKFT.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_HKFT.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"

using namespace std;

namespace Cantera
{

VPSSMgr_Water_HKFT::VPSSMgr_Water_HKFT(VPStandardStateTP* vp_ptr,
                                       SpeciesThermo* spth) :
    VPSSMgr(vp_ptr, spth),
    m_waterSS(0),
    m_tlastRef(-1.0)
{
    m_useTmpRefStateStorage      = true;
    m_useTmpStandardStateStorage = true;
}

VPSSMgr_Water_HKFT::VPSSMgr_Water_HKFT(const VPSSMgr_Water_HKFT& right) :
    VPSSMgr(right.m_vptp_ptr, right.m_spthermo),
    m_waterSS(0),
    m_tlastRef(-1.0)
{
    m_useTmpRefStateStorage = true;
    m_useTmpStandardStateStorage = true;
    *this = right;
}


VPSSMgr_Water_HKFT&
VPSSMgr_Water_HKFT::operator=(const VPSSMgr_Water_HKFT& b)
{
    if (&b == this) {
        return *this;
    }
    VPSSMgr::operator=(b);
    m_waterSS = (PDSS_Water*) m_vptp_ptr->providePDSS(0);
    m_tlastRef = -1.0;
    return *this;
}

VPSSMgr*
VPSSMgr_Water_HKFT::duplMyselfAsVPSSMgr() const
{
    return new VPSSMgr_Water_HKFT(*this);
}

void
VPSSMgr_Water_HKFT::getEnthalpy_RT_ref(doublereal* hrt) const
{
    updateRefStateThermo();
    copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
}

void
VPSSMgr_Water_HKFT::getGibbs_RT_ref(doublereal* grt) const
{
    updateRefStateThermo();
    copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
}

void
VPSSMgr_Water_HKFT::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    doublereal RT = GasConstant * m_tlast;
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT;
    }
}

void
VPSSMgr_Water_HKFT::getEntropy_R_ref(doublereal* sr) const
{
    updateRefStateThermo();
    copy(m_s0_R.begin(), m_s0_R.end(), sr);
}

void
VPSSMgr_Water_HKFT::getCp_R_ref(doublereal* cpr) const
{
    updateRefStateThermo();
    copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
}

void
VPSSMgr_Water_HKFT::getStandardVolumes_ref(doublereal* vol) const
{
    updateRefStateThermo();
    copy(m_V0.begin(), m_V0.end(), vol);
}

void VPSSMgr_Water_HKFT::setState_P(doublereal pres)
{
    if (m_plast != pres) {
        m_plast = pres;
        _updateStandardStateThermo();
    }
}

void VPSSMgr_Water_HKFT::setState_T(doublereal temp)
{
    if (m_tlast != temp) {
        m_tlast = temp;
        _updateStandardStateThermo();
    }
}

void VPSSMgr_Water_HKFT::setState_TP(doublereal temp, doublereal pres)
{
    if (m_tlast != temp) {
        m_tlast = temp;
        m_plast = pres;
        _updateStandardStateThermo();
    } else if (m_plast != pres) {
        m_plast = pres;
        _updateStandardStateThermo();
    }
}

void VPSSMgr_Water_HKFT::updateRefStateThermo() const
{
    if (m_tlastRef != m_tlast) {
        m_tlastRef = m_tlast;
        _updateRefStateThermo();
    }
}

void VPSSMgr_Water_HKFT::_updateRefStateThermo() const
{
    m_p0 = m_waterSS->pref_safe(m_tlast);
    doublereal RT = GasConstant * m_tlast;
    m_waterSS->setState_TP(m_tlast, m_p0);
    m_h0_RT[0] = (m_waterSS->enthalpy_mole())/ RT;
    m_s0_R[0]  = (m_waterSS->entropy_mole()) / GasConstant;
    m_cp0_R[0] = (m_waterSS->cp_mole()) / GasConstant;
    m_g0_RT[0] = (m_hss_RT[0] - m_sss_R[0]);
    m_V0[0]    = (m_waterSS->density()) / m_vptp_ptr->molecularWeight(0);
    PDSS_HKFT* ps;
    for (size_t k = 1; k < m_kk; k++) {
        ps = (PDSS_HKFT*) m_vptp_ptr->providePDSS(k);
        ps->setState_TP(m_tlast, m_p0);
        m_cp0_R[k]  = ps->cp_R();
        m_s0_R[k]   = ps->entropy_mole() / GasConstant;
        m_g0_RT[k]  = ps->gibbs_RT();
        m_h0_RT[k]  = m_g0_RT[k] + m_s0_R[k];
#ifdef DEBUG_MODE_NOT
        double h = ps->enthalpy_RT();
        if (fabs(m_h0_RT[k] - h) > 1.0E-4) {
            printf(" VPSSMgr_Water_HKFT::_updateRefStateThermo:: we have a discrepancy\n");
        }
#endif
        m_V0[k]     = ps->molarVolume();

    }
    m_waterSS->setState_TP(m_tlast, m_plast);
    for (size_t k = 1; k < m_kk; k++) {
        ps = (PDSS_HKFT*) m_vptp_ptr->providePDSS(k);
        ps->setState_TP(m_tlast, m_plast);
    }
}

void VPSSMgr_Water_HKFT::_updateStandardStateThermo()
{
    doublereal RT = GasConstant * m_tlast;
    // Do the water
    m_waterSS->setState_TP(m_tlast, m_plast);
    m_hss_RT[0] = (m_waterSS->enthalpy_mole())/ RT;
    m_sss_R[0]  = (m_waterSS->entropy_mole()) / GasConstant;
    m_cpss_R[0] = (m_waterSS->cp_mole())      / GasConstant;
    m_gss_RT[0] = (m_hss_RT[0] - m_sss_R[0]);
    m_Vss[0]    = (m_vptp_ptr->molecularWeight(0)) / (m_waterSS->density());

    for (size_t k = 1; k < m_kk; k++) {
        PDSS_HKFT* ps = (PDSS_HKFT*) m_vptp_ptr->providePDSS(k);
        ps->setState_TP(m_tlast, m_plast);
        m_cpss_R[k]  = ps->cp_R();
        m_sss_R[k]   = ps->entropy_R();
        m_gss_RT[k]  = ps->gibbs_RT();
        m_hss_RT[k]  = m_gss_RT[k] + m_sss_R[k];
        m_Vss[k]     = ps->molarVolume();
    }
}

void VPSSMgr_Water_HKFT::initThermo()
{
    VPSSMgr::initThermo();
}


void
VPSSMgr_Water_HKFT::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    VPSSMgr::initThermoXML(phaseNode, id);

    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &phaseNode.root());
    const vector<string> &sss = m_vptp_ptr->speciesNames();

    m_waterSS->setState_TP(300., OneAtm);
    m_Vss[0] = (m_waterSS->density())      / m_vptp_ptr->molecularWeight(0);

    for (size_t k = 1; k < m_kk; k++) {
        const XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
        if (!s) {
            throw CanteraError("VPSSMgr_Water_HKFT::initThermoXML",
                               "No species Node for species " + sss[k]);
        }
        const XML_Node* ss = s->findByName("standardState");
        if (!ss) {
            throw CanteraError("VPSSMgr_Water_HKFT::initThermoXML",
                               "No standardState Node for species " + sss[k]);
        }
        std::string model = lowercase((*ss)["model"]);
        if (model != "hkft") {
            throw CanteraError("VPSSMgr_Water_HKFT::initThermoXML",
                               "Standard state model for a solute species isn't "
                               "the HKFT standard state model: " + sss[k]);
        }
    }
}

PDSS*
VPSSMgr_Water_HKFT::createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                      const XML_Node* const phaseNode_ptr)
{
    PDSS* kPDSS = 0;

    const XML_Node* ss = speciesNode.findByName("standardState");
    if (!ss) {
        std::string sName = speciesNode["name"];
        throw CanteraError("VPSSMgr_Water_HKFT::installSpecies",
                           "No standardState Node for species " + sName);
    }
    // Will have to do something for water
    // -> make sure it's species 0
    // -> make sure it's designated as a real water EOS
    if (k == 0) {
        string xn = speciesNode["name"];
        if (xn != "H2O(L)") {
            throw CanteraError("VPSSMgr_Water_HKFT::installSpecies",
                               "h2o wrong name: " + xn);
        }

        std::string model = (*ss)["model"];
        if (model != "waterIAPWS" && model != "waterPDSS") {
            throw CanteraError("VPSSMgr_Water_HKFT::installSpecies",
                               "wrong SS mode: " + model);
        }
        //VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);
        delete m_waterSS;
        m_waterSS = new PDSS_Water(m_vptp_ptr, 0);

        GeneralSpeciesThermo* genSpthermo = dynamic_cast<GeneralSpeciesThermo*>(m_spthermo);
        if (!genSpthermo) {
            throw CanteraError("VPSSMgr_Water_HKFT::installSpecies",
                               "failed dynamic cast");
        }
        genSpthermo->installPDSShandler(k, m_waterSS, this);

        kPDSS = m_waterSS;
    } else {
        std::string model = (*ss)["model"];
        if (model != "HKFT") {
            std::string sName = speciesNode["name"];
            throw CanteraError("VPSSMgr_Water_HKFT::initThermoXML",
                               "standardState model for species isn't "
                               "HKFT: " + sName);
        }

        kPDSS = new PDSS_HKFT(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);

        GeneralSpeciesThermo* genSpthermo = dynamic_cast<GeneralSpeciesThermo*>(m_spthermo);
        if (!genSpthermo) {
            throw CanteraError("VPSSMgr_Water_HKFT::installSpecies",
                               "failed dynamic cast");
        }
        genSpthermo->installPDSShandler(k, kPDSS, this);
    }
    return kPDSS;
}

void
VPSSMgr_Water_HKFT::initAllPtrs(VPStandardStateTP* vp_ptr,
                                    SpeciesThermo* sp_ptr)
{
    VPSSMgr::initAllPtrs(vp_ptr, sp_ptr);
    m_waterSS = dynamic_cast<PDSS_Water*>(m_vptp_ptr->providePDSS(0));
    if (!m_waterSS) {
        throw CanteraError("VPSSMgr_Water_ConstVol::initAllPtrs",
                           "bad dynamic cast");
    }
}

PDSS_enumType VPSSMgr_Water_HKFT::reportPDSSType(int k) const
{
    return cPDSS_UNDEF;
}

VPSSMgr_enumType VPSSMgr_Water_HKFT::reportVPSSMgrType() const
{
    return cVPSSMGR_WATER_HKFT;
}
}
