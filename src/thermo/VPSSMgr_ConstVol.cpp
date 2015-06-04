/**
 *  @file VPSSMgr_ConstVol.cpp
 *  Definition file for a derived class that handles the calculation
 *  of standard state thermo properties for
 *  a set of species which have a constant molar volume pressure
 *  dependence (see \ref thermoprops and
 * class \link Cantera::VPSSMgr_ConstVol VPSSMgr_ConstVol\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/VPSSMgr_ConstVol.h"
#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

VPSSMgr_ConstVol::VPSSMgr_ConstVol(VPStandardStateTP* vp_ptr, SpeciesThermo* spth) :
    VPSSMgr(vp_ptr, spth)
{
    m_useTmpRefStateStorage      = true;
    m_useTmpStandardStateStorage = true;
}

VPSSMgr_ConstVol::VPSSMgr_ConstVol(const VPSSMgr_ConstVol& right) :
    VPSSMgr(right.m_vptp_ptr, right.m_spthermo)
{
    m_useTmpRefStateStorage = true;
    m_useTmpStandardStateStorage = true;
    *this = right;
}


VPSSMgr_ConstVol& VPSSMgr_ConstVol::operator=(const VPSSMgr_ConstVol& b)
{
    if (&b == this) {
        return *this;
    }
    VPSSMgr::operator=(b);
    return *this;
}

VPSSMgr* VPSSMgr_ConstVol::duplMyselfAsVPSSMgr() const
{
    return new VPSSMgr_ConstVol(*this);
}

/*
 * Note, this is equal to the reference state entropies
 * due to the zero volume expansivity:
 * i.e., (dS/dp)_T = (dV/dT)_P = 0.0
 */
void VPSSMgr_ConstVol::_updateStandardStateThermo()
{

    doublereal del_pRT = (m_plast - m_p0) / (GasConstant * m_tlast);

    for (size_t k = 0; k < m_kk; k++) {
        m_hss_RT[k]  = m_h0_RT[k] + del_pRT * m_Vss[k];
        m_cpss_R[k]  = m_cp0_R[k];
        m_sss_R[k]   = m_s0_R[k];
        m_gss_RT[k]  = m_hss_RT[k] - m_sss_R[k];
        // m_Vss[k] constant
    }
}

void VPSSMgr_ConstVol::getGibbs_RT_ref(doublereal* grt) const
{
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    } else {
        throw CanteraError("VPSSMgr_ConstVol::getGibbs_RT_ref",
                           "unimplemented without m_useTmpRefStateStorage");
    }
}

void VPSSMgr_ConstVol::getStandardVolumes_ref(doublereal* vol) const
{
    if (m_useTmpStandardStateStorage) {
        std::copy(m_Vss.begin(), m_Vss.end(), vol);
    } else {
        throw CanteraError("VPSSMgr_ConstVol::getStandardVolumes_ref",
                           "unimplemented without m_useTmpRefStateStorage");
    }
}

void VPSSMgr_ConstVol::initThermo()
{
    VPSSMgr::initThermo();
}

void
VPSSMgr_ConstVol::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    VPSSMgr::initThermoXML(phaseNode, id);

    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &phaseNode.root());

    for (size_t k = 0; k < m_kk; k++) {
        const XML_Node* s = speciesDB->findByAttr("name", m_vptp_ptr->speciesName(k));
        if (!s) {
            throw CanteraError("VPSSMgr_ConstVol::initThermoXML",
                               "no species Node for species " + m_vptp_ptr->speciesName(k));
        }
        const XML_Node* ss = s->findByName("standardState");
        if (!ss) {
            throw CanteraError("VPSSMgr_ConstVol::initThermoXML",
                               "no standardState Node for species " + s->attrib("name"));
        }
        std::string model = ss->attrib("model");
        if (model != "constant_incompressible" && model != "constantVolume") {
            throw CanteraError("VPSSMgr_ConstVol::initThermoXML",
                               "standardState model for species isn't constant_incompressible: " + s->attrib("name"));
        }
        m_Vss[k] = getFloat(*ss, "molarVolume", "toSI");
    }
}

PDSS*
VPSSMgr_ConstVol::createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr)
{
    const XML_Node* ss = speciesNode.findByName("standardState");
    if (!ss) {
        throw CanteraError("VPSSMgr_ConstVol::createInstallPDSS",
                           "no standardState Node for species " + speciesNode["name"]);
    }
    std::string model = ss->attrib("model");
    if (model != "constant_incompressible" && model != "constantVolume") {
        throw CanteraError("VPSSMgr_ConstVol::createInstallPDSS",
                           "standardState model for species isn't "
                           "constant_incompressible: " + speciesNode["name"]);
    }
    if (m_Vss.size() < k+1) {
        m_Vss.resize(k+1, 0.0);
    }
    m_Vss[k] = getFloat(*ss, "molarVolume", "toSI");

    installSTSpecies(k, speciesNode, phaseNode_ptr);
    return new PDSS_ConstVol(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
}

PDSS_enumType VPSSMgr_ConstVol::reportPDSSType(int k) const
{
    return cPDSS_CONSTVOL;
}

VPSSMgr_enumType VPSSMgr_ConstVol::reportVPSSMgrType() const
{
    return  cVPSSMGR_CONSTVOL;
}
}
