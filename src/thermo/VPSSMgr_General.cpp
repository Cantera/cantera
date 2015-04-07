/**
 *  @file VPSSMgr_General.cpp
 *  Definition file for a derived class that handles the calculation
 *  of standard state thermo properties for
 *  a set of species belonging to a single phase in a completely general
 *  but slow way (see \ref thermoprops and
 *  class \link Cantera::VPSSMgr_General VPSSMgr_General\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/VPSSMgr_General.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/base/xml.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/PDSS_HKFT.h"
#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"

using namespace std;

namespace Cantera
{

VPSSMgr_General::VPSSMgr_General(VPStandardStateTP* vp_ptr,
                                 SpeciesThermo* spth) :
    VPSSMgr(vp_ptr, spth)
{
    // Might want to do something other than holding this true.
    //    However, for the sake of getting this all up and running,
    //    will not go there for now.
    m_useTmpStandardStateStorage = true;
    m_useTmpRefStateStorage = true;
}

VPSSMgr_General::VPSSMgr_General(const VPSSMgr_General& right) :
    VPSSMgr(right.m_vptp_ptr, right.m_spthermo)
{
    m_useTmpStandardStateStorage = true;
    m_useTmpRefStateStorage = true;
    *this = right;
}

VPSSMgr_General& VPSSMgr_General::operator=(const VPSSMgr_General& b)
{
    if (&b == this) {
        return *this;
    }
    VPSSMgr::operator=(b);
    /*
     *  Must fill in the shallow pointers. These must have already been transfered
     *  and stored in the owning VPStandardStateTP class.  Note we are aware that at this point
     *  m_vptr_ptr may refer back to the wrong ThermoPhase object. However, the shallow copy
     *  performed here is consistent with the assignment operator's general functionality.
     */
    m_PDSS_ptrs.resize(m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        m_PDSS_ptrs[k] = m_vptp_ptr->providePDSS(k);
    }
    return *this;
}

VPSSMgr* VPSSMgr_General::duplMyselfAsVPSSMgr() const
{
    return new VPSSMgr_General(*this);
}

void VPSSMgr_General::initAllPtrs(VPStandardStateTP* vp_ptr, SpeciesThermo* sp_ptr)
{
    VPSSMgr::initAllPtrs(vp_ptr, sp_ptr);
    /*
     *  Must fill in the shallow pointers. These must have already been transfered
     *  and stored in the owning VPStandardStateTP class.
     */
    m_PDSS_ptrs.resize(m_kk);
    for (size_t k = 0; k < m_kk; k++) {
        m_PDSS_ptrs[k] = m_vptp_ptr->providePDSS(k);
    }
}

void VPSSMgr_General::_updateRefStateThermo() const
{
    if (m_useTmpRefStateStorage) {
        for (size_t k = 0; k < m_kk; k++) {
            PDSS* kPDSS = m_PDSS_ptrs[k];
            kPDSS->setState_TP(m_tlast, m_plast);
            m_h0_RT[k] = kPDSS->enthalpy_RT_ref();
            m_s0_R[k]  = kPDSS->entropy_R_ref();
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
            m_cp0_R[k] = kPDSS->cp_R_ref();
            m_V0[k]    = kPDSS->molarVolume_ref();
        }
    }
}

void VPSSMgr_General::_updateStandardStateThermo()
{
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_ptrs[k];
        kPDSS->setState_TP(m_tlast, m_plast);
        m_hss_RT[k] = kPDSS->enthalpy_RT();
        m_sss_R[k]  = kPDSS->entropy_R();
        m_gss_RT[k] = m_hss_RT[k] - m_sss_R[k];
        m_cpss_R[k] = kPDSS->cp_R();
        m_Vss[k]    = kPDSS->molarVolume();
    }
}

void VPSSMgr_General::initThermo()
{
    initLengths();
}

void VPSSMgr_General::getGibbs_ref(doublereal* g) const
{
    doublereal _rt = GasConstant * m_tlast;
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), g);
        scale(g, g+m_kk, g, _rt);
    } else {
        for (size_t k = 0; k < m_kk; k++) {
            PDSS* kPDSS = m_PDSS_ptrs[k];
            kPDSS->setState_TP(m_tlast, m_plast);
            double h0_RT = kPDSS->enthalpy_RT_ref();
            double s0_R  = kPDSS->entropy_R_ref();
            g[k] = _rt * (h0_RT - s0_R);
        }
    }
}

void
VPSSMgr_General::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    VPSSMgr::initThermoXML(phaseNode, id);
}

PDSS*
VPSSMgr_General::returnPDSS_ptr(size_t k, const XML_Node& speciesNode,
                                const XML_Node* const phaseNode_ptr, bool& doST)
{
    PDSS* kPDSS = 0;
    doST = true;
    GeneralSpeciesThermo* genSpthermo = dynamic_cast<GeneralSpeciesThermo*>(m_spthermo);


    const XML_Node* const ss = speciesNode.findByName("standardState");
    if (!ss) {
        VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);
        kPDSS = new PDSS_IdealGas(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
        return kPDSS;
    }
    std::string model = (*ss)["model"];
    if (model == "constant_incompressible") {
        VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);
        kPDSS = new PDSS_ConstVol(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
        if (!kPDSS) {
            throw CanteraError("VPSSMgr_General::returnPDSS_ptr", "new PDSS_ConstVol failed");
        }
    } else if (model == "waterIAPWS" || model == "waterPDSS") {
        // VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);
        kPDSS = new PDSS_Water(m_vptp_ptr, 0);
        if (!genSpthermo) {
            throw CanteraError("VPSSMgr_General::returnPDSS_ptr",
                               "failed dynamic cast");
        }
        genSpthermo->installPDSShandler(k, kPDSS, this);
        m_useTmpRefStateStorage = false;
    } else if (model == "HKFT") {
        doST = false;
        kPDSS = new PDSS_HKFT(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
        if (!genSpthermo) {
            throw CanteraError("VPSSMgr_General::returnPDSS_ptr",
                               "failed dynamic cast");
        }
        genSpthermo->installPDSShandler(k, kPDSS, this);

    } else if (model == "IonFromNeutral") {
        if (!genSpthermo) {
            throw CanteraError("VPSSMgr_General::returnPDSS_ptr",
                               "failed dynamic cast");
        }
        doST = false;
        kPDSS = new PDSS_IonsFromNeutral(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
        if (!kPDSS) {
            throw CanteraError("VPSSMgr_General::returnPDSS_ptr",
                               "new PDSS_IonsFromNeutral failed");
        }
        genSpthermo->installPDSShandler(k, kPDSS, this);

    } else if (model == "constant" || model == "temperature_polynomial" || model == "density_temperature_polynomial") {
        VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);
        kPDSS = new PDSS_SSVol(m_vptp_ptr, k, speciesNode, *phaseNode_ptr, true);
        if (!kPDSS) {
            throw CanteraError("VPSSMgr_General::returnPDSS_ptr", "new PDSS_SSVol failed");
        }
    } else {
        throw CanteraError("VPSSMgr_General::returnPDSS_ptr",
                           "unknown standard state formulation: " + model);
    }
    return kPDSS;
}

PDSS*
VPSSMgr_General::createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                   const XML_Node* const phaseNode_ptr)
{
    bool doST;
    PDSS* kPDSS = returnPDSS_ptr(k, speciesNode, phaseNode_ptr, doST);
    // VPSSMgr::installSTSpecies(k, speciesNode, phaseNode_ptr);
    if (m_PDSS_ptrs.size() < k+1) {
        m_PDSS_ptrs.resize(k+1, 0);
    }
    m_PDSS_ptrs[k] = kPDSS;
    if ((k+1) >= m_kk) {
        m_kk = k+1;
    }

    doublereal minTemp = kPDSS->minTemp();
    if (minTemp > m_minTemp) {
        m_minTemp = minTemp;
    }

    doublereal maxTemp = kPDSS->maxTemp();
    if (maxTemp < m_maxTemp) {
        m_maxTemp = maxTemp;
    }

    doublereal p0 = kPDSS->refPressure();
    if (k == 0) {
        m_p0 = p0;
    }
    return kPDSS;
}

PDSS_enumType VPSSMgr_General::reportPDSSType(int k) const
{
    PDSS* kPDSS = m_PDSS_ptrs[k];
    return  kPDSS->reportPDSSType();
}

VPSSMgr_enumType  VPSSMgr_General::reportVPSSMgrType() const
{
    return cVPSSMGR_GENERAL;
}
}
