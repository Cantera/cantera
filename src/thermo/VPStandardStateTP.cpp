/**
 *  @file VPStandardStateTP.cpp
 * Definition file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::VPStandardStateTP VPStandardStateTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/PDSS.h"

using namespace std;

namespace Cantera
{

VPStandardStateTP::VPStandardStateTP() :
    m_Pcurrent(OneAtm),
    m_Tlast_ss(-1.0),
    m_Plast_ss(-1.0)
{
}

int VPStandardStateTP::standardStateConvention() const
{
    return cSS_CONVENTION_VPSS;
}

void VPStandardStateTP::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= 1.0 / RT();
    }
}

// ----- Thermodynamic Values for the Species Standard States States ----

void VPStandardStateTP::getStandardChemPotentials(doublereal* g) const
{
    getGibbs_RT(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT();
    }
}

void VPStandardStateTP::getEnthalpy_RT(doublereal* hrt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEnthalpy_RT(hrt);
}

void VPStandardStateTP::getEntropy_R(doublereal* srt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEntropy_R(srt);
}

void VPStandardStateTP::getGibbs_RT(doublereal* grt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getGibbs_RT(grt);
}

void VPStandardStateTP::getPureGibbs(doublereal* g) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getStandardChemPotentials(g);
}

void VPStandardStateTP::getIntEnergy_RT(doublereal* urt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getIntEnergy_RT(urt);
}

void VPStandardStateTP::getCp_R(doublereal* cpr) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getCp_R(cpr);
}

void VPStandardStateTP::getStandardVolumes(doublereal* vol) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getStandardVolumes(vol);
}
const vector_fp& VPStandardStateTP::getStandardVolumes() const
{
    updateStandardStateThermo();
    return m_VPSS_ptr->getStandardVolumes();
}

// ----- Thermodynamic Values for the Species Reference States ----

void VPStandardStateTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEnthalpy_RT_ref(hrt);
}

void VPStandardStateTP::getGibbs_RT_ref(doublereal* grt) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getGibbs_RT_ref(grt);
}

void VPStandardStateTP::getGibbs_ref(doublereal* g) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getGibbs_ref(g);
}

const vector_fp& VPStandardStateTP::Gibbs_RT_ref() const
{
    updateStandardStateThermo();
    return m_VPSS_ptr->Gibbs_RT_ref();
}

void VPStandardStateTP::getEntropy_R_ref(doublereal* er) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getEntropy_R_ref(er);
}

void VPStandardStateTP::getCp_R_ref(doublereal* cpr) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getCp_R_ref(cpr);
}

void VPStandardStateTP::getStandardVolumes_ref(doublereal* vol) const
{
    updateStandardStateThermo();
    m_VPSS_ptr->getStandardVolumes_ref(vol);
}

void VPStandardStateTP::initThermo()
{
    ThermoPhase::initThermo();
    m_VPSS_ptr->initThermo();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k].get();
        if (kPDSS) {
            kPDSS->initThermo();
        }
    }
}

void VPStandardStateTP::setVPSSMgr(VPSSMgr* vp_ptr)
{
    m_VPSS_ptr.reset(vp_ptr);
}

bool VPStandardStateTP::addSpecies(shared_ptr<Species> spec)
{
    // Specifically skip ThermoPhase::addSpecies since the Species object
    // doesn't have an associated SpeciesThermoInterpType object
    return Phase::addSpecies(spec);
}

void VPStandardStateTP::setTemperature(const doublereal temp)
{
    setState_TP(temp, m_Pcurrent);
    updateStandardStateThermo();
}

void VPStandardStateTP::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
    updateStandardStateThermo();
}

void VPStandardStateTP::calcDensity()
{
    throw NotImplementedError("VPStandardStateTP::calcDensity() called, "
                              "but EOS for phase is not known");
}

void VPStandardStateTP::setState_TP(doublereal t, doublereal pres)
{
    // A pretty tricky algorithm is needed here, due to problems involving
    // standard states of real fluids. For those cases you need to combine the T
    // and P specification for the standard state, or else you may venture into
    // the forbidden zone, especially when nearing the triple point. Therefore,
    // we need to do the standard state thermo calc with the (t, pres) combo.
    Phase::setTemperature(t);
    m_Pcurrent = pres;
    updateStandardStateThermo();

    // Now, we still need to do the calculations for general ThermoPhase
    // objects. So, we switch back to a virtual function call, setTemperature,
    // and setPressure to recalculate stuff for child ThermoPhase objects of the
    // VPStandardStateTP object. At this point, we haven't touched m_tlast or
    // m_plast, so some calculations may still need to be done at the
    // ThermoPhase object level.
    calcDensity();
}

void VPStandardStateTP::createInstallPDSS(size_t k, const XML_Node& s,
                                          const XML_Node* phaseNode_ptr)
{
    if (m_PDSS_storage.size() < k+1) {
        m_PDSS_storage.resize(k+1);
    }
    m_PDSS_storage[k].reset(m_VPSS_ptr->createInstallPDSS(k, s, phaseNode_ptr));
}

PDSS* VPStandardStateTP::providePDSS(size_t k)
{
    return m_PDSS_storage[k].get();
}

const PDSS* VPStandardStateTP::providePDSS(size_t k) const
{
    return m_PDSS_storage[k].get();
}

void VPStandardStateTP::invalidateCache()
{
    ThermoPhase::invalidateCache();
    m_Tlast_ss += 0.0001234;
}

void VPStandardStateTP::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k].get();
        AssertTrace(kPDSS != 0);
        if (kPDSS) {
            kPDSS->initThermoXML(phaseNode, id);
        }
    }
    m_VPSS_ptr->initThermoXML(phaseNode, id);
    ThermoPhase::initThermoXML(phaseNode, id);
}

VPSSMgr* VPStandardStateTP::provideVPSSMgr()
{
    return m_VPSS_ptr.get();
}

void VPStandardStateTP::_updateStandardStateThermo() const
{
    double Tnow = temperature();
    m_Plast_ss = m_Pcurrent;
    m_Tlast_ss = Tnow;
    AssertThrowMsg(m_VPSS_ptr != 0, "VPStandardStateTP::_updateStandardStateThermo()",
                   "Probably indicates that ThermoPhase object wasn't initialized correctly");
    m_VPSS_ptr->setState_TP(Tnow, m_Pcurrent);
}

void VPStandardStateTP::updateStandardStateThermo() const
{
    double Tnow = temperature();
    if (Tnow != m_Tlast_ss || m_Pcurrent != m_Plast_ss) {
        _updateStandardStateThermo();
    }
}
}
