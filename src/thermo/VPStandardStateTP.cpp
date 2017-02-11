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
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/PDSS_HKFT.h"
#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/base/utilities.h"
#include "cantera/base/ctml.h"

using namespace std;

namespace Cantera
{

VPStandardStateTP::VPStandardStateTP() :
    m_Pcurrent(OneAtm),
    m_Tlast_ss(-1.0),
    m_Plast_ss(-1.0),
    m_useTmpRefStateStorage(true)
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
    std::copy(m_hss_RT.begin(), m_hss_RT.end(), hrt);
}

void VPStandardStateTP::getEntropy_R(doublereal* sr) const
{
    updateStandardStateThermo();
    std::copy(m_sss_R.begin(), m_sss_R.end(), sr);
}

void VPStandardStateTP::getGibbs_RT(doublereal* grt) const
{
    updateStandardStateThermo();
    std::copy(m_gss_RT.begin(), m_gss_RT.end(), grt);
}

void VPStandardStateTP::getPureGibbs(doublereal* g) const
{
    updateStandardStateThermo();
    std::copy(m_gss_RT.begin(), m_gss_RT.end(), g);
    scale(g, g+m_kk, g, RT());
}

void VPStandardStateTP::getIntEnergy_RT(doublereal* urt) const
{
    updateStandardStateThermo();
    std::copy(m_hss_RT.begin(), m_hss_RT.end(), urt);
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] -= m_Plast_ss / RT() * m_Vss[k];
    }
}

void VPStandardStateTP::getCp_R(doublereal* cpr) const
{
    updateStandardStateThermo();
    std::copy(m_cpss_R.begin(), m_cpss_R.end(), cpr);
}

void VPStandardStateTP::getStandardVolumes(doublereal* vol) const
{
    updateStandardStateThermo();
    std::copy(m_Vss.begin(), m_Vss.end(), vol);
}
const vector_fp& VPStandardStateTP::getStandardVolumes() const
{
    updateStandardStateThermo();
    return m_Vss;
}

// ----- Thermodynamic Values for the Species Reference States ----

void VPStandardStateTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    updateStandardStateThermo();
    if (m_useTmpRefStateStorage) {
        std::copy(m_h0_RT.begin(), m_h0_RT.end(), hrt);
    } else {
        throw NotImplementedError("VPStandardStateTP::getEnthalpy_RT_ref");
    }
}

void VPStandardStateTP::getGibbs_RT_ref(doublereal* grt) const
{
    updateStandardStateThermo();
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), grt);
    } else {
        throw NotImplementedError("VPStandardStateTP::getGibbs_RT_ref");
    }
}

void VPStandardStateTP::getGibbs_ref(doublereal* g) const
{
    updateStandardStateThermo();
    if (m_useTmpRefStateStorage) {
        std::copy(m_g0_RT.begin(), m_g0_RT.end(), g);
        scale(g, g+m_kk, g, RT());
    } else {
        for (size_t k = 0; k < m_kk; k++) {
            PDSS* kPDSS = m_PDSS_storage[k].get();
            kPDSS->setState_TP(m_tlast, m_Plast_ss);
            double h0_RT = kPDSS->enthalpy_RT_ref();
            double s0_R = kPDSS->entropy_R_ref();
            g[k] = RT() * (h0_RT - s0_R);
        }
    }
}

const vector_fp& VPStandardStateTP::Gibbs_RT_ref() const
{
    updateStandardStateThermo();
    if (m_useTmpRefStateStorage) {
        return m_g0_RT;
    } else {
        throw NotImplementedError("VPStandardStateTP::getGibbs_RT_ref");
    }
}

void VPStandardStateTP::getEntropy_R_ref(doublereal* sr) const
{
    updateStandardStateThermo();
    if (m_useTmpRefStateStorage) {
        std::copy(m_s0_R.begin(), m_s0_R.end(), sr);
    } else {
        throw NotImplementedError("VPStandardStateTP::getEntropy_R_ref");
    }
}

void VPStandardStateTP::getCp_R_ref(doublereal* cpr) const
{
    updateStandardStateThermo();
    if (m_useTmpRefStateStorage) {
        std::copy(m_cp0_R.begin(), m_cp0_R.end(), cpr);
    } else {
        throw NotImplementedError("VPStandardStateTP::getCp_R_ref");
    }
}

void VPStandardStateTP::getStandardVolumes_ref(doublereal* vol) const
{
    updateStandardStateThermo();
    std::copy(m_Vss.begin(), m_Vss.end(), vol);
}

void VPStandardStateTP::initThermo()
{
    ThermoPhase::initThermo();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k].get();
        if (kPDSS) {
            kPDSS->initThermo();
        }
    }
}

bool VPStandardStateTP::addSpecies(shared_ptr<Species> spec)
{
    // Specifically skip ThermoPhase::addSpecies since the Species object
    // doesn't have an associated SpeciesThermoInterpType object
    bool added = Phase::addSpecies(spec);
    if (!added) {
        return false;
    }
    m_h0_RT.push_back(0.0);
    m_cp0_R.push_back(0.0);
    m_g0_RT.push_back(0.0);
    m_s0_R.push_back(0.0);
    m_V0.push_back(0.0);
    m_hss_RT.push_back(0.0);
    m_cpss_R.push_back(0.0);
    m_gss_RT.push_back(0.0);
    m_sss_R.push_back(0.0);
    m_Vss.push_back(0.0);
    return true;
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
                                          const XML_Node* phaseNode)
{
    if (m_PDSS_storage.size() < k+1) {
        m_PDSS_storage.resize(k+1);
    }
    PDSS* kPDSS = nullptr;
    bool use_STITbyPDSS;

    const XML_Node* const ss = s.findByName("standardState");
    if (!ss) {
        use_STITbyPDSS = false;
        kPDSS = new PDSS_IdealGas(this, k, s, *phaseNode, true);
    } else {
        std::string model = ss->attrib("model");
        if (model == "constant_incompressible") {
            kPDSS = new PDSS_ConstVol(this, k, s, *phaseNode, true);
            use_STITbyPDSS = false;
        } else if (model == "waterIAPWS" || model == "waterPDSS") {
            kPDSS = new PDSS_Water(this, 0);
            use_STITbyPDSS = true;
            m_useTmpRefStateStorage = false;
        } else if (model == "HKFT") {
            kPDSS = new PDSS_HKFT(this, k, s, *phaseNode, true);
            use_STITbyPDSS = true;
        } else if (model == "IonFromNeutral") {
            kPDSS = new PDSS_IonsFromNeutral(this, k, s, *phaseNode, true);
            use_STITbyPDSS = true;
        } else if (model == "constant" || model == "temperature_polynomial" || model == "density_temperature_polynomial") {
            kPDSS = new PDSS_SSVol(this, k, s, *phaseNode, true);
            use_STITbyPDSS = false;
        } else {
            throw CanteraError("VPStandardStateTP::createInstallPDSS",
                               "unknown standard state formulation: " + model);
        }
    }

    if (use_STITbyPDSS) {
        auto stit = make_shared<STITbyPDSS>(kPDSS);
        m_spthermo->install_STIT(k, stit);
    } else {
        shared_ptr<SpeciesThermoInterpType> stit(
            newSpeciesThermoInterpType(s.child("thermo")));
        stit->validate(s["name"]);
        m_spthermo->install_STIT(k, stit);
    }

    m_PDSS_storage[k].reset(kPDSS);
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
    ThermoPhase::initThermoXML(phaseNode, id);
}

void VPStandardStateTP::_updateStandardStateThermo() const
{
    double Tnow = temperature();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k].get();
        kPDSS->setState_TP(Tnow, m_Pcurrent);
        // reference state thermo
        if (Tnow != m_tlast && m_useTmpRefStateStorage) {
            m_h0_RT[k] = kPDSS->enthalpy_RT_ref();
            m_s0_R[k] = kPDSS->entropy_R_ref();
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
            m_cp0_R[k] = kPDSS->cp_R_ref();
            m_V0[k] = kPDSS->molarVolume_ref();
        }
        // standard state thermo
        m_hss_RT[k] = kPDSS->enthalpy_RT();
        m_sss_R[k] = kPDSS->entropy_R();
        m_gss_RT[k] = m_hss_RT[k] - m_sss_R[k];
        m_cpss_R[k] = kPDSS->cp_R();
        m_Vss[k] = kPDSS->molarVolume();
    }
    m_Plast_ss = m_Pcurrent;
    m_Tlast_ss = Tnow;
    m_tlast = Tnow;
}

void VPStandardStateTP::updateStandardStateThermo() const
{
    double Tnow = temperature();
    if (Tnow != m_Tlast_ss || Tnow != m_tlast || m_Pcurrent != m_Plast_ss) {
        _updateStandardStateThermo();
    }
}
}
