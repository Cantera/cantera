/**
 *  @file VPStandardStateTP.cpp
 * Definition file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see @ref thermoprops and
 * class @link Cantera::VPStandardStateTP VPStandardStateTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/VPStandardStateTP.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/utilities.h"
#include "cantera/base/global.h"

namespace Cantera
{

VPStandardStateTP::VPStandardStateTP()
{
    // Defined in .cpp to limit dependence on PDSS.h via vector<unique_ptr<PDSS>>
}

VPStandardStateTP::~VPStandardStateTP()
{
    // Defined in .cpp to limit dependence on PDSS.h
}

int VPStandardStateTP::standardStateConvention() const
{
    return cSS_CONVENTION_VPSS;
}

// ----- Thermodynamic Values for the Species Standard States States ----

void VPStandardStateTP::getStandardChemPotentials(span<double> g) const
{
    getGibbs_RT(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT();
    }
}

void VPStandardStateTP::getEnthalpy_RT(span<double> hrt) const
{
    checkArraySize("VPStandardStateTP::getEnthalpy_RT", hrt.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_hss_RT.begin(), m_hss_RT.end(), hrt.begin());
}

void VPStandardStateTP::getEntropy_R(span<double> sr) const
{
    checkArraySize("VPStandardStateTP::getEntropy_R", sr.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_sss_R.begin(), m_sss_R.end(), sr.begin());
}

void VPStandardStateTP::getGibbs_RT(span<double> grt) const
{
    checkArraySize("VPStandardStateTP::getGibbs_RT", grt.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_gss_RT.begin(), m_gss_RT.end(), grt.begin());
}

void VPStandardStateTP::getIntEnergy_RT(span<double> urt) const
{
    checkArraySize("VPStandardStateTP::getIntEnergy_RT", urt.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_hss_RT.begin(), m_hss_RT.end(), urt.begin());
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] -= m_Plast_ss / RT() * m_Vss[k];
    }
}

void VPStandardStateTP::getCp_R(span<double> cpr) const
{
    checkArraySize("VPStandardStateTP::getCp_R", cpr.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_cpss_R.begin(), m_cpss_R.end(), cpr.begin());
}

void VPStandardStateTP::getStandardVolumes(span<double> vol) const
{
    checkArraySize("VPStandardStateTP::getStandardVolumes", vol.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_Vss.begin(), m_Vss.end(), vol.begin());
}
span<const double> VPStandardStateTP::getStandardVolumes() const
{
    updateStandardStateThermo();
    return m_Vss;
}

// ----- Thermodynamic Values for the Species Reference States ----

void VPStandardStateTP::getEnthalpy_RT_ref(span<double> hrt) const
{
    checkArraySize("VPStandardStateTP::getEnthalpy_RT_ref", hrt.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_h0_RT.begin(), m_h0_RT.end(), hrt.begin());
}

void VPStandardStateTP::getGibbs_RT_ref(span<double> grt) const
{
    checkArraySize("VPStandardStateTP::getGibbs_RT_ref", grt.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_g0_RT.begin(), m_g0_RT.end(), grt.begin());
}

void VPStandardStateTP::getGibbs_ref(span<double> g) const
{
    checkArraySize("VPStandardStateTP::getGibbs_ref", g.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_g0_RT.begin(), m_g0_RT.end(), g.begin());
    scale(g.begin(), g.end(), g.begin(), RT());
}

void VPStandardStateTP::getEntropy_R_ref(span<double> sr) const
{
    checkArraySize("VPStandardStateTP::getEntropy_R_ref", sr.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_s0_R.begin(), m_s0_R.end(), sr.begin());
}

void VPStandardStateTP::getCp_R_ref(span<double> cpr) const
{
    checkArraySize("VPStandardStateTP::getCp_R_ref", cpr.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_cp0_R.begin(), m_cp0_R.end(), cpr.begin());
}

void VPStandardStateTP::getStandardVolumes_ref(span<double> vol) const
{
    checkArraySize("VPStandardStateTP::getStandardVolumes_ref", vol.size(), m_kk);
    updateStandardStateThermo();
    std::copy(m_Vss.begin(), m_Vss.end(), vol.begin());
}

void VPStandardStateTP::initThermo()
{
    ThermoPhase::initThermo();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k].get();
        if (kPDSS == 0) {
            throw CanteraError("VPStandardStateTP::initThermo",
                "No PDSS object for species {}", k);
        }
        kPDSS->initThermo();
    }
}

void VPStandardStateTP::getSpeciesParameters(const string& name,
                                             AnyMap& speciesNode) const
{
    AnyMap eos;
    providePDSS(speciesIndex(name, true))->getParameters(eos);
    speciesNode["equation-of-state"].getMapWhere(
        "model", eos.getString("model", ""), true) = std::move(eos);
}

bool VPStandardStateTP::addSpecies(shared_ptr<Species> spec)
{
    // Specifically skip ThermoPhase::addSpecies since the Species object
    // doesn't have an associated SpeciesThermoInterpType object
    bool added = Phase::addSpecies(spec);
    if (!added) {
        return false;
    }

    // VPStandardState does not use m_spthermo - install a dummy object
    m_spthermo.install_STIT(m_kk-1, make_shared<SpeciesThermoInterpType>());
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

void VPStandardStateTP::setTemperature(const double temp)
{
    setState_TP(temp, m_Pcurrent);
    updateStandardStateThermo();
}

void VPStandardStateTP::setPressure(double p)
{
    setState_TP(temperature(), p);
    updateStandardStateThermo();
}

void VPStandardStateTP::calcDensity()
{
    getPartialMolarVolumes(m_workS);
    double dd = meanMolecularWeight() / mean_X(m_workS);
    Phase::assignDensity(dd);
}

void VPStandardStateTP::setState_TP(double t, double pres)
{
    if (pres < 0) {
        throw CanteraError("VPStandardStateTP::setState_TP",
            "Pressure must be positive. Specified state was T = {}, P = {}.", t, pres);
    }
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

void VPStandardStateTP::installPDSS(size_t k, unique_ptr<PDSS>&& pdss)
{
    pdss->setParent(this, k);
    pdss->setMolecularWeight(molecularWeight(k));
    Species& spec = *species(k);
    if (spec.thermo) {
        pdss->setReferenceThermo(spec.thermo);
        spec.thermo->validate(spec.name);
    }
    m_minTemp = std::max(m_minTemp, pdss->minTemp());
    m_maxTemp = std::min(m_maxTemp, pdss->maxTemp());

    if (m_PDSS_storage.size() < k+1) {
        m_PDSS_storage.resize(k+1);
    }
    m_PDSS_storage[k].swap(pdss);
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

void VPStandardStateTP::_updateStandardStateThermo() const
{
    double Tnow = temperature();
    for (size_t k = 0; k < m_kk; k++) {
        PDSS* kPDSS = m_PDSS_storage[k].get();
        kPDSS->setState_TP(Tnow, m_Pcurrent);
        // reference state thermo
        if (Tnow != m_tlast) {
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

double VPStandardStateTP::minTemp(size_t k) const
{
    if (k == npos) {
        return m_minTemp;
    } else {
        return m_PDSS_storage.at(k)->minTemp();
    }
}

double VPStandardStateTP::maxTemp(size_t k) const
{
    if (k == npos) {
        return m_maxTemp;
    } else {
        return m_PDSS_storage.at(k)->maxTemp();
    }
}

}
