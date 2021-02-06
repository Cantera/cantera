/**
 *  @file SurfPhase.cpp
 *  Definitions for a simple thermodynamic model of a surface phase
 *  derived from ThermoPhase, assuming an ideal solution model
 *  (see \ref thermoprops and class
 *  \link Cantera::SurfPhase SurfPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/EdgePhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
SurfPhase::SurfPhase(doublereal n0):
    m_press(OneAtm)
{
    setSiteDensity(n0);
    setNDim(2);
}

SurfPhase::SurfPhase(const std::string& infile, const std::string& id_) :
    m_press(OneAtm)
{
    initThermoFile(infile, id_);
}

SurfPhase::SurfPhase(XML_Node& xmlphase) :
    m_press(OneAtm)
{
    importPhase(xmlphase, this);
}

doublereal SurfPhase::enthalpy_mole() const
{
    if (m_n0 <= 0.0) {
        return 0.0;
    }
    _updateThermo();
    return mean_X(m_h0);
}

doublereal SurfPhase::intEnergy_mole() const
{
    return enthalpy_mole();
}

doublereal SurfPhase::entropy_mole() const
{
    _updateThermo();
    doublereal s = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        s += moleFraction(k) * (m_s0[k] -
            GasConstant * log(std::max(concentration(k) * size(k)/m_n0, SmallNumber)));
    }
    return s;
}

doublereal SurfPhase::cp_mole() const
{
    _updateThermo();
    return mean_X(m_cp0);
}

doublereal SurfPhase::cv_mole() const
{
    return cp_mole();
}

void SurfPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT();
    }
}

void SurfPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= GasConstant;
    }
}

void SurfPhase::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

// HKM 9/1/11  The partial molar volumes returned here are really partial molar areas.
//             Partial molar volumes for this phase should actually be equal to zero.
void SurfPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

void SurfPhase::getStandardChemPotentials(doublereal* mu0) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu0);
}

void SurfPhase::getChemPotentials(doublereal* mu) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), mu);
    getActivityConcentrations(m_work.data());
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += RT() * (log(m_work[k]) - logStandardConc(k));
    }
}

void SurfPhase::getActivityConcentrations(doublereal* c) const
{
    getConcentrations(c);
}

doublereal SurfPhase::standardConcentration(size_t k) const
{
    return m_n0/size(k);
}

doublereal SurfPhase::logStandardConc(size_t k) const
{
    return m_logn0 - m_logsize[k];
}

void SurfPhase::setParameters(int n, doublereal* const c)
{
    if (n != 1) {
        throw CanteraError("SurfPhase::setParameters",
                           "Bad value for number of parameter");
    }
    setSiteDensity(c[0]);
}

void SurfPhase::getPureGibbs(doublereal* g) const
{
    _updateThermo();
    copy(m_mu0.begin(), m_mu0.end(), g);
}

void SurfPhase::getGibbs_RT(doublereal* grt) const
{
    _updateThermo();
    scale(m_mu0.begin(), m_mu0.end(), grt, 1.0/RT());
}

void SurfPhase::getEnthalpy_RT(doublereal* hrt) const
{
    _updateThermo();
    scale(m_h0.begin(), m_h0.end(), hrt, 1.0/RT());
}

void SurfPhase::getEntropy_R(doublereal* sr) const
{
    _updateThermo();
    scale(m_s0.begin(), m_s0.end(), sr, 1.0/GasConstant);
}

void SurfPhase::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    scale(m_cp0.begin(), m_cp0.end(), cpr, 1.0/GasConstant);
}

void SurfPhase::getStandardVolumes(doublereal* vol) const
{
    _updateThermo();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = 1.0/standardConcentration(k);
    }
}

void SurfPhase::getGibbs_RT_ref(doublereal* grt) const
{
    getGibbs_RT(grt);
}

void SurfPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    getEnthalpy_RT(hrt);
}

void SurfPhase::getEntropy_R_ref(doublereal* sr) const
{
    getEntropy_R(sr);
}

void SurfPhase::getCp_R_ref(doublereal* cprt) const
{
    getCp_R(cprt);
}

bool SurfPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        m_h0.push_back(0.0);
        m_s0.push_back(0.0);
        m_cp0.push_back(0.0);
        m_mu0.push_back(0.0);
        m_work.push_back(0.0);
        m_speciesSize.push_back(spec->size);
        m_logsize.push_back(log(spec->size));
        if (m_kk == 1) {
            vector_fp cov{1.0};
            setCoverages(cov.data());
        }
    }
    return added;
}

void SurfPhase::setSiteDensity(doublereal n0)
{
    if (n0 <= 0.0) {
        throw CanteraError("SurfPhase::setSiteDensity",
                           "Site density must be positive. Got {}", n0);
    }
    m_n0 = n0;
    m_logn0 = log(m_n0);
}

void SurfPhase::setCoverages(const doublereal* theta)
{
    double sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        sum += theta[k];
    }
    if (sum <= 0.0) {
        throw CanteraError("SurfPhase::setCoverages",
                           "Sum of Coverage fractions is zero or negative");
    }
    for (size_t k = 0; k < m_kk; k++) {
        m_work[k] = m_n0*theta[k]/(sum*size(k));
    }
    // Call the Phase:: class function setConcentrations.
    setConcentrations(m_work.data());
}

void SurfPhase::setCoveragesNoNorm(const doublereal* theta)
{
    for (size_t k = 0; k < m_kk; k++) {
        m_work[k] = m_n0*theta[k]/size(k);
    }
    setConcentrationsNoNorm(m_work.data());
}

void SurfPhase::getCoverages(doublereal* theta) const
{
    getConcentrations(theta);
    for (size_t k = 0; k < m_kk; k++) {
        theta[k] *= size(k)/m_n0;
    }
}

void SurfPhase::setCoveragesByName(const std::string& cov)
{
    setCoveragesByName(parseCompString(cov, speciesNames()));
}

void SurfPhase::setCoveragesByName(const compositionMap& cov)
{
    vector_fp cv(m_kk, 0.0);
    bool ifound = false;
    for (size_t k = 0; k < m_kk; k++) {
        double c = getValue(cov, speciesName(k), 0.0);
        if (c > 0.0) {
            ifound = true;
            cv[k] = c;
        }
    }
    if (!ifound) {
        throw CanteraError("SurfPhase::setCoveragesByName",
                           "Input coverages are all zero or negative");
    }
    setCoverages(cv.data());
}

void SurfPhase::setState(const AnyMap& state) {
    if (state.hasKey("coverages")) {
        if (state["coverages"].is<string>()) {
            setCoveragesByName(state["coverages"].asString());
        } else {
            setCoveragesByName(state["coverages"].asMap<double>());
        }
    }
    ThermoPhase::setState(state);
}

void SurfPhase::_updateThermo(bool force) const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow || force) {
        m_spthermo.update(tnow, m_cp0.data(), m_h0.data(), m_s0.data());
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_h0[k] *= GasConstant * tnow;
            m_s0[k] *= GasConstant;
            m_cp0[k] *= GasConstant;
            m_mu0[k] = m_h0[k] - tnow*m_s0[k];
        }
        m_tlast = tnow;
    }
}

void SurfPhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","Surface");
    doublereal n = getFloat(eosdata, "site_density", "toSI");
    setSiteDensity(n);
}

void SurfPhase::initThermo()
{
    if (m_input.hasKey("site-density")) {
        // Units are kmol/m^2 for surface phases or kmol/m for edge phases
        setSiteDensity(m_input.convert("site-density",
            Units(1.0, 0, -static_cast<double>(m_ndim), 0, 0, 0, 1)));
    }
}

void SurfPhase::setStateFromXML(const XML_Node& state)
{
    double t;
    if (getOptionalFloat(state, "temperature", t, "temperature")) {
        setTemperature(t);
    }

    if (state.hasChild("coverages")) {
        string comp = getChildValue(state,"coverages");
        setCoveragesByName(comp);
    }
}

EdgePhase::EdgePhase(doublereal n0) : SurfPhase(n0)
{
    setNDim(1);
}

void EdgePhase::setParametersFromXML(const XML_Node& eosdata)
{
    eosdata._require("model","Edge");
    doublereal n = getFloat(eosdata, "site_density", "toSI");
    setSiteDensity(n);
}

}
