/**
 *  @file IdealGasPhase.cpp
 *   ThermoPhase object for the ideal gas equation of
 * state - workhorse for %Cantera (see \ref thermoprops
 * and class \link Cantera::IdealGasPhase IdealGasPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

IdealGasPhase::IdealGasPhase() :
    m_p0(-1.0)
{
}

IdealGasPhase::IdealGasPhase(const std::string& inputFile, const std::string& id_) :
    m_p0(-1.0)
{
    initThermoFile(inputFile, id_);
}

IdealGasPhase::IdealGasPhase(XML_Node& phaseRef, const std::string& id_) :
    m_p0(-1.0)
{
    importPhase(phaseRef, this);
}

// Molar Thermodynamic Properties of the Solution ------------------

doublereal IdealGasPhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx() - std::log(pressure() / refPressure()));
}

doublereal IdealGasPhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

doublereal IdealGasPhase::cv_mole() const
{
    return cp_mole() - GasConstant;
}

doublereal IdealGasPhase::standardConcentration(size_t k) const
{
    return pressure() / RT();
}

void IdealGasPhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

void IdealGasPhase::getStandardChemPotentials(doublereal* muStar) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), muStar, RT());
    double tmp = log(pressure() / refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp; // add RT*ln(P/P_0)
    }
}

//  Partial Molar Properties of the Solution --------------

void IdealGasPhase::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT() * log(xx);
    }
}

void IdealGasPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    scale(_h.begin(), _h.end(), hbar, RT());
}

void IdealGasPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    scale(_s.begin(), _s.end(), sbar, GasConstant);
    doublereal logp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (-logp - log(xx));
    }
}

void IdealGasPhase::getPartialMolarIntEnergies(doublereal* ubar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
}

void IdealGasPhase::getPartialMolarCp(doublereal* cpbar) const
{
    const vector_fp& _cp = cp_R_ref();
    scale(_cp.begin(), _cp.end(), cpbar, GasConstant);
}

void IdealGasPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    double vol = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --

void IdealGasPhase::getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

void IdealGasPhase::getEntropy_R(doublereal* sr) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void IdealGasPhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void IdealGasPhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), gpure, RT());
    double tmp = log(pressure() / refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] += tmp;
    }
}

void IdealGasPhase::getIntEnergy_RT(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

void IdealGasPhase::getCp_R(doublereal* cpr) const
{
    const vector_fp& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

void IdealGasPhase::getStandardVolumes(doublereal* vol) const
{
    double tmp = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Thermodynamic Values for the Species Reference States ---------

void IdealGasPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

void IdealGasPhase::getGibbs_RT_ref(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
}

void IdealGasPhase::getGibbs_ref(doublereal* g) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g, RT());
}

void IdealGasPhase::getEntropy_R_ref(doublereal* er) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er);
}

void IdealGasPhase::getIntEnergy_RT_ref(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

void IdealGasPhase::getCp_R_ref(doublereal* cprt) const
{
    const vector_fp& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cprt);
}

void IdealGasPhase::getStandardVolumes_ref(doublereal* vol) const
{
    doublereal tmp = RT() / m_p0;
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

bool IdealGasPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            m_p0 = refPressure();
        }
        m_h0_RT.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_expg0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_s0_R.push_back(0.0);
        m_pp.push_back(0.0);
    }
    return added;
}

void IdealGasPhase::setToEquilState(const doublereal* mu_RT)
{
    const vector_fp& grt = gibbs_RT_ref();

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    doublereal pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        double tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 300.0) {
            double tmp2 = tmp / 300.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(300.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setState_PX(pres, &m_pp[0]);
}

void IdealGasPhase::_updateThermo() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    doublereal tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tnow) {
        m_spthermo.update(tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        cached.state1 = tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
}
}
