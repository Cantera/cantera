/**
 *  @file IdealGasPhase.cpp
 *   ThermoPhase object for the ideal gas equation of
 * state - workhorse for %Cantera (see @ref thermoprops
 * and class @link Cantera::IdealGasPhase IdealGasPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

IdealGasPhase::IdealGasPhase(const string& inputFile, const string& id_)
{
    initThermoFile(inputFile, id_);
}

// Molar Thermodynamic Properties of the Solution ------------------

double IdealGasPhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx() - std::log(pressure() / refPressure()));
}

double IdealGasPhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

double IdealGasPhase::cv_mole() const
{
    return cp_mole() - GasConstant;
}

double IdealGasPhase::soundSpeed() const {
    return sqrt(
        cp_mole() / cv_mole() * GasConstant / meanMolecularWeight() * temperature()
    );
}

double IdealGasPhase::standardConcentration(size_t k) const
{
    return pressure() / RT();
}

void IdealGasPhase::getActivityCoefficients(span<double> ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

void IdealGasPhase::getStandardChemPotentials(span<double> muStar) const
{
    getGibbs_ref(muStar);
    double tmp = log(pressure() / refPressure()) * RT();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp; // add RT*ln(P/P_0)
    }
}

//  Partial Molar Properties of the Solution --------------

void IdealGasPhase::getChemPotentials(span<double> mu) const
{
    getStandardChemPotentials(mu);
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT() * log(xx);
    }
}

void IdealGasPhase::getPartialMolarEnthalpies(span<double> hbar) const
{
    auto _h = enthalpy_RT_ref();
    scale(_h.begin(), _h.end(), hbar.begin(), RT());
}

void IdealGasPhase::getPartialMolarEntropies(span<double> sbar) const
{
    auto _s = entropy_R_ref();
    scale(_s.begin(), _s.end(), sbar.begin(), GasConstant);
    double logp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (-logp - log(xx));
    }
}

void IdealGasPhase::getPartialMolarIntEnergies(span<double> ubar) const
{
    auto _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = RT() * (_h[k] - 1.0);
    }
}

void IdealGasPhase::getPartialMolarCp(span<double> cpbar) const
{
    auto _cp = cp_R_ref();
    scale(_cp.begin(), _cp.end(), cpbar.begin(), GasConstant);
}

void IdealGasPhase::getPartialMolarVolumes(span<double> vbar) const
{
    double vol = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --

void IdealGasPhase::getEnthalpy_RT(span<double> hrt) const
{
    auto _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt.begin());
}

void IdealGasPhase::getEntropy_R(span<double> sr) const
{
    auto _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr.begin());
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void IdealGasPhase::getGibbs_RT(span<double> grt) const
{
    auto gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt.begin());
    double tmp = log(pressure() / refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void IdealGasPhase::getIntEnergy_RT(span<double> urt) const
{
    getIntEnergy_RT_ref(urt);
}

void IdealGasPhase::getCp_R(span<double> cpr) const
{
    auto _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr.begin());
}

void IdealGasPhase::getStandardVolumes(span<double> vol) const
{
    double tmp = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Thermodynamic Values for the Species Reference States ---------

void IdealGasPhase::getEnthalpy_RT_ref(span<double> hrt) const
{
    auto _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt.begin());
}

void IdealGasPhase::getGibbs_RT_ref(span<double> grt) const
{
    auto gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt.begin());
}

void IdealGasPhase::getGibbs_ref(span<double> g) const
{
    auto gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g.begin(), RT());
}

void IdealGasPhase::getEntropy_R_ref(span<double> er) const
{
    auto _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er.begin());
}

void IdealGasPhase::getIntEnergy_RT_ref(span<double> urt) const
{
    auto _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

void IdealGasPhase::getCp_R_ref(span<double> cprt) const
{
    auto _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cprt.begin());
}

void IdealGasPhase::getStandardVolumes_ref(span<double> vol) const
{
    double tmp = RT() / m_p0;
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

void IdealGasPhase::setToEquilState(span<const double> mu_RT)
{
    auto grt = gibbs_RT_ref();

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    double pres = 0.0;
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
    setMoleFractions(m_pp);
    setPressure(pres);
}

void IdealGasPhase::updateThermo() const
{
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tnow) {
        m_spthermo.update(tnow, m_cp0_R, m_h0_RT, m_s0_R);
        cached.state1 = tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
}
}
