/**
 *  @file IdealGasPhase.cpp
 *   ThermoPhase object for the ideal gas equation of
 * state - workhorse for %Cantera (see \ref thermoprops
 * and class \link Cantera::IdealGasPhase IdealGasPhase\endlink).
 */

#include "cantera/base/ct_defs.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/SpeciesThermo.h"

using namespace std;

namespace Cantera
{

IdealGasPhase::IdealGasPhase() :
    m_p0(-1.0),
    m_tlast(0.0),
    m_logc0(0.0)
{
}

IdealGasPhase::IdealGasPhase(const std::string& inputFile, const std::string& id_) :
    m_p0(-1.0),
    m_tlast(0.0),
    m_logc0(0.0)
{
    initThermoFile(inputFile, id_);
}

IdealGasPhase::IdealGasPhase(XML_Node& phaseRef, const std::string& id_) :
    m_p0(-1.0),
    m_tlast(0.0),
    m_logc0(0.0)
{
    initThermoXML(phaseRef, id_);
}

IdealGasPhase::IdealGasPhase(const IdealGasPhase& right) :
    m_p0(right.m_p0),
    m_tlast(right.m_tlast),
    m_logc0(right.m_logc0)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = right;
}

IdealGasPhase& IdealGasPhase::operator=(const IdealGasPhase& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_p0 = right.m_p0;
        m_tlast = right.m_tlast;
        m_logc0 = right.m_logc0;
        m_h0_RT = right.m_h0_RT;
        m_cp0_R = right.m_cp0_R;
        m_g0_RT = right.m_g0_RT;
        m_s0_R = right.m_s0_R;
        m_expg0_RT = right.m_expg0_RT;
        m_pp = right.m_pp;
    }
    return *this;
}

ThermoPhase* IdealGasPhase::duplMyselfAsThermoPhase() const
{
    return new IdealGasPhase(*this);
}

// Molar Thermodynamic Properties of the Solution ------------------

doublereal IdealGasPhase::intEnergy_mole() const
{
    return GasConstant * temperature() * (mean_X(&enthalpy_RT_ref()[0]) - 1.0);
}

doublereal IdealGasPhase::entropy_mole() const
{
    return GasConstant * (mean_X(&entropy_R_ref()[0]) - sum_xlogx() - std::log(pressure() / m_spthermo->refPressure()));
}

doublereal IdealGasPhase::gibbs_mole() const
{
    return enthalpy_mole() - temperature() * entropy_mole();
}

doublereal IdealGasPhase::cp_mole() const
{
    return GasConstant * mean_X(&cp_R_ref()[0]);
}

doublereal IdealGasPhase::cv_mole() const
{
    return cp_mole() - GasConstant;
}

doublereal IdealGasPhase::cv_tr(doublereal atomicity) const
{
    // k is the species number
    int dum = 0;
    int type = m_spthermo->reportType();
    doublereal c[12];
    doublereal minTemp_;
    doublereal maxTemp_;
    doublereal refPressure_;

    if (type != 111) {
        throw CanteraError("Error in IdealGasPhase.cpp", "cv_tr only supported for StatMech!. \n\n");
    }

    m_spthermo->reportParams(dum, type, c, minTemp_, maxTemp_, refPressure_);

    // see reportParameters for specific details
    return c[3];
}

doublereal IdealGasPhase::cv_trans() const
{
    return 1.5 * GasConstant;
}

doublereal IdealGasPhase::cv_rot(double atom) const
{
    return std::max(cv_tr(atom) - cv_trans(), 0.);
}

doublereal IdealGasPhase::cv_vib(const int k, const doublereal T) const
{

    // k is the species number
    int dum = 0;
    int type = m_spthermo->reportType();
    doublereal c[12];
    doublereal minTemp_;
    doublereal maxTemp_;
    doublereal refPressure_;

    c[0] = temperature();

    // basic sanity check
    if (type != 111) {
        throw CanteraError("Error in IdealGasPhase.cpp", "cv_vib only supported for StatMech!. \n\n");
    }

    m_spthermo->reportParams(dum, type, c, minTemp_, maxTemp_, refPressure_);

    // see reportParameters for specific details
    return c[4];

}

doublereal IdealGasPhase::standardConcentration(size_t k) const
{
    double p = pressure();
    return p / (GasConstant * temperature());
}

doublereal IdealGasPhase::logStandardConc(size_t k) const
{
    _updateThermo();
    double p = pressure();
    return std::log(p / (GasConstant * temperature()));
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
    scale(gibbsrt.begin(), gibbsrt.end(), muStar, _RT());
    double tmp = log(pressure() / m_spthermo->refPressure());
    tmp *= GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp; // add RT*ln(P/P_0)
    }
}

//  Partial Molar Properties of the Solution --------------

void IdealGasPhase::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
    //doublereal logp = log(pressure()/m_spthermo->refPressure());
    doublereal xx;
    doublereal rt = temperature() * GasConstant;
    //const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += rt * (log(xx));
    }
}

void IdealGasPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    scale(_h.begin(), _h.end(), hbar, rt);
}

void IdealGasPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    doublereal r = GasConstant;
    scale(_s.begin(), _s.end(), sbar, r);
    doublereal logp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += r * (-logp - log(xx));
    }
}

void IdealGasPhase::getPartialMolarIntEnergies(doublereal* ubar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] = rt * (_h[k] - 1.0);
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
    double tmp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

void IdealGasPhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

void IdealGasPhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), gpure, _RT());
    double tmp = log(pressure() / m_spthermo->refPressure());
    tmp *= _RT();
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
    scale(gibbsrt.begin(), gibbsrt.end(), g, _RT());
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
    doublereal tmp = _RT() / m_p0;
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

void IdealGasPhase::initThermo()
{
    m_p0 = refPressure();
    m_h0_RT.resize(m_kk);
    m_g0_RT.resize(m_kk);
    m_expg0_RT.resize(m_kk);
    m_cp0_R.resize(m_kk);
    m_s0_R.resize(m_kk);
    m_pp.resize(m_kk);
}

void IdealGasPhase::setToEquilState(const doublereal* mu_RT)
{
    double tmp, tmp2;
    const vector_fp& grt = gibbs_RT_ref();

    /*
     * Within the method, we protect against inf results if the
     * exponent is too high.
     *
     * If it is too low, we set
     * the partial pressure to zero. This capability is needed
     * by the elemental potential method.
     */
    doublereal pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 300.0) {
            tmp2 = tmp / 300.;
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
    doublereal tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (m_tlast != tnow) {
        m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        m_tlast = tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_logc0 = log(m_p0 / (GasConstant * tnow));
        m_tlast = tnow;
    }
}
}
