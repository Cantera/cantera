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

namespace Cantera {

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase() :
        ThermoPhase<ValAndDerivType>(),
        m_p0(-1.0),
        m_tlast(0.0),
        m_logc0(0.0)
{
}

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase(const std::string& inputFile, const std::string& id) :
        ThermoPhase<ValAndDerivType>(),
        m_p0(-1.0),
        m_tlast(0.0),
        m_logc0(0.0)
{
    this->initThermoFile(inputFile, id);
}

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase(XML_Node& phaseRef, const std::string& id) :
        ThermoPhase<ValAndDerivType>(),
        m_p0(-1.0),
        m_tlast(0.0),
        m_logc0(0.0)
{
    this->initThermoXML(phaseRef, id);
}

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::~IdealGasPhase()
{
}

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase(const IdealGasPhase& right) :
        ThermoPhase<ValAndDerivType>(),
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

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>& IdealGasPhase<ValAndDerivType>::operator=(const IdealGasPhase<ValAndDerivType>& right)
{
    if (&right != this) {
        ThermoPhase<ValAndDerivType>::operator=(right);
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

template<typename ValAndDerivType>
ThermoPhase<ValAndDerivType>* IdealGasPhase<ValAndDerivType>::duplMyselfAsThermoPhase() const
{
    return new IdealGasPhase(*this);
}

// Molar Thermodynamic Properties of the Solution ------------------

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::intEnergy_mole() const
{
    return GasConstant * this->temperature() * (mean_X(&enthalpy_RT_ref()[0]) - 1.0);
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::entropy_mole() const
{
    return GasConstant
            * (mean_X(&entropy_R_ref()[0]) - this->sum_xlogx() - std::log(pressure() / (this->m_spthermo)->refPressure()));
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::gibbs_mole() const
{
    return enthalpy_mole() - this->temperature() * entropy_mole();
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cp_mole() const
{
    return GasConstant * mean_X(&cp_R_ref()[0]);
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_mole() const
{
    return cp_mole() - GasConstant;
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_tr(doublereal atomicity) const
{
    // k is the species number
    int dum = 0;
    int type = 0;
    doublereal c[12];
    doublereal minTemp;
    doublereal maxTemp;
    doublereal refPressure;

    (this->m_spthermo)->reportParams(dum, type, c, minTemp, maxTemp, refPressure);

    if (type != 111) {
        throw CanteraError("Error in IdealGasPhase.cpp", "cv_tr only supported for StatMech!. \n\n");

    }

    // see reportParameters for specific details
    return c[3];
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_trans() const
{
    return 1.5 * GasConstant;
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_rot(double atom) const
{
    return std::max(cv_tr(atom) - cv_trans(), 0.);
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_vib(const int k, const doublereal T) const
{

    // k is the species number
    int dum = 0;
    int type = 0;
    doublereal c[12];
    doublereal minTemp;
    doublereal maxTemp;
    doublereal refPressure;

    c[0] = this->temperature();

    (this->m_spthermo)->reportParams(dum, type, c, minTemp, maxTemp, refPressure);

    // basic sanity check
    if (type != 111) {
        throw CanteraError("Error in IdealGasPhase.cpp", "cv_vib only supported for StatMech!. \n\n");

    }

    // see reportParameters for specific details
    return c[4];

}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::standardConcentration(size_t k) const
{
    ValAndDerivType p = pressure();
    return p / (GasConstant * this->temperature());
}

template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::logStandardConc(size_t k) const
{
    _updateThermo();
    ValAndDerivType p = pressure();
    ValAndDerivType lc = std::log(p / (GasConstant * this->temperature()));
    return lc;
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getActivityCoefficients(ValAndDerivType* ac) const
{
    for (size_t k = 0; k < this->m_kk; k++) {
        ac[k] = 1.0;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getStandardChemPotentials(ValAndDerivType* muStar) const
{
    const vector_ValAndDeriv& gibbsrt = gibbs_RT_ref();
    //scale(gibbsrt.begin(), gibbsrt.end(), muStar, this->_RT());
    ValAndDerivType rt = this->_RT();
    for (size_t k = 0; k < this->m_kk; k++) {
        muStar[k] = rt * gibbsrt[k];
    }
    ValAndDerivType tmp = log(pressure() / (this->m_spthermo)->refPressure());
    tmp *= GasConstant * this->temperature();
    for (size_t k = 0; k < this->m_kk; k++) {
        muStar[k] += tmp; // add RT*ln(P/P_0)
    }
}

//  Partial Molar Properties of the Solution --------------

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getChemPotentials(ValAndDerivType* mu) const
{
    getStandardChemPotentials(mu);
    //doublereal logp = log(pressure()/m_spthermo->refPressure());
    ValAndDerivType xx;
    ValAndDerivType rt = this->temperature() * GasConstant;
    //const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, this->moleFraction(k));
        mu[k] += rt * (log(xx));
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarEnthalpies(ValAndDerivType* hbar) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    ValAndDerivType rt = GasConstant * this->temperature();
    //scale(_h.begin(), _h.end(), hbar, rt);
    for (size_t k = 0; k < this->m_kk; k++) {
        hbar[k] = rt * _h[k];
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarEntropies(ValAndDerivType* sbar) const
{
    const vector_ValAndDeriv& _s = entropy_R_ref();
    doublereal r = GasConstant;
    // scale(_s.begin(), _s.end(), sbar, r);
    for (size_t k = 0; k < this->m_kk; k++) {
        sbar[k] = r * _s[k];
    }
    ValAndDerivType logp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < this->m_kk; k++) {
        ValAndDerivType xx = std::max(SmallNumber, this->moleFraction(k));
        sbar[k] += r * (-logp - log(xx));
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarIntEnergies(ValAndDerivType* ubar) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    ValAndDerivType rt = GasConstant * this->temperature();
    for (size_t k = 0; k < this->m_kk; k++) {
        ubar[k] = rt * (_h[k] - 1.0);
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarCp(ValAndDerivType* cpbar) const
{
    const vector_ValAndDeriv& _cp = cp_R_ref();
    // scale(_cp.begin(), _cp.end(), cpbar, GasConstant);
    for (size_t k = 0; k < this->m_kk; k++) {
        cpbar[k] = GasConstant * _cp[k];
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarVolumes(ValAndDerivType* vbar) const
{
    ValAndDerivType vol = 1.0 / this->molarDensity();
    for (size_t k = 0; k < this->m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEnthalpy_RT(ValAndDerivType* hrt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEntropy_R(ValAndDerivType* sr) const
{
    const vector_ValAndDeriv& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    ValAndDerivType tmp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < this->m_kk; k++) {
        sr[k] -= tmp;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getGibbs_RT(ValAndDerivType* grt) const
{
    const vector_ValAndDeriv& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    ValAndDerivType tmp = log(pressure() / m_spthermo->refPressure());
    for (size_t k = 0; k < this->m_kk; k++) {
        grt[k] += tmp;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPureGibbs(ValAndDerivType* gpure) const
{
    const vector_ValAndDeriv& gibbsrt = gibbs_RT_ref();
    // scale(gibbsrt.begin(), gibbsrt.end(), gpure, this->_RT());
    ValAndDerivType rt = this->_RT();
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = rt * gibbsrt[k];
    }
    ValAndDerivType tmp = log(pressure() / m_spthermo->refPressure());
    tmp *= _RT();
    for (size_t k = 0; k < this->m_kk; k++) {
        gpure[k] += tmp;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getIntEnergy_RT(ValAndDerivType* urt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getCp_R(ValAndDerivType* cpr) const
{
    const vector_ValAndDeriv& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getStandardVolumes(ValAndDerivType* vol) const
{
    ValAndDerivType tmp = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Thermodynamic Values for the Species Reference States ---------

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEnthalpy_RT_ref(ValAndDerivType* hrt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getGibbs_RT_ref(ValAndDerivType* grt) const
{
    const vector_ValAndDeriv& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getGibbs_ref(ValAndDerivType* g) const
{
    const vector_ValAndDeriv& gibbsrt = gibbs_RT_ref();
    // scale(gibbsrt.begin(), gibbsrt.end(), g, _RT());
    ValAndDerivType rt = this->_RT();
    for (size_t k = 0; k < m_kk; k++) {
        g[k] = rt * gibbsrt[k];
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEntropy_R_ref(ValAndDerivType* er) const
{
    const vector_ValAndDeriv& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getIntEnergy_RT_ref(ValAndDerivType* urt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getCp_R_ref(ValAndDerivType* cprt) const
{
    const vector_ValAndDeriv& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cprt);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getStandardVolumes_ref(ValAndDerivType* vol) const
{
    ValAndDerivType tmp = _RT() / m_p0;
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::initThermo()
{
    m_p0 = this->refPressure();
    m_h0_RT.resize(m_kk);
    m_g0_RT.resize(m_kk);
    m_expg0_RT.resize(m_kk);
    m_cp0_R.resize(m_kk);
    m_s0_R.resize(m_kk);
    m_pp.resize(m_kk);
}

template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::setToEquilState(const doublereal* mu_RT)
{
    ValAndDerivType tmp;
    ValAndDerivType tmp2;
    const vector_ValAndDeriv& grt = gibbs_RT_ref();

    /*
     * Within the method, we protect against inf results if the
     * exponent is too high.
     *
     * If it is too low, we set
     * the partial pressure to zero. This capability is needed
     * by the elemental potential method.
     */
    ValAndDerivType pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    this->setState_PX(pres, &m_pp[0]);
}

template<>
void IdealGasPhase<doubleFAD>::setToEquilState(const doublereal* mu_RT)
{
    doublereal tmp;
    doublereal tmp2;
    const vector_ValAndDeriv& grt = gibbs_RT_ref();

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
        tmp = - (grt[k]).val() + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    this->setState_PX(pres, &m_pp[0]);
}

/// This method is called each time a thermodynamic property is
/// requested, to check whether the internal species properties
/// within the object need to be updated.
/// Currently, this updates the species thermo polynomial values
/// for the current value of the temperature. A check is made
/// to see if the temperature has changed since the last
/// evaluation. This object does not contain any persistent
/// data that depends on the concentration, that needs to be
/// updated. The state object modifies its concentration
/// dependent information at the time the setMoleFractions()
/// (or equivalent) call is made.
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::_updateThermo() const
{
    doublereal tnow = this->temperature();

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

// Explicit Instantiation Section
template class IdealGasPhase<doublereal> ;
#ifdef INDEPENDENT_VARIABLE_DERIVATIVES
#ifdef HAS_SACADO
template class IdealGasPhase<doubleFAD> ;
#endif
#endif

}
