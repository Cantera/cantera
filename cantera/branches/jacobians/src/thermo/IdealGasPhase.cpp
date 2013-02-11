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
// Default empty Constructor
template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase() :
        m_p0(-1.0),
        m_tlast(0.0),
        m_logc0(0.0)
{
}

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase(const std::string& inputFile, const std::string& id) :
    m_p0(-1.0),
    m_tlast(0.0),
    m_logc0(0.0)
{
    this->initThermoFile(inputFile, id);
}

template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase(XML_Node& phaseRef, const std::string& id) :
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

// Copy Constructor
template<typename ValAndDerivType>
IdealGasPhase<ValAndDerivType>::IdealGasPhase(const IdealGasPhase& right) :
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

// Assignment operator
/*
 * Assignment operator for the object. Constructed
 * object will be a clone of this object, but will
 * also own all of its data.
 *
 * @param right Object to be copied.
 */
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

// Duplicator from the %ThermoPhase parent class
/*
 * Given a pointer to a %ThermoPhase object, this function will
 * duplicate the %ThermoPhase object and all underlying structures.
 * This is basically a wrapper around the copy constructor.
 *
 * @return returns a pointer to a %ThermoPhase
 */
template<typename ValAndDerivType>
ThermoPhase<ValAndDerivType>* IdealGasPhase<ValAndDerivType>::duplMyselfAsThermoPhase() const
{
    return new IdealGasPhase(*this);
}

// Molar Thermodynamic Properties of the Solution ------------------

/*
 * Molar internal energy. J/kmol. For an ideal gas mixture,
 * \f[
 * \hat u(T) = \sum_k X_k \hat h^0_k(T) - \hat R T,
 * \f]
 * and is a function only of temperature.
 * The reference-state pure-species enthalpies
 * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic
 * property manager.
 * @see SpeciesThermo
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::intEnergy_mole() const
{
    return GasConstant * this->temperature() * (mean_X(&enthalpy_RT_ref()[0]) - 1.0);
}

/*
 * Molar entropy. Units: J/kmol/K.
 * For an ideal gas mixture,
 * \f[
 * \hat s(T, P) = \sum_k X_k \hat s^0_k(T) - \hat R \log (P/P^0).
 * \f]
 * The reference-state pure-species entropies
 * \f$ \hat s^0_k(T) \f$ are computed by the species thermodynamic
 * property manager.
 * @see SpeciesThermo
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::entropy_mole() const
{
    return GasConstant
            * (mean_X(&entropy_R_ref()[0]) - this->sum_xlogx() - std::log(pressure() / (this->m_spthermo)->refPressure()));
}

/*
 * Molar Gibbs free Energy for an ideal gas.
 * Units =  J/kmol.
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::gibbs_mole() const
{
    return enthalpy_mole() - this->temperature() * entropy_mole();
}

/*
 * Molar heat capacity at constant pressure. Units: J/kmol/K.
 * For an ideal gas mixture,
 * \f[
 * \hat c_p(t) = \sum_k \hat c^0_{p,k}(T).
 * \f]
 * The reference-state pure-species heat capacities
 * \f$ \hat c^0_{p,k}(T) \f$ are computed by the species thermodynamic
 * property manager.
 * @see SpeciesThermo
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cp_mole() const
{
    return GasConstant * mean_X(&cp_R_ref()[0]);
}

/*
 * Molar heat capacity at constant volume. Units: J/kmol/K.
 * For an ideal gas mixture,
 * \f[ \hat c_v = \hat c_p - \hat R. \f]
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_mole() const
{
    return cp_mole() - GasConstant;
}

/**
 * @returns species translational/rotational specific heat at
 * constant volume.
 *
 * Either: $5/2 R_s$ or $3/2 R_s$ for molecules/atoms.
 *
 */
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

/**
 * @returns species translational specific heat at constant volume.
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_trans() const
{
    return 1.5 * GasConstant;
}

/**
 * @returns species rotational specific heat at constant volume.
 *
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::cv_rot(double atom) const
{
    return std::max(cv_tr(atom) - cv_trans(), 0.);
}

/**
 * @returns species vibrational specific heat at
 * constant volume.
 *
 * C^{vib}_{v,s} = \frac{\partial e^{vib}_{v,s} }{\partial T}
 *
 * The species vibration energy ($e^{vib}_{v,s}$) is:
 *
 * 0: atom
 *
 * Diatomic:
 * \f[
 * \frac{R_s \theta_{v,s}}{e^{\theta_{v,s}/T}-1}
 * \f]
 *
 * General Molecules:
 * \f[
 * \sum_i \frac{R_s \theta_{v,s,i}}{e^{\theta_{v,s,i}/T}-1}
 * \f]
 *
 */
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

// Mechanical Equation of State ----------------------------
// Chemical Potentials and Activities ----------------------

/*
 * Returns the standard concentration \f$ C^0_k \f$, which is used to normalize
 * the generalized concentration.
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::standardConcentration(size_t k) const
{
    ValAndDerivType p = pressure();
    return p / (GasConstant * this->temperature());
}

/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
template<typename ValAndDerivType>
ValAndDerivType IdealGasPhase<ValAndDerivType>::logStandardConc(size_t k) const
{
    _updateThermo();
    ValAndDerivType p = pressure();
    ValAndDerivType lc = std::log(p / (GasConstant * this->temperature()));
    return lc;
}

/*
 * Get the array of non-dimensional activity coefficients
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getActivityCoefficients(ValAndDerivType* ac) const
{
    for (size_t k = 0; k < this->m_kk; k++) {
        ac[k] = 1.0;
    }
}

/*
 * Get the array of chemical potentials at unit activity \f$
 * \mu^0_k(T,P) \f$.
 */
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

/*
 * Get the array of partial molar enthalpies of the species
 * units = J / kmol
 */
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

/*
 * Get the array of partial molar entropies of the species
 * units = J / kmol / K
 */
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

/*
 * Get the array of partial molar internal energies of the species
 * units = J / kmol
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarIntEnergies(ValAndDerivType* ubar) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    ValAndDerivType rt = GasConstant * this->temperature();
    for (size_t k = 0; k < this->m_kk; k++) {
        ubar[k] = rt * (_h[k] - 1.0);
    }
}

/*
 * Get the array of partial molar heat capacities
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarCp(ValAndDerivType* cpbar) const
{
    const vector_ValAndDeriv& _cp = cp_R_ref();
    // scale(_cp.begin(), _cp.end(), cpbar, GasConstant);
    for (size_t k = 0; k < this->m_kk; k++) {
        cpbar[k] = GasConstant * _cp[k];
    }
}

/*
 * Get the array of partial molar volumes
 * units = m^3 / kmol
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getPartialMolarVolumes(ValAndDerivType* vbar) const
{
    ValAndDerivType vol = 1.0 / this->molarDensity();
    for (size_t k = 0; k < this->m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --

/*
 * Get the nondimensional Enthalpy functions for the species
 * at their standard states at the current T and P of the
 * solution
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEnthalpy_RT(ValAndDerivType* hrt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

/*
 * Get the array of nondimensional entropy functions for the
 * standard state species
 * at the current <I>T</I> and <I>P</I> of the solution.
 */
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

/*
 * Get the nondimensional gibbs function for the species
 * standard states at the current T and P of the solution.
 */
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

/*
 * get the pure Gibbs free energies of each species assuming
 * it is in its standard state. This is the same as
 * getStandardChemPotentials().
 */
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

/*
 *  Returns the vector of nondimensional
 *  internal Energies of the standard state at the current temperature
 *  and pressure of the solution for each species.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getIntEnergy_RT(ValAndDerivType* urt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

/*
 * Get the nondimensional heat capacity at constant pressure
 * function for the species
 * standard states at the current T and P of the solution.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getCp_R(ValAndDerivType* cpr) const
{
    const vector_ValAndDeriv& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

/*
 *  Get the molar volumes of the species standard states at the current
 *  <I>T</I> and <I>P</I> of the solution.
 *  units = m^3 / kmol
 *
 * @param vol     Output vector containing the standard state volumes.
 *                Length: m_kk.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getStandardVolumes(ValAndDerivType* vol) const
{
    ValAndDerivType tmp = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vol[k] = tmp;
    }
}

// Thermodynamic Values for the Species Reference States ---------

/*
 *  Returns the vector of nondimensional
 *  enthalpies of the reference state at the current temperature
 *  and reference pressure.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEnthalpy_RT_ref(ValAndDerivType* hrt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

/*
 *  Returns the vector of nondimensional
 *  enthalpies of the reference state at the current temperature
 *  and reference pressure.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getGibbs_RT_ref(ValAndDerivType* grt) const
{
    const vector_ValAndDeriv& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
}

/*
 *  Returns the vector of the
 *  gibbs function of the reference state at the current temperature
 *  and reference pressure.
 *  units = J/kmol
 */
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

/*
 *  Returns the vector of nondimensional
 *  entropies of the reference state at the current temperature
 *  and reference pressure.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getEntropy_R_ref(ValAndDerivType* er) const
{
    const vector_ValAndDeriv& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er);
}

/*
 *  Returns the vector of nondimensional
 *  internal Energies of the reference state at the current temperature
 *  of the solution and the reference pressure for each species.
 */
template<typename ValAndDerivType>
void IdealGasPhase<ValAndDerivType>::getIntEnergy_RT_ref(ValAndDerivType* urt) const
{
    const vector_ValAndDeriv& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

/*
 *  Returns the vector of nondimensional
 *  constant pressure heat capacities of the reference state
 *   at the current temperature and reference pressure.
 */
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

// new methods defined here -------------------------------

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

/*
 * Set mixture to an equilibrium state consistent with specified
 * chemical potentials and temperature. This method is needed by
 * the ChemEquil equilibrium solver.
 */
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
        tmp = -(grt[k]).val() + mu_RT[k];
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

