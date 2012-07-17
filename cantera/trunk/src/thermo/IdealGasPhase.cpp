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
// Default empty Constructor
IdealGasPhase::IdealGasPhase():
    m_mm(0),
    m_tmin(0.0),
    m_tmax(0.0),
    m_p0(-1.0),
    m_tlast(0.0),
    m_logc0(0.0)
{
}

// Copy Constructor
IdealGasPhase::IdealGasPhase(const IdealGasPhase& right):
    m_mm(right.m_mm),
    m_tmin(right.m_tmin),
    m_tmax(right.m_tmax),
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
IdealGasPhase& IdealGasPhase::
operator=(const IdealGasPhase& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_mm      = right.m_mm;
        m_tmin    = right.m_tmin;
        m_tmax    = right.m_tmax;
        m_p0      = right.m_p0;
        m_tlast   = right.m_tlast;
        m_logc0   = right.m_logc0;
        m_h0_RT   = right.m_h0_RT;
        m_cp0_R   = right.m_cp0_R;
        m_g0_RT   = right.m_g0_RT;
        m_s0_R    = right.m_s0_R;
        m_expg0_RT= right.m_expg0_RT;
        m_pp      = right.m_pp;
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
ThermoPhase* IdealGasPhase::duplMyselfAsThermoPhase() const
{
    ThermoPhase* igp = new IdealGasPhase(*this);
    return (ThermoPhase*) igp;
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
doublereal IdealGasPhase::intEnergy_mole() const
{
    return GasConstant * temperature()
           * (mean_X(&enthalpy_RT_ref()[0]) - 1.0);
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
doublereal IdealGasPhase::entropy_mole() const
{
    return GasConstant * (mean_X(&entropy_R_ref()[0]) -
                          sum_xlogx() - std::log(pressure()/m_spthermo->refPressure()));
}

/*
 * Molar Gibbs free Energy for an ideal gas.
 * Units =  J/kmol.
 */
doublereal IdealGasPhase::gibbs_mole() const
{
    return enthalpy_mole() - temperature() * entropy_mole();
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
doublereal IdealGasPhase::cp_mole() const
{
    return GasConstant * mean_X(&cp_R_ref()[0]);
}

/*
 * Molar heat capacity at constant volume. Units: J/kmol/K.
 * For an ideal gas mixture,
 * \f[ \hat c_v = \hat c_p - \hat R. \f]
 */
doublereal IdealGasPhase::cv_mole() const
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
  doublereal IdealGasPhase::cv_tr (doublereal atomicity) const
  { 
    // k is the species number
    int dum  = 0;
    int type = 0;
    doublereal c[12];
    doublereal minTemp;
    doublereal maxTemp;
    doublereal refPressure;

    m_spthermo->reportParams(dum,type,c,minTemp,maxTemp,refPressure);
    
    if(type != 111)
      {
	throw CanteraError("Error in IdealGasPhase.cpp",
			   "cv_tr only supported for StatMech!. \n\n");
	
      }
    
    // see reportParameters for specific details
    return c[3];
  }
      
  /**
   * @returns species translational specific heat at constant volume.
   */
  doublereal IdealGasPhase::cv_trans () const
  { return 1.5*GasConstant; }

  /**
   * @returns species rotational specific heat at constant volume.
   *
   */
  doublereal IdealGasPhase::cv_rot (double atom) const
  { return std::max(cv_tr(atom) - cv_trans(), 0.); }
     
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
  doublereal IdealGasPhase::cv_vib (const int k, const doublereal T) const
  {

    // k is the species number
    int dum  = 0;
    int type = 0;
    doublereal c[12];
    doublereal minTemp;
    doublereal maxTemp;
    doublereal refPressure;

    c[0] = temperature();

    m_spthermo->reportParams(dum,type,c,minTemp,maxTemp,refPressure);

    // basic sanity check 
    if(type != 111)
      {
	throw CanteraError("Error in IdealGasPhase.cpp",
			   "cv_vib only supported for StatMech!. \n\n");

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
doublereal IdealGasPhase::standardConcentration(size_t k) const
{
    double p = pressure();
    return p/(GasConstant * temperature());
}

/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal IdealGasPhase::logStandardConc(size_t k) const
{
    _updateThermo();
    double p = pressure();
    double lc = std::log(p / (GasConstant * temperature()));
    return lc;
}

/*
 * Get the array of non-dimensional activity coefficients
 */
void IdealGasPhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

/*
 * Get the array of chemical potentials at unit activity \f$
 * \mu^0_k(T,P) \f$.
 */
void IdealGasPhase::getStandardChemPotentials(doublereal* muStar) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), muStar, _RT());
    double tmp = log(pressure() /m_spthermo->refPressure());
    tmp *=  GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        muStar[k] += tmp;  // add RT*ln(P/P_0)
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
        mu[k] += rt*(log(xx));
    }
}

/*
 * Get the array of partial molar enthalpies of the species
 * units = J / kmol
 */
void IdealGasPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    scale(_h.begin(), _h.end(), hbar, rt);
}

/*
 * Get the array of partial molar entropies of the species
 * units = J / kmol / K
 */
void IdealGasPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    doublereal r = GasConstant;
    scale(_s.begin(), _s.end(), sbar, r);
    doublereal logp = log(pressure()/m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += r * (- logp - log(xx));
    }
}

/*
 * Get the array of partial molar internal energies of the species
 * units = J / kmol
 */
void IdealGasPhase::getPartialMolarIntEnergies(doublereal* ubar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    for (size_t k = 0; k < m_kk; k++) {
        ubar[k] =  rt * (_h[k] - 1.0);
    }
}

/*
 * Get the array of partial molar heat capacities
 */
void IdealGasPhase::getPartialMolarCp(doublereal* cpbar) const
{
    const vector_fp& _cp = cp_R_ref();
    scale(_cp.begin(), _cp.end(), cpbar, GasConstant);
}

/*
 * Get the array of partial molar volumes
 * units = m^3 / kmol
 */
void IdealGasPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    double vol = 1.0 / molarDensity();
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] = vol;
    }
}

// Properties of the Standard State of the Species in the Solution --

/*
 * Get the nondimensional Enthalpy functions for the species
 * at their standard states at the current T and P of the
 * solution
 */
void IdealGasPhase::getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

/*
 * Get the array of nondimensional entropy functions for the
 * standard state species
 * at the current <I>T</I> and <I>P</I> of the solution.
 */
void IdealGasPhase::getEntropy_R(doublereal* sr) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
    double tmp = log(pressure() /m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        sr[k] -= tmp;
    }
}

/*
 * Get the nondimensional gibbs function for the species
 * standard states at the current T and P of the solution.
 */
void IdealGasPhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
    double tmp = log(pressure() /m_spthermo->refPressure());
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] += tmp;
    }
}

/*
 * get the pure Gibbs free energies of each species assuming
 * it is in its standard state. This is the same as
 * getStandardChemPotentials().
 */
void IdealGasPhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), gpure, _RT());
    double tmp = log(pressure() /m_spthermo->refPressure());
    tmp *= _RT();
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] += tmp;
    }
}

/*
 *  Returns the vector of nondimensional
 *  internal Energies of the standard state at the current temperature
 *  and pressure of the solution for each species.
 */
void IdealGasPhase::getIntEnergy_RT(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

/*
 * Get the nondimensional heat capacity at constant pressure
 * function for the species
 * standard states at the current T and P of the solution.
 */
void IdealGasPhase::getCp_R(doublereal* cpr) const
{
    const vector_fp& _cpr = cp_R_ref();
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
void IdealGasPhase::getStandardVolumes(doublereal* vol) const
{
    double tmp = 1.0 / molarDensity();
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
void IdealGasPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    copy(_h.begin(), _h.end(), hrt);
}

/*
 *  Returns the vector of nondimensional
 *  enthalpies of the reference state at the current temperature
 *  and reference pressure.
 */
void IdealGasPhase::getGibbs_RT_ref(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    copy(gibbsrt.begin(), gibbsrt.end(), grt);
}

/*
 *  Returns the vector of the
 *  gibbs function of the reference state at the current temperature
 *  and reference pressure.
 *  units = J/kmol
 */
void IdealGasPhase::getGibbs_ref(doublereal* g) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), g, _RT());
}

/*
 *  Returns the vector of nondimensional
 *  entropies of the reference state at the current temperature
 *  and reference pressure.
 */
void IdealGasPhase::getEntropy_R_ref(doublereal* er) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), er);
}

/*
 *  Returns the vector of nondimensional
 *  internal Energies of the reference state at the current temperature
 *  of the solution and the reference pressure for each species.
 */
void IdealGasPhase::getIntEnergy_RT_ref(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - 1.0;
    }
}

/*
 *  Returns the vector of nondimensional
 *  constant pressure heat capacities of the reference state
 *   at the current temperature and reference pressure.
 */
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


// new methods defined here -------------------------------


void IdealGasPhase::initThermo()
{

    m_mm = nElements();
    doublereal tmin = m_spthermo->minTemp();
    doublereal tmax = m_spthermo->maxTemp();
    if (tmin > 0.0) {
        m_tmin = tmin;
    }
    if (tmax > 0.0) {
        m_tmax = tmax;
    }
    m_p0 = refPressure();

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
    setState_PX(pres, &m_pp[0]);
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
void IdealGasPhase::_updateThermo() const
{
    doublereal tnow = temperature();

    // If the temperature has changed since the last time these
    // properties were computed, recompute them.
    if (m_tlast != tnow) {
        m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0],
                           &m_s0_R[0]);
        m_tlast = tnow;

        // update the species Gibbs functions
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_logc0 = log(m_p0/(GasConstant * tnow));
        m_tlast = tnow;
    }
}
}

