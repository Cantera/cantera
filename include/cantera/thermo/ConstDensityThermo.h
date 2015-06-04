/**
 *  @file ConstDensityThermo.h
 * Header for a Thermo manager for incompressible ThermoPhases
 * (see \ref thermoprops and \link Cantera::ConstDensityThermo ConstDensityThermo\endlink).
 */
/*
 *  Copyright 2002 California Institute of Technology
 */

#ifndef CT_CONSTRHOTHERMO_H
#define CT_CONSTRHOTHERMO_H

#include "ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

//!   Overloads the virtual methods of class ThermoPhase to implement the
//!   incompressible equation of state.
/**
 * <b> Specification of Solution Thermodynamic Properties </b>
 *
 * The density is assumed to be constant, no matter what the concentration of the solution.
 *
 * @ingroup thermoprops
 */
class ConstDensityThermo : public ThermoPhase
{
public:
    //! Constructor.
    ConstDensityThermo() {}

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    ConstDensityThermo(const ConstDensityThermo& right);

    //! Assignment Operator
    /*!
     * @param right Object to be copied
     */
    ConstDensityThermo& operator=(const ConstDensityThermo& right);

    //! Duplication routine for objects which inherit from ThermoPhase
    /*!
     *  This virtual routine can be used to duplicate objects
     *  derived from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Returns a constant corresponding to this class's equation of state
    virtual int eosType() const;

    /// Molar enthalpy. Units: J/kmol.
    virtual doublereal enthalpy_mole() const;

    /// Molar entropy. Units: J/kmol/K.
    virtual doublereal entropy_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    virtual doublereal cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    virtual doublereal pressure() const;

    //! Set the internally stored pressure (Pa) at constant
    //! temperature and composition
    /*!
     *  @param p input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in the
     * thermodynamic functions.  These activity (or generalized)
     * concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. Note that they may
     * or may not have units of concentration --- they might be
     * partial pressures, mole fractions, or surface coverages,
     * for example.
     *
     * @param c Output array of generalized concentrations. The
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;

    //! Get the array of chemical potentials at unit activity for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu0     Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration. In many cases, this quantity
     * will be the same for all species in a phase - for example,
     * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
     * reason, this method returns a single value, instead of an
     * array.  However, for phases in which the standard
     * concentration is species-specific (e.g. surface species of
     * different sizes), this method may be called with an
     * optional parameter indicating the species.
     *
     * @param k Optional parameter indicating the species. The default
     *         is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of  standard state Gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const {
        const vector_fp& gibbsrt = gibbs_RT();
        scale(gibbsrt.begin(), gibbsrt.end(), gpure, _RT());
    }

    //! Get the nondimensional Enthalpy functions for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    void getEnthalpy_RT(doublereal* hrt) const {
        const vector_fp& _h = enthalpy_RT();
        std::copy(_h.begin(), _h.end(), hrt);
    }

    //! Get the array of nondimensional Entropy functions for the
    //! standard state species at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    void getEntropy_R(doublereal* sr) const {
        const vector_fp& _s = entropy_R();
        std::copy(_s.begin(), _s.end(), sr);
    }

    //! Get the nondimensional Gibbs functions for the species
    //! in their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state Gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const {
        const vector_fp& gibbsrt = gibbs_RT();
        std::copy(gibbsrt.begin(), gibbsrt.end(), grt);
    }

    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the species standard states
    //! at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat capacities
     *              Length: m_kk.
     */
    void getCp_R(doublereal* cpr) const {
        const vector_fp& _cpr = cp_R();
        std::copy(_cpr.begin(), _cpr.end(), cpr);
    }

    //!  Returns a reference to the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    const vector_fp& enthalpy_RT() const {
        _updateThermo();
        return m_h0_RT;
    }

    //!  Returns a reference to the vector of nondimensional
    //!  Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    const vector_fp& gibbs_RT() const {
        _updateThermo();
        return m_g0_RT;
    }

    //!  Returns a reference to the vector of nondimensional
    //!  entropies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    const vector_fp& entropy_R() const {
        _updateThermo();
        return m_s0_R;
    }

    //!  Returns a reference to the vector of nondimensional
    //!  constant pressure heat capacities of the reference state
    //!  at the current temperature of the solution
    //!  and reference pressure for each species.
    const vector_fp& cp_R() const {
        _updateThermo();
        return m_cp0_R;
    }

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from function importPhase().
     */
    virtual void initThermo();

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT Input vector of dimensionless element potentials
     *                  The length is equal to nElements().
     */
    virtual void setToEquilState(const doublereal* lambda_RT);

    //! Set the equation of state parameters
    /*!
     * @internal
     *  The number and meaning of these depends on the subclass.
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void setParameters(int n, doublereal* const c) {
        setDensity(c[0]);
    }

    //! Get the equation of state parameters in a vector
    /*!
     * @internal
     * The number and meaning of these depends on the subclass.
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     */
    virtual void getParameters(int& n, doublereal* const c) const {
        double d = density();
        c[0] = d;
        n = 1;
    }

    //! Set equation of state parameter values from XML entries.
    /*!
     *
     * This method is called by function importPhase()  when processing a phase
     * definition in an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialized with elements and/or species.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);

protected:
    //! Temporary storage for dimensionless reference state enthalpies
    mutable vector_fp      m_h0_RT;

    //! Temporary storage for dimensionless reference state heat capacities
    mutable vector_fp      m_cp0_R;

    //! Temporary storage for dimensionless reference state Gibbs energies
    mutable vector_fp      m_g0_RT;

    //! Temporary storage for dimensionless reference state entropies
    mutable vector_fp      m_s0_R;

    //! Temporary array containing internally calculated partial pressures
    mutable vector_fp      m_pp;

    //! Current pressure (Pa)
    doublereal m_press;

private:
    //! Function to update the reference state thermo functions
    void _updateThermo() const;
};
}

#endif
