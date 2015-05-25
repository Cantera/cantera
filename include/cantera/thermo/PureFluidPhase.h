/**
 *  @file PureFluidPhase.h
 *
 *   Header for a ThermoPhase class for a pure fluid phase consisting of
 *   gas, liquid, mixed-gas-liquid and supercrit fluid (see \ref thermoprops
 *   and class \link Cantera::PureFluidPhase PureFluidPhase\endlink).
 *
 * It inherits from ThermoPhase, but is built on top of the tpx package.
 */

//  Copyright 2003 California Institute of Technology

#ifndef CT_EOS_TPX_H
#define CT_EOS_TPX_H

#include "ThermoPhase.h"
#include "mix_defs.h"
#include "cantera/tpx/Sub.h"

namespace Cantera
{
//!   This phase object consists of a single component that can be a
//!   gas, a liquid, a mixed gas-liquid fluid, or a fluid beyond its
//!   critical point
/*!
 *  The object inherits from ThermoPhase. However, it's built on top
 *  of the tpx package.
 *
 * @ingroup thermoprops
 */
class PureFluidPhase  : public ThermoPhase
{
public:

    //! Empty Base Constructor
    PureFluidPhase();

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    PureFluidPhase(const PureFluidPhase& right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    PureFluidPhase& operator=(const PureFluidPhase& right);

    //! Destructor
    virtual ~PureFluidPhase();

    //! Duplication function
    /*!
     * This virtual function is used to create a duplicate of the
     * current phase. It's used to duplicate the phase when given
     * a ThermoPhase pointer to the phase.
     *
     * @return It returns a ThermoPhase pointer.
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Equation of state type
    virtual int eosType() const {
        return cPureFluid;
    }

    /// Molar enthalpy. Units: J/kmol.
    virtual doublereal enthalpy_mole() const;

    /// Molar internal energy. Units: J/kmol.
    virtual doublereal intEnergy_mole() const;

    /// Molar entropy. Units: J/kmol/K.
    virtual doublereal entropy_mole() const;

    /// Molar Gibbs function. Units: J/kmol.
    virtual doublereal gibbs_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    virtual doublereal cv_mole() const;

    //! Return the thermodynamic pressure (Pa).
    /*!
     * This method calculates the current pressure consistent with the
     * independent variables, T, rho.
     */
    virtual doublereal pressure() const;

    //! sets the thermodynamic pressure (Pa).
    /*!
     * This method calculates the density that is consistent with the
     * desired pressure, given the temperature.
     *
     * @param p  Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const {
        mu[0] = gibbs_mole();
    }

    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture. Units (J/kmol)
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Return an array of partial molar internal energies for the
    //! species in the mixture.  Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Return an array of partial molar heat capacities for the
    //! species in the mixture.  Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat
     *                capacities at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

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
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard concentration. The units are by definition
     *   dependent on the ThermoPhase and kinetics manager representation.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Get the array of non-dimensional activities at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * Note, for molality based formulations, this returns the
     * molality based activities.
     *
     * We resolve this function at this level by calling
     * on the activityConcentration function. However,
     * derived classes may want to override this default
     * implementation.
     *
     * @param a   Output vector of activities. Length: m_kk.
     */
    virtual void getActivities(doublereal* a) const;

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const;

    //! Returns a reference to the substance object
    tpx::Substance& TPX_Substance();

    //@}
    /// @name Properties of the Standard State of the Species in the Solution
    /*!
     *  The standard state of the pure fluid is defined as the real properties
     *  of the pure fluid at the most stable state of the fluid at the current
     *  temperature and pressure of the solution. With this definition, the
     *  activity of the fluid is always then defined to be equal to one.
     */
    //@{

    //! Get the array of chemical potentials at unit activity for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu      Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

    //! Get the nondimensional Enthalpy functions for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the
    //! standard state species at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Gibbs functions for the species
    //! in their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state Gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //@}

    /// @name Thermodynamic Values for the Species Reference States
    /*!
     *    The species reference state for pure fluids is defined as an ideal gas at the
     *    reference pressure and current temperature of the fluid.
     */
    //@{

    //!  Returns the vector of nondimensional enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference state enthalpies
     *                Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    //!  Returns the vector of nondimensional Gibbs Free Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    //!  Returns the vector of the Gibbs function of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     *  units = J/kmol
     *
     * @param g       Output vector containing the reference state
     *                Gibbs Free energies.  Length: m_kk. Units: J/kmol.
     */
    virtual void getGibbs_ref(doublereal* g) const;

    //!  Returns the vector of nondimensional entropies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param er      Output vector containing the nondimensional reference state
     *                entropies.  Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal* er) const;

    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic state.
     * @{
     */

    //! Set the internally stored specific enthalpy (J/kg) and pressure (Pa) of the phase.
    /*!
     * @param h     Specific enthalpy (J/kg)
     * @param p     Pressure (Pa)
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_HP(doublereal h, doublereal p,
                             doublereal tol = 1.e-8);

    //! Set the specific internal energy (J/kg) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific internal energy and specific volume have the value of the input parameters.
     *
     * @param u    specific internal energy (J/kg)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_UV(doublereal u, doublereal v,
                             doublereal tol = 1.e-8);

    //! Set the specific entropy (J/kg/K) and specific volume (m^3/kg).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific entropy and specific volume have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param v    specific volume (m^3/kg).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_SV(doublereal s, doublereal v,
                             doublereal tol = 1.e-8);

    //! Set the specific entropy (J/kg/K) and pressure (Pa).
    /*!
     * This function fixes the internal state of the phase so that
     * the specific entropy and the pressure have the value of the input parameters.
     *
     * @param s    specific entropy (J/kg/K)
     * @param p    specific pressure (Pa).
     * @param tol  Optional parameter setting the tolerance of the
     *             calculation.
     */
    virtual void setState_SP(doublereal s, doublereal p,
                             doublereal tol = 1.e-8);
    //@}

    //! @name Critical State Properties
    //@{

    //! critical temperature
    virtual doublereal critTemperature() const;

    //! critical pressure
    virtual doublereal critPressure() const;

    //! critical density
    virtual doublereal critDensity() const;

    //@}

    //! @name Saturation properties.
    //@{

    //! saturation temperature
    /*!
     * @param p  Pressure (Pa)
     */
    virtual doublereal satTemperature(doublereal p) const;

    //! Return the saturation pressure given the temperature
    /*!
     * @param t Temperature (Kelvin)
     */
    virtual doublereal satPressure(doublereal t);

    //! Return the fraction of vapor at the current conditions
    virtual doublereal vaporFraction() const;

    //! Set the state to a saturated system at a particular temperature
    /*!
     * @param t  Temperature (kelvin)
     * @param x  Fraction of vapor
     */
    virtual void setState_Tsat(doublereal t, doublereal x);

    //! Set the state to a saturated system at a particular pressure
    /*!
     * @param p  Pressure (Pa)
     * @param x  Fraction of vapor
     */
    virtual void setState_Psat(doublereal p, doublereal x);
    //@}

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

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() when processing a phase
     * definition in an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialized with elements and/or species.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);

    //! returns a summary of the state of the phase as a string
    /*!
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     * @param threshold   Unused in this subclass
     */
    virtual std::string report(bool show_thermo=true,
                               doublereal threshold=1e-14) const;

protected:

    //! Main call to the tpx level to set the state of the system
    /*!
     * @param n  Integer indicating which 2 thermo components are held constant
     * @param x  Value of the first component
     * @param y  Value of the second component
     */
    void Set(tpx::PropertyPair::type n, double x, double y) const;

    //! Sets the state using a TPX::TV call
    void setTPXState() const;

private:
    //! Pointer to the underlying tpx object Substance that does the work
    mutable tpx::Substance* m_sub;

    //! Int indicating the type of the fluid
    /*!
     * The tpx package uses an int to indicate what fluid is being sought.
     */
    int m_subflag;

    //! Molecular weight of the substance (kg kmol-1)
    doublereal m_mw;

    //! flag to turn on some printing.
    bool m_verbose;
};

}

#endif
