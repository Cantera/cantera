/**
 *  @file SingleSpeciesTP.h
 *  Header for the SingleSpeciesTP class, which is a filter class for ThermoPhase,
 *  that eases the construction of single species phases
 *  ( see \ref thermoprops and class \link Cantera::SingleSpeciesTP SingleSpeciesTP\endlink).
 *
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_SINGLESPECIESTP_H
#define CT_SINGLESPECIESTP_H

#include "ThermoPhase.h"


namespace Cantera
{

/**
 * @ingroup thermoprops
 *
 *  The SingleSpeciesTP class is a filter class for ThermoPhase.
 *  What it does is to simplify the construction of ThermoPhase
 *  objects by assuming that the phase consists of one and
 *  only one type of species. In other words, it's a stoichiometric
 *  phase. However, no assumptions are made concerning the
 *  thermodynamic functions or the equation of state of the
 *  phase. Therefore it's an incomplete description of
 *  the thermodynamics. The complete description must be
 *  made in a derived class of SingleSpeciesTP.
 *
 *  Several different groups of thermodynamic functions are resolved
 *  at this level by this class. For example, All partial molar property
 *  routines call their single species standard state equivalents.
 *  All molar solution thermodynamic routines call the single species
 *  standard state equivalents.
 *  Activities routines are resolved at this level, as there is only
 *  one species.
 *
 *  It is assumed that the reference state thermodynamics may be
 *  obtained by a pointer to a populated species thermodynamic property
 *  manager class (see ThermoPhase::m_spthermo). How to relate pressure
 *  changes to the reference state thermodynamics is again left open
 *  to implementation.
 *
 *  Mole fraction and Mass fraction vectors are assumed to be equal
 *  to x[0] = 1 y[0] = 1, respectively. Simplifications to the interface
 *  of setState_TPY() and setState_TPX() functions result and are made
 *  within the class.
 *
 *  Note, this class can handle the thermodynamic description of one
 *  phase of one species. It can not handle the description of phase
 *  equilibrium between two phases of a stoichiometric compound
 *  (e.g. water liquid and water vapor, below the critical point).
 *  However, it may be used to describe the thermodynamics of one phase
 *  of such a compound even past the phase equilibrium point, up to the
 *  point where the phase itself ceases to be a stable phase.
 *
 *  This class doesn't do much at the initialization level.
 *  Its SingleSpeciesTP::initThermo()
 *  member does check that one and only one species has been defined
 *  to occupy the phase.
 */
class SingleSpeciesTP : public ThermoPhase
{
public:
    //! Base empty constructor.
    SingleSpeciesTP();

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    SingleSpeciesTP(const SingleSpeciesTP&  right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    SingleSpeciesTP& operator=(const SingleSpeciesTP& right);

    //! Duplication function
    /*!
     * This virtual function is used to create a duplicate of the
     * current phase. It's used to duplicate the phase when given
     * a ThermoPhase pointer to the phase.
     *
     * @return It returns a ThermoPhase pointer.
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    /**
     * Returns the equation of state type flag.
     * This is a modified base class.
     * Therefore, if not overridden in derived classes,
     * this call will throw an exception.
     */
    virtual int eosType() const;

    /**
     *  @name  Molar Thermodynamic Properties of the Solution
     *
     *  These functions are resolved at this level, by reference
     *  to the partial molar functions and standard state
     *  functions for species 0. Derived classes don't need
     *  to supply entries for these functions.
     * @{
     */

    /// Molar enthalpy. Units: J/kmol.
    /*!
     * This function is resolved here by calling the standard state
     * thermo function.
     */
    doublereal enthalpy_mole() const;

    /// Molar internal energy. Units: J/kmol.
    /*!
     * This function is resolved here by calling the standard state
     * thermo function.
     */
    doublereal intEnergy_mole() const;

    /// Molar entropy. Units: J/kmol/K.
    /*!
     * This function is resolved here by calling the standard state
     * thermo function.
     */
    doublereal entropy_mole() const;

    /// Molar Gibbs function. Units: J/kmol.
    /*!
     * This function is resolved here by calling the standard state
     * thermo function.
     */
    doublereal gibbs_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    /*!
     * This function is resolved here by calling the standard state
     * thermo function.
     */
    doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    /*!
     * This function is resolved here by calling the standard state
     * thermo function.
     */
    doublereal cv_mole() const;

    /**
     * @}
     * @name Activities, Standard State, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature.
     * @{
     */

    /**
     * Get the array of non-dimensional activities at
     * the current solution temperature, pressure, and
     * solution concentration.
     *
     * We redefine this function to just return 1.0 here.
     *
     * @param a   Output vector of activities. Length: 1.
     */
    virtual void getActivities(doublereal* a) const {
        a[0] = 1.0;
    }

    /**
     * Get the array of non-dimensional activity coefficients at
     * the current solution temperature, pressure, and
     * solution concentration.
     *
     * @param ac Output vector of activity coefficients. Length: 1.
     */
    virtual void getActivityCoefficients(doublereal* ac) const {
        ac[0] = 1.0;
    }

    //@}
    /// @name  Partial Molar Properties of the Solution
    ///
    ///  These functions are resolved at this level, by reference
    ///  to the partial molar functions and standard state
    ///  functions for species 0. Derived classes don't need
    ///  to supply entries for these functions.
    //@{

    //!  Get the array of non-dimensional species chemical potentials
    //! These are partial molar Gibbs free energies.
    /*!
     *  These are the phase, partial molar, and the standard state
     *  dimensionless chemical potentials.
     *  \f$ \mu_k / \hat R T \f$.
     *
     * Units: unitless
     *
     *  @param murt   On return, Contains the chemical potential / RT of the single species
     *                and the phase. Units are unitless. Length = 1
     */
    void getChemPotentials_RT(doublereal* murt) const;

    //! Get the array of chemical potentials
    /*!
     * These are the phase, partial molar, and the standard state chemical potentials.
     *     \f$ \mu(T,P) = \mu^0_k(T,P) \f$.
     *
     *  @param mu   On return, Contains the chemical potential of the single species
     *              and the phase. Units are J / kmol . Length = 1
     */
    void getChemPotentials(doublereal* mu) const;

    //! Get the species electrochemical potentials. Units: J/kmol.
    /*!
     * This method adds a term \f$ Fz_k \phi_k \f$ to
     * each chemical potential.
     *
     * This is resolved here. A single  species phase
     * is not allowed to have anything other than a zero charge.
     *
     *  @param mu   On return, Contains the electrochemical potential of the single species
     *                and the phase. Units J/kmol . Length = 1
     */
    void getElectrochemPotentials(doublereal* mu) const;

    //!  Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * These are the phase enthalpies.  \f$ h_k \f$.
     *
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: 1. units are J/kmol.
     */
    void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Get the species partial molar internal energies. Units: J/kmol.
    /*!
     * These are the phase internal energies.  \f$ u_k \f$.
     *
     *  @param ubar On return, Contains the internal energy of the single species
     *              and the phase. Units are J / kmol . Length = 1
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Get the species partial molar entropy. Units: J/kmol K.
    /*!
     * This is the phase entropy.  \f$ s(T,P) = s_o(T,P) \f$.
     *
     *  @param sbar On return, Contains the entropy of the single species
     *              and the phase. Units are J / kmol / K . Length = 1
     */
    void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the species partial molar Heat Capacities. Units: J/ kmol /K.
    /*!
     * This is the phase heat capacity.  \f$ Cp(T,P) = Cp_o(T,P) \f$.
     *
     *  @param cpbar On return, Contains the heat capacity of the single species
     *              and the phase. Units are J / kmol / K . Length = 1
     */
    void getPartialMolarCp(doublereal* cpbar) const;

    //! Get the species partial molar volumes. Units: m^3/kmol.
    /*!
     * This is the phase molar volume.  \f$ V(T,P) = V_o(T,P) \f$.
     *
     *  @param vbar On return, Contains the molar volume of the single species
     *              and the phase. Units are m^3 / kmol. Length = 1
     */
    void getPartialMolarVolumes(doublereal* vbar) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    /// These functions are the primary way real properties are
    /// supplied to derived thermodynamics classes of SingleSpeciesTP.
    /// These functions must be supplied in derived classes. They
    /// are not resolved at the SingleSpeciesTP level.
    //@{

    /**
     * Get the dimensional Gibbs functions for the standard
     * state of the species at the current T and P.
     *
     * @param gpure returns a vector of size 1, containing the Gibbs function
     *              Units: J/kmol.
     */
    void getPureGibbs(doublereal* gpure) const;

    //! Get the molar volumes of each species in their standard
    //! states at the current  <I>T</I> and <I>P</I> of the solution.
    /*!
     *   units = m^3 / kmol
     *
     * We resolve this function at this level, by assigning
     * the molecular weight divided by the phase density
     *
     * @param vbar On output this contains the standard volume of the species
     *             and phase (m^3/kmol). Vector of length 1
     */
    void getStandardVolumes(doublereal* vbar) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference State
    ///
    /// Almost all functions in this group are resolved by this
    /// class. It is assumed that the m_spthermo species thermo
    /// pointer is populated and yields the reference state thermodynamics
    /// The internal energy function is not given by this
    /// class, since it would involve a specification of the
    /// equation of state.
    //@{

    /*!
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param hrt     Output vector containing the nondimensional reference state enthalpies
     *                Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /*!
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    /*!
     *  Returns the vector of the
     *  Gibbs function of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *  units = J/kmol
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param g       Output vector containing the  reference state
     *                Gibbs Free energies.  Length: m_kk. Units: J/kmol.
     */
    virtual void  getGibbs_ref(doublereal* g) const;

    /*!
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for each species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param er      Output vector containing the nondimensional reference state
     *                entropies.  Length: m_kk.
       */
    virtual void getEntropy_R_ref(doublereal* er) const;

    /*!
     *  Returns the vector of nondimensional
     *  constant pressure heat capacities of the reference state
     *  at the current temperature of the solution
     *  and reference pressure for each species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    /**
     * @name Setting the State
     *
     * These methods set all or part of the thermodynamic state.
     * @{
     */

    //! Mass fractions are fixed, with Y[0] = 1.0.
    void setMassFractions(const doublereal* const y) {};

    //! Mole fractions are fixed, with x[0] = 1.0.
    void setMoleFractions(const doublereal* const x) {};

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
    //@}

    /**
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.   When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase().
     *
     * Inheriting objects should call this function
     *
     * This version sets the mole fraction vector to x[0] = 1.0, and then
     * calls the ThermoPhase::initThermo() function.
     */
    virtual void initThermo();

protected:
    //! The current pressure of the solution (Pa)
    /*!
     * It gets initialized to 1 atm.
     */
    doublereal m_press;

    /*!
     * Reference pressure (Pa) must be the same for all species
     * - defaults to 1 atm.
     */
    doublereal m_p0;

    //! Dimensionless enthalpy at the (mtlast, m_p0)
    mutable vector_fp       m_h0_RT;
    //! Dimensionless heat capacity at the (mtlast, m_p0)
    mutable vector_fp       m_cp0_R;
    //! Dimensionless entropy at the (mtlast, m_p0)
    mutable vector_fp       m_s0_R;

protected:
    /**
     * @internal
     *        This crucial internal routine calls the species thermo
     *        update program to calculate new species Cp0, H0, and
     *        S0 whenever the temperature has changed.
     */
    void _updateThermo() const;
};

}

#endif
