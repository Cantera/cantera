/**
 * @file IdealSolidSolnPhase.h Header file for an ideal solid solution model
 *      with incompressible thermodynamics (see \ref thermoprops and
 *      \link Cantera::IdealSolidSolnPhase IdealSolidSolnPhase\endlink).
 *
 * This class inherits from the Cantera class ThermoPhase and implements an
 * ideal solid solution model with incompressible thermodynamics.
 */
/*
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 */

#ifndef CT_IDEALSOLIDSOLNPHASE_H
#define CT_IDEALSOLIDSOLNPHASE_H

#include "mix_defs.h"
#include "ThermoPhase.h"
#include "ThermoFactory.h"
#include "SpeciesThermo.h"

namespace Cantera
{

/*!
 * @name CONSTANTS - Models for the Standard State of IdealSolidSolnPhase's
 */
//@{
const int cIdealSolidSolnPhase0 = 5010;
const int cIdealSolidSolnPhase1 = 5011;
const int cIdealSolidSolnPhase2 = 5012;
//@}

/**
 * Class IdealSolidSolnPhase represents a condensed phase ideal
 * solution compound. The phase and the pure species phases which
 * comprise the standard states of the species are assumed to have
 * zero volume expansivity and zero isothermal compressibility.
 * Each species does, however, have constant but distinct partial
 * molar volumes equal to their pure species molar volumes.
 * The class derives from class ThermoPhase,
 * and overloads the virtual methods defined there with ones that
 * use expressions appropriate for ideal solution mixtures.
 *
 * The generalized concentrations can have three different forms
 * depending on the value of the member attribute #m_formGC, which
 * is supplied in the constructor and in the XML file.
 * The value and form of the generalized concentration will affect
 * reaction rate constants involving species in this phase.
 *
 * @ingroup thermoprops
 */
class IdealSolidSolnPhase : public ThermoPhase
{
public:
    /**
     * Constructor for IdealSolidSolnPhase.
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute #m_formGC, which
     * is supplied in the constructor or read from the xml data file.
     *
     * @param formCG This parameter initializes the #m_formGC variable.
     */
    IdealSolidSolnPhase(int formCG=0);

    //! Construct and initialize an IdealSolidSolnPhase ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * This constructor will also fully initialize the object.
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute #m_formGC, which
     * is supplied in the constructor or read from the xml data file.
     *
     * @param infile File name for the XML datafile containing information
     *               for this phase
     * @param id     The name of this phase. This is used to look up
     *               the phase in the XML datafile.
     * @param formCG This parameter initializes the #m_formGC variable.
     */
    IdealSolidSolnPhase(const std::string& infile, const std::string& id="", int formCG=0);

    //! Construct and initialize an IdealSolidSolnPhase ThermoPhase object
    //! directly from an XML database
    /*!
     * The generalized concentrations can have three different forms
     * depending on the value of the member attribute #m_formGC, which
     * is supplied in the constructor and/or read from the data file.
     *
     * @param root   XML tree containing a description of the phase.
     *               The tree must be positioned at the XML element
     *               named phase with id, "id", on input to this routine.
     * @param id     The name of this phase. This is used to look up
     *               the phase in the XML datafile.
     * @param formCG This parameter initializes the #m_formGC variable.
     */
    IdealSolidSolnPhase(XML_Node& root, const std::string& id="", int formCG=0);

    //! Copy Constructor
    IdealSolidSolnPhase(const IdealSolidSolnPhase&);

    //! Assignment operator
    IdealSolidSolnPhase& operator=(const IdealSolidSolnPhase&);

    /*!
     * Base Class Duplication Function
     *
     * Given a pointer to ThermoPhase, this function can duplicate the object.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    /**
     * Equation of state flag. Returns a value depending upon the value of
     * #m_formGC, which is defined at instantiation.
     */
    virtual int eosType() const;

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    /**
     * Molar enthalpy of the solution. Units: J/kmol.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity and
     * zero isothermal compressibility:
     * \f[
     * \hat h(T,P) = \sum_k X_k \hat h^0_k(T) + (P - P_{ref}) (\sum_k X_k \hat V^0_k)
     * \f]
     * The reference-state pure-species enthalpies at the reference pressure Pref
     * \f$ \hat h^0_k(T) \f$, are computed by the species thermodynamic
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     */
    virtual doublereal enthalpy_mole() const;

    /**
     * Molar internal energy of the solution. Units: J/kmol.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity and
     * zero isothermal compressibility:
     * \f[
     * \hat u(T,X) = \hat h(T,P,X) - p \hat V
     *         =  \sum_k X_k \hat h^0_k(T)  - P_{ref} (\sum_k{X_k \hat V^0_k})
     * \f]
     * and is a function only of temperature.
     * The reference-state pure-species enthalpies
     * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal intEnergy_mole() const;

    /**
     * Molar entropy of the solution. Units: J/kmol/K.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     * \f[
     * \hat s(T, P, X_k) = \sum_k X_k \hat s^0_k(T)  - \hat R \sum_k X_k log(X_k)
     * \f]
     * The reference-state pure-species entropies
     * \f$ \hat s^0_k(T,p_{ref}) \f$ are computed by the species thermodynamic
     * property manager. The pure species entropies are independent of
     * pressure since the volume expansivities are equal to zero.
     * @see SpeciesThermo
     */
    virtual doublereal entropy_mole() const;

    /**
     * Molar gibbs free energy of the solution. Units: J/kmol.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     * \f[
     * \hat g(T, P) = \sum_k X_k \hat g^0_k(T,P) + \hat R T \sum_k X_k log(X_k)
     * \f]
     * The reference-state pure-species gibbs free energies
     * \f$ \hat g^0_k(T) \f$ are computed by the species thermodynamic
     * property manager, while the standard state gibbs free energies
     * \f$ \hat g^0_k(T,P) \f$ are computed by the member function, gibbs_RT().
     * @see SpeciesThermo
     */
    virtual doublereal gibbs_mole() const;

    /**
     * Molar heat capacity at constant pressure of the solution.
     * Units: J/kmol/K.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     * \f[
     * \hat c_p(T,P) = \sum_k X_k \hat c^0_{p,k}(T) .
     * \f]
     * The heat capacity is independent of pressure.
     * The reference-state pure-species heat capacities
     * \f$ \hat c^0_{p,k}(T) \f$ are computed by the species thermodynamic
     * property manager.
     * @see SpeciesThermo
     */
    virtual doublereal cp_mole() const;

    /**
     * Molar heat capacity at constant volume of the solution.
     * Units: J/kmol/K.
     * For an ideal, constant partial molar volume solution mixture with
     * pure species phases which exhibit zero volume expansivity:
     * \f[ \hat c_v(T,P) = \hat c_p(T,P) \f]
     * The two heat capacities are equal.
     */
    virtual doublereal cv_mole() const {
        return cp_mole();
    }

    //@}
    /** @name Mechanical Equation of State Properties
     *
     *   In this equation of state implementation, the density is a
     *   function only of the mole fractions. Therefore, it can't be
     *   an independent variable. Instead, the pressure is used as the
     *   independent variable. Functions which try to set the thermodynamic
     *   state by calling setDensity() may cause an exception to be
     *   thrown.
     */
    //@{

    /**
     * Pressure. Units: Pa.
     * For this incompressible system, we return the internally stored
     * independent value of the pressure.
     */
    virtual doublereal pressure() const {
        return m_Pcurrent;
    }

    /**
     * Set the pressure at constant temperature. Units: Pa.
     * This method sets a constant within the object.
     * The mass density is not a function of pressure.
     *
     * @param p   Input Pressure (Pa)
     */
    virtual void setPressure(doublereal p);

    /**
     * Calculate the density of the mixture using the partial
     * molar volumes and mole fractions as input
     *
     * The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are
     * the molecular weights, and \f$V_k\f$ are the pure species
     * molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal
     * solution the partial molar volumes are equal to the pure
     * species molar volumes. We have additionally specified
     * in this class that the pure species molar volumes are
     * independent of temperature and pressure.
     */
    void calcDensity();

    /**
     * Overwritten setDensity() function is necessary because the
     * density is not an independent variable.
     *
     * This function will now throw an error condition
     *
     * @internal May have to adjust the strategy here to make
     * the eos for these materials slightly compressible, in order
     * to create a condition where the density is a function of
     * the pressure.
     *
     * @param rho  Input density
     */
    virtual void setDensity(const doublereal rho);

    /**
     * Overwritten setMolarDensity() function is necessary because the
     * density is not an independent variable.
     *
     * This function will now throw an error condition.
     *
     * @param rho   Input Density
     */
    virtual void setMolarDensity(const doublereal rho);

    //! Set the mole fractions
    /*!
     * @param x  Input vector of mole fractions.
     *           Length: m_kk.
     */
    virtual void setMoleFractions(const doublereal* const x);

    //! Set the mole fractions, but don't normalize them to one.
    /*!
     * @param x  Input vector of mole fractions.
     *           Length: m_kk.
     */
    virtual void setMoleFractions_NoNorm(const doublereal* const x);

    //! Set the mass fractions, and normalize them to one.
    /*!
     * @param y  Input vector of mass fractions.
     *           Length: m_kk.
     */
    virtual void setMassFractions(const doublereal* const y);

    //! Set the mass fractions, but don't normalize them to one
    /*!
     * @param y  Input vector of mass fractions.
     *           Length: m_kk.
     */
    virtual void setMassFractions_NoNorm(const doublereal* const y);

    //! Set the concentration,
    /*!
     * @param c  Input vector of concentrations.
     *           Length: m_kk.
     */
    virtual void setConcentrations(const doublereal* const c);

    //@}

    /**
     * @name Chemical Potentials and Activities
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by
     * \f[
     *  \mu_k(T,P,X_k) = \mu_k^0(T,P)
     * + \hat R T \log a_k.
     *  \f]
     * The quantity \f$\mu_k^0(T,P)\f$ is
     * the standard state chemical potential at unit activity.
     * It may depend on the pressure and the temperature. However,
     * it may not depend on the mole fractions of the species
     * in the solid solution.
     *
     * The activities are related to the generalized
     * concentrations, \f$\tilde C_k\f$, and standard
     * concentrations, \f$C^0_k\f$, by the following formula:
     *
     *  \f[
     *  a_k = \frac{\tilde C_k}{C^0_k}
     *  \f]
     * The generalized concentrations are used in the kinetics classes
     * to describe the rates of progress of reactions involving the
     * species. Their formulation depends upon the specification
     * of the rate constants for reaction, especially the units used
     * in specifying the rate constants. The bridge between the
     * thermodynamic equilibrium expressions that use a_k and the
     * kinetics expressions which use the generalized concentrations
     * is provided by the multiplicative factor of the
     * standard concentrations.
     * @{
     */

    /**
     * This method returns the array of generalized
     * concentrations. The generalized concentrations are used
     * in the evaluation of the rates of progress for reactions
     * involving species in this phase. The generalized
     * concentration divided by the standard concentration is also
     * equal to the activity of species.
     *
     * For this implementation the activity is defined to be the
     * mole fraction of the species. The generalized concentration
     * is defined to be equal to the mole fraction divided by
     * the partial molar volume. The generalized concentrations
     * for species in this phase therefore have units of
     * kmol m<SUP>-3</SUP>. Rate constants must reflect this fact.
     *
     * On a general note, the following must be true.
     * For an ideal solution, the generalized concentration  must consist
     * of the mole fraction multiplied by a constant. The constant may be
     * fairly arbitrarily chosen, with differences adsorbed into the
     * reaction rate expression. 1/V_N, 1/V_k, or 1 are equally good,
     * as long as the standard concentration is adjusted accordingly.
     * However, it must be a constant (and not the concentration, btw,
     * which is a function of the mole fractions) in order for the
     * ideal solution properties to hold at the same time having the
     * standard concentration to be independent of the mole fractions.
     *
     * In this implementation the form of the generalized concentrations
     * depend upon the member attribute, #m_formGC.
     *
     * HKM Note: We have absorbed the pressure dependence of the pure species
     *        state into the thermodynamics functions. Therefore the
     *        standard state on which the activities are based depend
     *        on both temperature and pressure. If we hadn't, it would have
     *        appeared in this function in a very awkward exp[] format.
     *
     * @param c  Pointer to array of doubles of length m_kk, which on exit
     *           will contain the generalized concentrations.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    /**
     * The standard concentration \f$ C^0_k \f$ used to normalize the
     * generalized concentration. In many cases, this quantity will be the
     * same for all species in a phase. However, for this case, we will return
     * a distinct concentration for each species. This is the inverse of the
     * species molar volume. Units for the standard concentration are kmol
     * m<SUP>-3</SUP>.
     *
     * @param k Species number: this is a require parameter,
     * a change from the ThermoPhase base class, where it was
     * an optional parameter.
     */
    virtual doublereal standardConcentration(size_t k) const;

    /**
     * The reference (ie standard) concentration \f$ C^0_k \f$ used to normalize
     * the generalized concentration. In many cases, this quantity
     * will be the same for all species in a phase.
     * However, for this case, we will return a distinct concentration
     * for each species. (clone of the standard concentration ->
     * suggest changing the name). This is the inverse of the species molar
     * volume.
     *
     * @param k  Species index.
     */
    virtual doublereal referenceConcentration(int k) const;

    /**
     * Returns the log of the standard concentration of the kth species
     *
     * @param k Species number: this is a require parameter,
     *          a change from the ThermoPhase base class, where it was
     *          an optional parameter.
     */
    virtual doublereal logStandardConc(size_t k) const;

    /**
     * Returns the units of the standard and general concentrations
     * Note they have the same units, as their divisor is
     * defined to be equal to the activity of the kth species
     * in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     * @param uA Output vector containing the units:
     *
     *     uA[0] = kmol units - default  = 1
     *     uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                   dimensions in the Phase class.
     *     uA[2] = kg   units - default  = 0;
     *     uA[3] = Pa(pressure) units - default = 0;
     *     uA[4] = Temperature units - default = 0;
     *     uA[5] = time units - default = 0
     *
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     *
     *  For EOS types other than cIdealSolidSolnPhase0, the default
     *  kmol/m3 holds for standard concentration units. For
     *  cIdealSolidSolnPhase0 type, the standard concentration is
     *  unitless.
     * @deprecated
     */
    virtual void getUnitsStandardConc(double* uA, int k = 0,
                                      int sizeUA = 6) const;

    //! Get the array of species activity coefficients
    /*!
     * @param ac output vector of activity coefficients. Length: m_kk
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    /**
     * Get the species chemical potentials. Units: J/kmol.
     *
     * This function returns a vector of chemical potentials of the
     * species in solution.
     * \f[
     *    \mu_k = \mu^{ref}_k(T) + V_k * (p - p_o) + R T ln(X_k)
     * \f]
     *  or another way to phrase this is
     * \f[
     *    \mu_k = \mu^o_k(T,p) + R T ln(X_k)
     * \f]
     *  where \f$ \mu^o_k(T,p) = \mu^{ref}_k(T) + V_k * (p - p_o)\f$
     *
     * @param mu  Output vector of chemical potentials.
     */
    virtual void getChemPotentials(doublereal* mu) const;

    /**
     * Get the array of non-dimensional species solution
     * chemical potentials at the current T and P
     * \f$\mu_k / \hat R T \f$.
     * \f[
     *  \mu^0_k(T,P) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k + RT ln(X_k)
     * \f]
     * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
     * \f$ \mu^{ref}_k(T)\f$ is the chemical potential of pure
     * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
     *
     * @param mu   Output vector of dimensionless chemical potentials. Length = m_kk.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Returns an array of partial molar enthalpies for the species in the mixture.
    /*!
     * Units (J/kmol)
     * For this phase, the partial molar enthalpies are equal to the
     * pure species enthalpies
     *  \f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * \f]
     * The reference-state pure-species enthalpies, \f$ \hat h^{ref}_k(T) \f$,
     * at the reference pressure,\f$ P_{ref} \f$,
     * are computed by the species thermodynamic
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     *
     * @param hbar Output vector containing partial molar enthalpies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K.
     * For this phase, the partial molar entropies are equal to the
     * pure species entropies plus the ideal solution contribution.
     *  \f[
     * \bar s_k(T,P) =  \hat s^0_k(T) - R log(X_k)
     * \f]
     * The reference-state pure-species entropies,\f$ \hat s^{ref}_k(T) \f$,
     * at the reference pressure, \f$ P_{ref} \f$,  are computed by the
     * species thermodynamic
     * property manager. They are polynomial functions of temperature.
     * @see SpeciesThermo
     *
     * @param sbar Output vector containing partial molar entropies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    /**
     * Returns an array of partial molar Heat Capacities at constant
     * pressure of the species in the
     * solution. Units: J/kmol/K.
     * For this phase, the partial molar heat capacities are equal
     * to the standard state heat capacities.
     *
     * @param cpbar  Output vector of partial heat capacities. Length: m_kk.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    /**
     * returns an array of partial molar volumes of the species
     * in the solution. Units: m^3 kmol-1.
     *
     * For this solution, thepartial molar volumes are equal to the
     * constant species molar volumes.
     *
     * @param vbar  Output vector of partial molar volumes. Length: m_kk.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

    /**
     *  Get the standard state chemical potentials of the species.
     *  This is the array of chemical potentials at unit activity
     *  \f$ \mu^0_k(T,P) \f$.
     *  We define these here as the chemical potentials of the pure
     *  species at the temperature and pressure of the solution.
     *  This function is used in the evaluation of the
     *  equilibrium constant Kc. Therefore, Kc will also depend
     *  on T and P. This is the norm for liquid and solid systems.
     *
     *  units = J / kmol
     *
     * @param mu0   Output vector of standard state chemical potentials.
     *              Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const {
        getPureGibbs(mu0);
    }

    //! Get the array of nondimensional Enthalpy functions for the standard state species
    //! at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * We assume an incompressible constant partial molar
     * volume here:
     * \f[
     *  h^0_k(T,P) = h^{ref}_k(T) + (P - P_{ref}) * V_k
     * \f]
     * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
     * \f$ h^{ref}_k(T)\f$ is the enthalpy of the pure
     * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
     *
     * @param hrt Vector of length m_kk, which on return hrt[k]
     *            will contain the nondimensional
     *            standard state enthalpy of species k.
     */
    void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the nondimensional Entropies for the species
    //! standard states at the current T and P of the solution.
    /*!
     * Note, this is equal to the reference state entropies
     * due to the zero volume expansivity:
     * i.e., (dS/dP)_T = (dV/dT)_P = 0.0
     *
     * @param sr Vector of length m_kk, which on return sr[k]
     *           will contain the nondimensional
     *           standard state entropy for species k.
     */
    void getEntropy_R(doublereal* sr) const;

    /**
     * Get the nondimensional gibbs function for the species
     * standard states at the current T and P of the solution.
     *
     *  \f[
     *  \mu^0_k(T,P) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
     * \f]
     * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
     * \f$ \mu^{ref}_k(T)\f$ is the chemical potential of pure
     * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
     *
     * @param grt Vector of length m_kk, which on return sr[k]
     *           will contain the nondimensional
     *           standard state gibbs function for species k.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    /**
     * Get the Gibbs functions for the pure species
     * at the current <I>T</I> and <I>P</I> of the solution.
     * We assume an incompressible constant partial molar
     * volume here:
     * \f[
     *  \mu^0_k(T,P) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
     * \f]
     * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
     * \f$ \mu^{ref}_k(T)\f$ is the chemical potential of pure
     * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
     *
     * @param gpure  Output vector of Gibbs functions for species
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    //! Returns the vector of nondimensional
    //!  internal Energies of the standard state at the current
    //! temperature and pressure of the solution for each species.
    /*!
     * @param urt  Output vector of standard state nondimensional internal energies.
     *             Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    /**
     * Get the nondimensional heat capacity at constant pressure function for
     * the species standard states at the current T and P of the solution.
     * \f[
     *  Cp^0_k(T,P) = Cp^{ref}_k(T)
     * \f]
     * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
     * \f$ Cp^{ref}_k(T)\f$ is the constant pressure heat capacity
     * of species <I>k</I> at the reference pressure, \f$p_{ref}\f$.
     *
     * @param cpr Vector of length m_kk, which on return cpr[k]
     *           will contain the nondimensional
     *           constant pressure heat capacity for species k.
     */
    void getCp_R(doublereal* cpr) const;

    /**
     * Get the molar volumes of each species in their standard
     * states at the current <I>T</I> and <I>P</I> of the solution.
     * units = m^3 / kmol
     *
     * @param vol  Output vector of standard state volumes.
     *             Length: m_kk.
     */
    virtual void getStandardVolumes(doublereal* vol) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    /**
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param hrt  Output vector containing reference nondimensional enthalpies.
     *             Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /**
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param grt  Output vector containing reference nondimensional Gibbs free energies.
     *             Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    /**
     *  Returns the vector of the
     *  gibbs function of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *  units = J/kmol
     *
     * @param g    Output vector containing reference Gibbs free energies.
     *             Length: m_kk.
     */
    virtual void getGibbs_ref(doublereal* g) const;

    /**
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param er  Output vector containing reference nondimensional entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R_ref(doublereal* er) const;

    /**
     *  Returns the vector of nondimensional
     *  internal Energies of the reference state at the current temperature
     *  of the solution and the reference pressure for each species.
     *
     * @param urt  Output vector containing reference nondimensional internal energies.
     *             Length: m_kk.
     */
    virtual void getIntEnergy_RT_ref(doublereal* urt) const;

    /**
     *  Returns the vector of nondimensional
     *  constant pressure heat capacities of the reference state
     *  at the current temperature of the solution
     *  and reference pressure for the species.
     *
     * @param cprt  Output vector containing reference nondimensional heat capacities.
     *             Length: m_kk.
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    /**
     *  Returns a reference to the vector of nondimensional
     *  enthalpies of the reference state at the current temperature.
     *  Real reason for its existence is that it also checks
     *  to see if a recalculation of the reference thermodynamics
     *  functions needs to be done.
     */
    const vector_fp& enthalpy_RT_ref() const;

    /**
     *  Returns a reference to the vector of nondimensional
     *  enthalpies of the reference state at the current temperature.
     *  Real reason for its existence is that it also checks
     *  to see if a recalculation of the reference thermodynamics
     *  functions needs to be done.
     */
    const vector_fp& gibbs_RT_ref() const {
        _updateThermo();
        return m_g0_RT;
    }

    /**
     *  Returns a reference to the vector of nondimensional
     *  enthalpies of the reference state at the current temperature.
     *  Real reason for its existence is that it also checks
     *  to see if a recalculation of the reference thermodynamics
     *  functions needs to be done.
     */
    const vector_fp& entropy_R_ref() const;

    /**
     *  Returns a reference to the vector of nondimensional
     *  enthalpies of the reference state at the current temperature.
     *  Real reason for its existence is that it also checks
     *  to see if a recalculation of the reference thermodynamics
     *  functions needs to be done.
     */
    const vector_fp& cp_R_ref() const {
        _updateThermo();
        return m_cp0_R;
    }

    virtual void setPotentialEnergy(int k, doublereal pe) {
        m_pe[k] = pe;
        _updateThermo();
    }

    virtual doublereal potentialEnergy(int k) const {
        return m_pe[k];
    }

    //@}
    /// @name Utility Functions
    //@{

    /**
     *  Initialization of an IdealSolidSolnPhase phase:
     *  Note this function is pretty much useless because it doesn't
     *  get the xml tree passed to it. Suggest a change.
     */
    virtual void initThermo();

    /**
      * @internal Import and initialize a ThermoPhase object using an XML
      *   tree. Here we read extra information about the XML description of a
      *   phase. Regular information about elements and species and their
      *   reference state thermodynamic information have already been read at
      *   this point. For example, we do not need to call this function for
      *   ideal gas equations of state. This function is called from
      *   importPhase() after the elements and the species are initialized
      *   with default ideal solution level data.
      *
      * @param phaseNode This object must be the phase node of a complete XML
      *             tree description of the phase, including all of the
      *             species data. In other words while "phase" must point to
      *             an XML phase object, it must have sibling nodes
      *             "speciesData" that describe the species in the phase.
      * @param id   ID of the phase. If nonnull, a check is done to see if
      *             phaseNode is pointing to the phase with the correct id.
      */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    /**
     * Set mixture to an equilibrium state consistent with specified
     * element potentials and the temperature.
     *
     * @param lambda_RT vector of non-dimensional element potentials
     * \f$ \lambda_m/RT \f$.
     */
    virtual void setToEquilState(const doublereal* lambda_RT);

    /**
     * Report the molar volume of species k
     *
     * units - \f$ m^3 kmol^-1 \f$
     *
     * @param  k species index
     */
    double speciesMolarVolume(int k) const;

    /**
     * Fill in a return vector containing the species molar volumes.
     *
     * units - \f$ m^3 kmol^-1 \f$
     *
     * @param smv  output vector containing species molar volumes.
     *             Length: m_kk.
     */
    void getSpeciesMolarVolumes(doublereal* smv) const;

    //@}

protected:
    /**
     *  Format for the generalized concentrations.
     *
     *  <TABLE>
     *  <TR><TD> m_formGC    </TD><TD> GeneralizedConc </TD><TD> StandardConc </TD></TR>
     *  <TR><TD> 0 (default) </TD><TD> X_k             </TD><TD> 1.0          </TD></TR>
     *  <TR><TD> 1           </TD><TD> X_k / V_k       </TD><TD> 1.0 / V_k    </TD></TR>
     *  <TR><TD> 2           </TD><TD> X_k / V_N       </TD><TD> 1.0 / V_N    </TD></TR>
     *  </TABLE>
     *
     *  The value and form of the generalized concentration will affect
     *  reaction rate constants involving species in this phase.
     */
    int m_formGC;

    /**
     * Value of the reference pressure for all species in this phase.
     * The T dependent polynomials are evaluated at the reference
     * pressure. Note, because this is a single value, all species
     * are required to have the same reference pressure.
     */
    doublereal m_Pref;

    /**
     * m_Pcurrent = The current pressure
     * Since the density isn't a function of pressure, but only of the
     * mole fractions, we need to independently specify the pressure.
     * The density variable which is inherited as part of the State class,
     * m_dens, is always kept current whenever T, P, or X[] change.
     */
    doublereal m_Pcurrent;

    //! Vector of molar volumes for each species in the solution
    /**
     * Species molar volumes \f$ m^3 kmol^-1 \f$
     */
    vector_fp   m_speciesMolarVolume;

    /**
     *  Value of the temperature at which the thermodynamics functions
     * for the reference state of the species were last evaluated.
     */
    mutable doublereal   m_tlast;

    //! Vector containing the species reference enthalpies at T = m_tlast
    mutable vector_fp      m_h0_RT;

    /**
     * Vector containing the species reference constant pressure
     * heat capacities at T = m_tlast
     */
    mutable vector_fp      m_cp0_R;

    //!  Vector containing the species reference Gibbs functions at T = m_tlast
    mutable vector_fp      m_g0_RT;

    //! Vector containing the species reference entropies at T = m_tlast
    mutable vector_fp      m_s0_R;

    /**
     * Vector containing the species reference exp(-G/RT) functions
     * at T = m_tlast
     */
    mutable vector_fp      m_expg0_RT;

    //! Vector of potential energies for the species.
    mutable vector_fp      m_pe;

    //! Temporary array used in equilibrium calculations
    mutable vector_fp      m_pp;

private:
    /// @name Utility Functions
    //@{
    /**
     * This function gets called for every call to functions in this
     * class. It checks to see whether the temperature has changed and
     * thus the reference thermodynamics functions for all of the species
     * must be recalculated.
     * If the temperature has changed, the species thermo manager is called
     * to recalculate G, Cp, H, and S at the current temperature.
     */
    void _updateThermo() const;

    //! This internal function adjusts the lengths of arrays
    void initLengths();

    //@}
};
}

#endif
