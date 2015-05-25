/**
 * @file MaskellSolidSolnPhase.h Header file for a solid solution model
 * following Maskell, Shaw, and Tye. Electrochimica Acta 1982
 *
 * This class inherits from the Cantera class ThermoPhase and implements a
 * non-ideal solid solution model with incompressible thermodynamics.
 */
/*
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 */

#ifndef CT_MASKELLSOLIDSOLNPHASE_H
#define CT_MASKELLSOLIDSOLNPHASE_H

#include "VPStandardStateTP.h"

namespace Cantera
{

/**
 * Class MaskellSolidSolnPhase represents a condensed phase
 * non-ideal solution with 2 species following the thermodynamic
 * model described in Maskell, Shaw, and Tye, Manganese Dioxide Electrode -- IX,
 * Electrochimica Acta 28(2) pp 231-235, 1983.
 *
 * @ingroup thermoprops
 */
class MaskellSolidSolnPhase : public VPStandardStateTP
{
public:
    MaskellSolidSolnPhase();

    //! Copy Constructor
    MaskellSolidSolnPhase(const MaskellSolidSolnPhase&);

    //! Assignment operator
    MaskellSolidSolnPhase& operator=(const MaskellSolidSolnPhase&);

    /*!
     * Base Class Duplication Function
     *
     * Given a pointer to ThermoPhase, this function can duplicate the object.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    /**
     * This method returns the array of generalized
     * concentrations. The generalized concentrations are used
     * in the evaluation of the rates of progress for reactions
     * involving species in this phase. The generalized
     * concentration divided by the standard concentration is also
     * equal to the activity of species.
     *
     * @param c  Pointer to array of doubles of length m_kk, which on exit
     *           will contain the generalized concentrations.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize the
     * generalized concentration.
     *
     * @param k Species number: this is an optional parameter,
     */
    virtual doublereal standardConcentration(size_t k=0) const { return 1.0; }

    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
     */
    virtual doublereal logStandardConc(size_t k=0) const { return 0.0; }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{
    /**
     * Molar enthalpy of the solution. Units: J/kmol.
     */
    virtual doublereal enthalpy_mole() const;

    /**
     * Molar entropy of the solution. Units: J/kmol/K.
     */
    virtual doublereal entropy_mole() const;

    /**
     * Molar heat capacity at constant pressure of the solution.
     * Units: J/kmol/K.
     */
    //virtual doublereal cp_mole() const;

    /**
     * Molar heat capacity at constant volume of the solution.
     * Units: J/kmol/K.
     */
    //virtual doublereal cv_mole() const {
    //    return cp_mole();
    //}

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

    virtual void calcDensity();

    /**
     * Overwritten setMolarDensity() function is necessary because the
     * density is not an independent variable.
     *
     * This function will now throw an error condition.
     *
     * @param rho   Input Density
     */
    virtual void setMolarDensity(const doublereal rho);

    //@}

    /**
     * @name Chemical Potentials and Activities
     * @{
     */

    //! Get the array of species activity coefficients
    /*!
     * @param ac output vector of activity coefficients. Length: m_kk
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

    /**
     * Get the species chemical potentials. Units: J/kmol.
     *
     * @param mu  Output vector of chemical potentials.
     */
    virtual void getChemPotentials(doublereal* mu) const;

    /**
     * Get the array of non-dimensional species solution
     * chemical potentials at the current T and P
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
     *
     * @param hbar Output vector containing partial molar enthalpies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K.
     *
     * @param sbar Output vector containing partial molar entropies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    /**
     * Returns an array of partial molar Heat Capacities at constant
     * pressure of the species in the
     * solution. Units: J/kmol/K.
     *
     * @param cpbar  Output vector of partial heat capacities. Length: m_kk.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    /**
     * returns an array of partial molar volumes of the species
     * in the solution. Units: m^3 kmol-1.
     *
     * @param vbar  Output vector of partial molar volumes. Length: m_kk.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of standard state Gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

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

    //@}
    /// @name Utility Functions
    //@{

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


    void set_h_mix(const doublereal hmix) { h_mixing = hmix; }
    //@}

private:
    /**
     * m_Pcurrent = The current pressure
     * Since the density isn't a function of pressure, but only of the
     * mole fractions, we need to independently specify the pressure.
     */
    doublereal m_Pcurrent;

    /**
     * Function to call through to m_spthermo->update and fill m_h0_RT,
     * m_cp0_R, m_g0_RT, m_s0_R.
     */
    void _updateThermo() const;

    //! Vector containing the species reference enthalpies at T = m_tlast
    mutable vector_fp m_h0_RT;

    /**
     * Vector containing the species reference constant pressure
     * heat capacities at T = m_tlast
     */
    mutable vector_fp m_cp0_R;

    //!  Vector containing the species reference Gibbs functions at T = m_tlast
    mutable vector_fp m_g0_RT;

    //! Vector containing the species reference entropies at T = m_tlast
    mutable vector_fp m_s0_R;

    //! Value of the enthalpy change on mixing due to protons changing from type B to type A configurations.
    doublereal h_mixing;

    //! Index of the species whose mole fraction defines the extent of reduction r
    int product_species_index;
    int reactant_species_index;

    // Functions to calculate some of the pieces of the mixing terms.
    doublereal s() const;
    doublereal fm(const doublereal r) const;
    doublereal p(const doublereal r) const;
};
}

#endif
