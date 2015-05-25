/**
 *  @file StoichSubstance.h
 * This file contains the class declarations for the StoichSubstance
 * ThermoPhase class.
 */

//  Copyright 2001 California Institute of Technology

#ifndef CT_STOICHSUBSTANCE_H
#define CT_STOICHSUBSTANCE_H

#include "mix_defs.h"
#include "ThermoPhase.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 *
 * Class StoichSubstance represents a stoichiometric (fixed composition)
 * incompressible substance.
 */
class StoichSubstance : public ThermoPhase
{
public:
    //! Default empty constructor
    StoichSubstance();

    //! Copy Constructor
    /*!
     * Copy constructor for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     * This is a wrapper around the assignment operator
     *
     * @param right Object to be copied.
     */
    StoichSubstance(const StoichSubstance& right);

    //! Assignment operator
    /*!
     * Assignment operator for the object. Constructed
     * object will be a clone of this object, but will
     * also own all of its data.
     *
     * @param right Object to be copied.
     */
    StoichSubstance& operator=(const StoichSubstance& right);

    //! Duplicator from the ThermoPhase parent class
    /*
     * Given a pointer to a ThermoPhase object, this function will
     * duplicate the ThermoPhase object and all underlying structures.
     * This is basically a wrapper around the copy constructor.
     *
     * @return returns a pointer to a ThermoPhase
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    /**
     * Equation of state flag. Returns the value cStoichSubstance,
     * defined in mix_defs.h.
     */
    virtual int eosType() const {
        return cStoichSubstance;
    }

    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    /**
     * Molar enthalpy. Units: J/kmol.  For an incompressible,
     * stoichiometric substance, the internal energy is
     * independent of pressure, and therefore the molar enthalpy
     * is \f[ \hat h(T, P) = \hat u(T) + P \hat v \f], where the
     * molar specific volume is constant.
     */
    virtual doublereal enthalpy_mole() const;

    /**
     * Molar internal energy. J/kmol.  For an incompressible,
     * stoichiometric substance, the molar internal energy is
     * independent of pressure. Since the thermodynamic properties
     * are specified by giving the standard-state enthalpy, the
     * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
     * enthalpy to compute the molar internal energy.
     */
    virtual doublereal intEnergy_mole() const;

    /**
     * Molar entropy. Units: J/kmol/K.  For an incompressible,
     * stoichiometric substance, the molar entropy depends only on
     * the temperature.
     */
    virtual doublereal entropy_mole() const;

    /**
     * Molar heat capacity at constant pressure. Units: J/kmol/K.
     * For an incompressible substance, \f$ \hat c_p = \hat c_v\f$.
     */
    virtual doublereal cp_mole() const;

    /**
     * Molar heat capacity at constant volume. Units: J/kmol/K.
     * For an incompressible substance, \f$ \hat c_p = \hat c_v\f$.
     */
    virtual doublereal cv_mole() const;

    //! @}
    //! @name Mechanical Equation of State
    //! @{

    //! Report the Pressure. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent
     * of pressure. This method simply returns the stored
     * pressure value.
     */
    virtual doublereal pressure() const;

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     * For an incompressible substance, the density is
     * independent of pressure. Therefore, this method only
     * stores the specified pressure value. It does not
     * modify the density.
     *
     * @param p Pressure (units - Pa)
     */
    virtual void setPressure(doublereal p);

    //! @}
    //! @name Chemical Potentials and Activities
    //! @{

    /**
     * This method returns the array of generalized
     * concentrations.  For a stoichiometric substance, there is
     * only one species, and the generalized concentration is 1.0.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    /**
     * The standard concentration. This is defined as the concentration
     * by which the generalized concentration is normalized to produce
     * the activity.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    /**
     * Returns the natural logarithm of the standard
     * concentration of the kth species
     */
    virtual doublereal logStandardConc(size_t k=0) const;

    /**
     * Get the array of chemical potentials at unit activity
     * \f$ \mu^0_k \f$.
     *
     * For a stoichiometric substance, there is no activity term in
     * the chemical potential expression, and therefore the
     * standard chemical potential and the chemical potential
     * are both equal to the molar Gibbs function.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    /**
     * Returns the units of the standard and generalized
     * concentrations. Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     *     uA[0] = kmol units - default  = 0
     *     uA[1] = m    units - default  = 0
     *     uA[2] = kg   units - default  = 0;
     *     uA[3] = Pa(pressure) units - default = 0;
     *     uA[4] = Temperature units - default = 0;
     *     uA[5] = time units - default = 0
     * @deprecated To be removed after Cantera 2.2.
     */
    virtual void getUnitsStandardConc(double* uA, int k = 0,
                                      int sizeUA = 6) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    /**
     * Get the array of non-dimensional chemical potentials
     * \f$ \mu_k / \hat R T \f$.
     */
    virtual void getChemPotentials_RT(doublereal* mu) const;

    /**
     * For a stoichiometric substance, there is only one species.
     * This method returns the molar Gibbs function in the
     * first element of array \c mu.
     */
    virtual void getChemPotentials(doublereal* mu) const;

    /**
     * Get the species electrochemical potentials. Units: J/kmol.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the
     * to each chemical potential.
     */
    void getElectrochemPotentials(doublereal* mu) const;

    /**
     * Returns an array of partial molar enthalpies for the species
     * in the mixture.
     * Units (J/kmol)
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    /**
     * Returns an array of partial molar entropies of the species in the
     * solution. Units: J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Get the partial molar heat capacities Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat capacities
     *                at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    /**
     * returns an array of partial molar volumes of the species
     * in the solution. Units: m^3 kmol-1.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{
    /**
     * Get the nondimensional Enthalpy functions for the species
     * at their standard states at the current
     * <I>T</I> and <I>P</I> of the solution.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    /**
     * Get the array of nondimensional Enthalpy functions for the
     * standard state species
     * at the current <I>T</I> and <I>P</I> of the solution.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    /**
     * Get the nondimensional Gibbs functions for the species
     * at their standard states of solution at the current T and P
     * of the solution.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the Gibbs functions for the standard
    //! state of the species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * Units are Joules/kmol
     * @param gpure  Output vector of  standard state Gibbs free energies
     *               Length: m_kk.
     */
    virtual void getPureGibbs(doublereal* gpure) const;

    /**
     * Get the nondimensional Heat Capacities at constant
     * pressure for the standard state of the species
     * at the current T and P.
     */
    virtual void getCp_R(doublereal* cpr) const;

    /**
     * Get the standard volumes for the standard state of the species
     * at the current T and P
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
     *  This function fills in its one entry in hrt[] by calling
     *  the underlying species thermo function for the
     *  dimensionless enthalpy.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /**
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     *  This function fills in its one entry in hrt[] by calling
     *  the underlying species thermo function for the
     *  dimensionless Gibbs free energy, calculated from the
     *  dimensionless enthalpy and entropy.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    /**
     *  Returns the vector of the
     *  Gibbs function of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *  units = J/kmol
     *
     *  This function fills in its one entry in g[] by calling
     *  the underlying species thermo functions for the
     *  Gibbs free energy, calculated from enthalpy and the
     *  entropy, and the multiplying by RT.
     */
    virtual void  getGibbs_ref(doublereal* g) const;

    /**
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     *  This function fills in its one entry in hrt[] by calling
     *  the underlying species thermo function for the
     *  dimensionless entropy.
     */
    virtual void getEntropy_R_ref(doublereal* er) const;

    //!  Returns the vector of nondimensional
    //!  constant pressure heat capacities of the reference state
    //!  at the current temperature of the solution
    //!  and reference pressure for each species.
    /*!
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const;
    //! @}

    virtual void initThermo();

    virtual void setParameters(int n, double* const c);

    virtual void getParameters(int& n, double* const c) const;

    virtual void setParametersFromXML(const XML_Node& eosdata);

protected:
    doublereal m_press;
    doublereal m_p0;

    mutable vector_fp      m_h0_RT;
    mutable vector_fp      m_cp0_R;
    mutable vector_fp      m_s0_R;

private:
    void _updateThermo() const;
};

}

#endif
