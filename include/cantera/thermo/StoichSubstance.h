/**
 * @file StoichSubstance.h
 * Header file for the StoichSubstance class, which represents a fixed-composition
 * incompressible substance (see @ref thermoprops and
 * class @link Cantera::StoichSubstance StoichSubstance@endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STOICHSUBSTANCE_H
#define CT_STOICHSUBSTANCE_H

#include "SingleSpeciesTP.h"

namespace Cantera
{

//! Class StoichSubstance represents a stoichiometric (fixed composition)
//! incompressible substance.
/*!
 * This class internally changes the independent degree of freedom from density
 * to pressure. This is necessary because the phase is incompressible. It uses a
 * constant volume approximation.
 *
 * ## Specification of Species Standard State Properties
 *
 * This class inherits from SingleSpeciesTP. It is assumed that the reference
 * state thermodynamics may be obtained by a pointer to a populated species
 * thermodynamic property manager class (see ThermoPhase::m_spthermo). How to
 * relate pressure changes to the reference state thermodynamics is resolved at
 * this level.
 *
 * For an incompressible, stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties are specified by
 * giving the standard-state enthalpy, the term @f$ P_0 \hat v @f$ is subtracted
 * from the specified molar enthalpy to compute the molar internal energy. The
 * entropy is assumed to be independent of the pressure.
 *
 * The enthalpy function is given by the following relation.
 *
 * @f[
 *              h^o_k(T,P) =
 *                  h^{ref}_k(T) + \tilde v \left( P - P_{ref} \right)
 * @f]
 *
 * For an incompressible, stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties are specified by
 * giving the standard-state enthalpy, the term @f$ P_{ref} \tilde v @f$ is
 * subtracted from the specified reference molar enthalpy to compute the molar
 * internal energy.
 *
 * @f[
 *            u^o_k(T,P) = h^{ref}_k(T) - P_{ref} \tilde v
 * @f]
 *
 * The standard state heat capacity and entropy are independent of pressure. The
 * standard state Gibbs free energy is obtained from the enthalpy and entropy
 * functions.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * All solution properties are obtained from the standard state species
 * functions, since there is only one species in the phase.
 *
 * ## Application within Kinetics Managers
 *
 * The standard concentration is equal to 1.0. This means that the kinetics
 * operator works on an (activities basis). Since this is a stoichiometric
 * substance, this means that the concentration of this phase drops out of
 * kinetics expressions.
 *
 * An example of a reaction using this is a sticking coefficient reaction of a
 * substance in an ideal gas phase on a surface with a bulk phase species in
 * this phase. In this case, the rate of progress for this reaction,
 * @f$ R_s @f$, may be expressed via the following equation:
 *   @f[
 *    R_s = k_s C_{gas}
 *   @f]
 * where the units for @f$ R_s @f$ are kmol m-2 s-1. @f$ C_{gas} @f$ has units
 * of kmol m-3. Therefore, the kinetic rate constant, @f$ k_s @f$, has units of
 * m s-1. Nowhere does the concentration of the bulk phase appear in the rate
 * constant expression, since it's a stoichiometric phase and the activity is
 * always equal to 1.0.
 *
 * @ingroup thermoprops
 */
class StoichSubstance : public SingleSpeciesTP
{
public:
    //! Construct and initialize a StoichSubstance ThermoPhase object directly
    //! from an input file
    /*!
     * @param infile name of the input file. If blank, an empty phase will be
     *               created.
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    explicit StoichSubstance(const string& infile="", const string& id="");

    string type() const override {
        return "fixed-stoichiometry";
    }

    bool isCompressible() const override {
        return false;
    }

    //! @name Mechanical Equation of State
    //! @{

    //! Report the Pressure. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent of pressure.
     * This method simply returns the stored pressure value.
     */
    double pressure() const override;

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent of pressure.
     * Therefore, this method only stores the specified pressure value. It does
     * not modify the density.
     *
     * @param p Pressure (units - Pa)
     */
    void setPressure(double p) override;

    double isothermalCompressibility() const override;
    double thermalExpansionCoeff() const override;

    //! @}
    //! @name Activities, Standard States, and Activity Concentrations
    //!
    //! This section is largely handled by parent classes, since there
    //! is only one species. Therefore, the activity is equal to one.
    //! @{

    Units standardConcentrationUnits() const override;

    //! This method returns an array of generalized concentrations
    /*!
     * @f$ C^a_k @f$ are defined such that @f$ a_k = C^a_k / C^0_k, @f$ where
     * @f$ C^0_k @f$ is a standard concentration defined below and @f$ a_k @f$
     * are activities used in the thermodynamic functions. These activity (or
     * generalized) concentrations are used by kinetics manager classes to
     * compute the forward and reverse rates of elementary reactions.
     *
     * For a stoichiometric substance, there is only one species, and the
     * generalized concentration is 1.0.
     *
     * @param c Output array of generalized concentrations. The
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    void getActivityConcentrations(double* c) const override;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration @f$ C^0_k @f$ used to normalize the activity
     * (that is, generalized) concentration. This phase assumes that the kinetics
     * operator works on an dimensionless basis. Thus, the standard
     * concentration is equal to 1.0.
     *
     * @param k Optional parameter indicating the species. The default
     *         is to assume this refers to species 0.
     * @return
     *   Returns The standard Concentration as 1.0
     */
    double standardConcentration(size_t k=0) const override;
    double logStandardConc(size_t k=0) const override;

    //! Get the array of chemical potentials at unit activity for the species at
    //! their standard states at the current *T* and *P* of the solution.
    /*!
     * For a stoichiometric substance, there is no activity term in the chemical
     * potential expression, and therefore the standard chemical potential and
     * the chemical potential are both equal to the molar Gibbs function.
     *
     * These are the standard state chemical potentials @f$ \mu^0_k(T,P) @f$.
     * The values are evaluated at the current temperature and pressure of the
     * solution
     *
     * @param mu0     Output vector of chemical potentials.
     *                Length: m_kk.
     */
    void getStandardChemPotentials(double* mu0) const override;

    //! @}
    //! @name  Properties of the Standard State of the Species in the Solution
    //! @{

    void getEnthalpy_RT(double* hrt) const override;
    void getEntropy_R(double* sr) const override;
    void getGibbs_RT(double* grt) const override;
    void getCp_R(double* cpr) const override;

    //! Returns the vector of nondimensional Internal Energies of the standard
    //! state species at the current *T* and *P* of the solution
    /*!
     * For an incompressible, stoichiometric substance, the molar internal
     * energy is independent of pressure. Since the thermodynamic properties
     * are specified by giving the standard-state enthalpy, the term
     * @f$ P_{ref} \hat v @f$ is subtracted from the specified reference molar
     * enthalpy to compute the standard state molar internal energy.
     *
     * @param urt  output vector of nondimensional standard state
     *             internal energies of the species. Length: m_kk.
     */
    void getIntEnergy_RT(double* urt) const override;

    //! @}
    //! @name Thermodynamic Values for the Species Reference States
    //! @{

    //! Returns the vector of nondimensional internal Energies of the reference
    //! state at the current temperature of the solution and the reference
    //! pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state internal
     *               energies of the species. Length: m_kk
     */
    void getIntEnergy_RT_ref(double* urt) const override;
    //! @}

    void initThermo() override;
    void getSpeciesParameters(const string& name, AnyMap& speciesNode) const override;
};

}

#endif
