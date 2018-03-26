/**
 *  @file MineralEQ3.h
 * Header file for the MineralEQ3 class, which represents a fixed-composition
 * incompressible substance based on EQ3's parameterization (see \ref thermoprops and
 * class \link Cantera::MineralEQ3 MineralEQ3\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_MINERALEQ3_H
#define CT_MINERALEQ3_H

#include "StoichSubstance.h"

namespace Cantera
{

//! Class MineralEQ3 represents a stoichiometric (fixed composition)
//! incompressible substance based on EQ3's parameterization
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * @deprecated To be removed after Cantera 2.4
 *
 * This class inherits from SingleSpeciesTP class. EQ's parameterization is
 * mapped onto the Shomate polynomial class.
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
 * giving the standard-state enthalpy, the term \f$ P_0 \hat v\f$ is subtracted
 * from the specified molar enthalpy to compute the molar internal energy. The
 * entropy is assumed to be independent of the pressure.
 *
 * The enthalpy function is given by the following relation.
 *
 *       \f[
 *              h^o_k(T,P) =
 *                  h^{ref}_k(T) + \tilde v \left( P - P_{ref} \right)
 *       \f]
 *
 * For an incompressible, stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties are specified by
 * giving the standard-state enthalpy, the term \f$ P_{ref} \tilde v\f$ is
 * subtracted from the specified reference molar enthalpy to compute the molar
 * internal energy.
 *
 *       \f[
 *            u^o_k(T,P) = h^{ref}_k(T) - P_{ref} \tilde v
 *       \f]
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
 * ## %Application within Kinetics Managers
 *
 * The standard concentration is equal to 1.0. This means that the kinetics
 * operator works on an (activities basis). Since this is a stoichiometric
 * substance, this means that the concentration of this phase drops out of
 * kinetics expressions.
 *
 * An example of a reaction using this is a sticking coefficient reaction of a
 * substance in an ideal gas phase on a surface with a bulk phase species in
 * this phase. In this case, the rate of progress for this reaction,
 * \f$ R_s \f$, may be expressed via the following equation:
 *   \f[
 *    R_s = k_s C_{gas}
 *   \f]
 * where the units for \f$ R_s \f$ are kmol m-2 s-1. \f$ C_{gas} \f$ has units
 * of kmol m-3. Therefore, the kinetic rate constant, \f$ k_s \f$, has units of
 * m s-1. Nowhere does the concentration of the bulk phase appear in the rate
 * constant expression, since it's a stoichiometric phase and the activity is
 * always equal to 1.0.
 *
 * @ingroup thermoprops
 */
class MineralEQ3 : public StoichSubstance
{
public:
    //! Default constructor for the MineralEQ3 class
    MineralEQ3() {
        warn_deprecated("Class MineralEQ3", "To be removed after Cantera 2.4");
    }

    //! Construct and initialize a MineralEQ3 ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    MineralEQ3(const std::string& infile, const std::string& id = "");

    //! Construct and initialize a MineralEQ3 ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML node pointing to a MineralEQ3 description
     *  @param id       Id of the phase.
     */
    MineralEQ3(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "MineralEQ3";
    }

    //! @name Mechanical Equation of State
    //! @{

    //! Report the Pressure. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent of pressure.
     * This method simply returns the stored pressure value.
     */
    virtual doublereal pressure() const;

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent of pressure.
     * Therefore, this method only stores the specified pressure value. It does
     * not modify the density.
     *
     * @param p Pressure (units - Pa)
     */
    virtual void setPressure(doublereal p);

    virtual doublereal isothermalCompressibility() const;
    virtual doublereal thermalExpansionCoeff() const;

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * This section is largely handled by parent classes, since there is only
     * one species. Therefore, the activity is equal to one.
     * @{
     */

    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k / C^0_k, \f$ where
     * \f$ C^0_k \f$ is a standard concentration defined below and \f$ a_k \f$
     * are activities used in the thermodynamic functions.  These activity (or
     * generalized) concentrations are used by kinetics manager classes to
     * compute the forward and reverse rates of elementary reactions.
     *
     * For a stoichiometric substance, there is only one species, and the
     * generalized concentration is 1.0.
     *
     * @param c Output array of generalized concentrations. The units depend
     *           upon the implementation of the reaction rate expressions within
     *           the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize the activity
     * (i.e., generalized) concentration. This phase assumes that the kinetics
     * operator works on an dimensionless basis. Thus, the standard
     * concentration is equal to 1.0.
     *
     * @param k Optional parameter indicating the species. The default
     *         is to assume this refers to species 0.
     * @return
     *   Returns The standard Concentration as 1.0
     */
    virtual doublereal standardConcentration(size_t k=0) const;
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Get the array of chemical potentials at unit activity for the species at
    //! their standard states at the current *T* and *P* of the solution.
    /*!
     * For a stoichiometric substance, there is no activity term in the chemical
     * potential expression, and therefore the standard chemical potential and
     * the chemical potential are both equal to the molar Gibbs function.
     *
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu0     Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

    virtual void getEnthalpy_RT(doublereal* hrt) const;
    virtual void getEntropy_R(doublereal* sr) const;
    virtual void getGibbs_RT(doublereal* grt) const;
    virtual void getCp_R(doublereal* cpr) const;

    //!  Returns the vector of nondimensional Internal Energies of the standard
    //!  state species at the current *T* and *P* of the solution
    /*!
     * For an incompressible, stoichiometric substance, the molar internal
     * energy is independent of pressure. Since the thermodynamic properties are
     * specified by giving the standard-state enthalpy, the term
     * \f$ P_{ref} \hat v\f$ is subtracted from the specified reference molar
     * enthalpy to compute the standard state molar internal energy.
     *
     * @param urt  output vector of nondimensional standard state internal
     *             energies of the species. Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    virtual void getIntEnergy_RT_ref(doublereal* urt) const;
    //! @}

    //! @copydoc ThermoPhase::initThermoXML
    /*!
     * This is the main routine for reading in activity coefficient parameters.
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Set the equation of state parameters
    /*!
     * @internal
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     *        c[0] = density of phase [ kg/m3 ]
     */
    virtual void setParameters(int n, doublereal* const c);

    //! Get the equation of state parameters in a vector
    /*!
     * @internal
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     *
     *  For this phase:
     *       -  n = 1
     *       -  c[0] = density of phase [ kg/m3 ]
     */
    virtual void getParameters(int& n, doublereal* const c) const;

    //! @copydoc ThermoPhase::setParametersFromXML
    /*!
     * For this phase, the density of the phase is specified in this block.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);
    doublereal LookupGe(const std::string& elemName);
    void convertDGFormation();

protected:
    //! Value of the Absolute Gibbs Free Energy NIST scale at T_r and P_r
    /*!
     * This is the NIST scale value of Gibbs free energy at T_r = 298.15
     * and P_r = 1 atm.
     *
     *  J kmol-1
     */
    doublereal m_Mu0_pr_tr;

    //! Input value of S_j at Tr and Pr    (cal gmol-1 K-1)
    /*!
     *  Tr = 298.15   Pr = 1 atm
     */
    doublereal m_Entrop_pr_tr;

    //! Input Value of deltaG of Formation at Tr and Pr    (cal gmol-1)
    /*!
     * Tr = 298.15   Pr = 1 atm
     *
     * This is the delta G for the formation reaction of the ion from elements
     * in their stable state at Tr, Pr.
     */
    doublereal m_deltaG_formation_pr_tr;

    //! Input Value of deltaH of Formation at Tr and Pr    (cal gmol-1)
    /*!
     * Tr = 298.15   Pr = 1 atm
     *
     * This is the delta H for the formation reaction of the ion from elements
     * in their stable state at Tr, Pr.
     */
    doublereal m_deltaH_formation_pr_tr;

    //! Input Value of the molar volume at T_r and P_r
    /*!
     *  cm^3 / gmol
     */
    doublereal m_V0_pr_tr;

    //! a coefficient (cal gmol-1 K-1)
    doublereal m_a;

    //! b coefficient (cal gmol-1 K-2) x 10^3
    doublereal m_b;

    //! c coefficient (cal K gmol-1 K) x 10^-5
    doublereal m_c;
};

}

#endif
