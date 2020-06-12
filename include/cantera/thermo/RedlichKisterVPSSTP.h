/**
 *  @file  RedlichKisterVPSSTP.h (see \ref thermoprops and class \link
 *      Cantera::RedlichKisterVPSSTP RedlichKisterVPSSTP\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REDLICHKISTERVPSSTP_H
#define CT_REDLICHKISTERVPSSTP_H

#include "cantera/thermo/GibbsExcessVPSSTP.h"

namespace Cantera
{

//! RedlichKisterVPSSTP is a derived class of GibbsExcessVPSSTP that employs the
//! Redlich-Kister approximation for the excess Gibbs free energy
/*!
 * RedlichKisterVPSSTP derives from class GibbsExcessVPSSTP which is derived
 * from VPStandardStateTP, and overloads the virtual methods defined there with
 * ones that use expressions appropriate for the Redlich Kister Excess Gibbs
 * free energy approximation.
 *
 * The independent unknowns are pressure, temperature, and mass fraction.
 *
 * ## Specification of Species Standard State Properties
 *
 * All species are defined to have standard states that depend upon both the
 * temperature and the pressure. The Redlich-Kister approximation assumes
 * symmetric standard states, where all of the standard state assume that the
 * species are in pure component states at the temperature and pressure of the
 * solution.  I don't think it prevents, however, some species from being dilute
 * in the solution.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The molar excess Gibbs free energy is given by the following formula which is
 * a sum over interactions *i*. Each of the interactions are binary interactions
 * involving two of the species in the phase, denoted, *Ai* and *Bi*. This is
 * the generalization of the Redlich-Kister formulation for a phase that has
 * more than 2 species.
 *
 * \f[
 *     G^E = \sum_{i} G^E_{i}
 * \f]
 *
 * where
 *
 * \f[
 *    G^E_{i} =   n X_{Ai} X_{Bi} \sum_m \left( A^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 * \f]
 *
 * and where we can break down the Gibbs free energy contributions into enthalpy and entropy contributions
 *
 * \f[
 *    H^E_i = n X_{Ai} X_{Bi} \sum_m \left( H^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 * \f]
 *
 * \f[
 *    S^E_i = n X_{Ai} X_{Bi} \sum_m \left( S^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 * \f]
 *
 * where n is the total moles in the solution. The activity of a species defined
 * in the phase is given by an excess Gibbs free energy formulation.
 *
 * \f[
 *      a_k = \gamma_k  X_k
 * \f]
 *
 * where
 *
 * \f[
 *      R T \ln( \gamma_k )= \frac{d(n G^E)}{d(n_k)}\Bigg|_{n_i}
 * \f]
 *
 * Taking the derivatives results in the following expression
 * \f[
 *      R T \ln( \gamma_k )=  \sum_i \delta_{Ai,k} (1 - X_{Ai}) X_{Bi} \sum_m \left( A^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 *                           + \sum_i \delta_{Ai,k} X_{Ai} X_{Bi} \sum_m \left(  A^{i}_0 +  A^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^{m-1} (1 - X_{Ai} + X_{Bi}) \right)
 * \f]
 *
 * This object inherits from the class VPStandardStateTP. Therefore, the
 * specification and calculation of all standard state and reference state
 * values are handled at that level. Various functional forms for the standard
 * state are permissible. The chemical potential for species *k* is equal to
 *
 * \f[
 *      \mu_k(T,P) = \mu^o_k(T, P) + R T \ln(\gamma_k X_k)
 * \f]
 *
 * The partial molar entropy for species *k* is given by the following relation,
 *
 * \f[
 *       \tilde{s}_k(T,P) =  s^o_k(T,P)  - R \ln( \gamma_k X_k )
 *              - R T \frac{d \ln(\gamma_k) }{dT}
 * \f]
 *
 * The partial molar enthalpy for species *k* is given by
 *
 * \f[
 *      \tilde{h}_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
 * \f]
 *
 * The partial molar volume for species *k* is
 *
 * \f[
 *        \tilde V_k(T,P)  = V^o_k(T,P)  + R T \frac{d \ln(\gamma_k) }{dP}
 * \f]
 *
 * The partial molar Heat Capacity for species *k* is
 *
 * \f[
 *      \tilde{C}_{p,k}(T,P) = C^o_{p,k}(T,P)   - 2 R T \frac{d \ln( \gamma_k )}{dT}
 *              - R T^2 \frac{d^2 \ln(\gamma_k) }{{dT}^2}
 * \f]
 *
 * ## %Application within Kinetics Managers
 *
 * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k / C^s_k, \f$ where
 * \f$ C^s_k \f$ is a standard concentration defined below and \f$ a_k \f$ are
 * activities used in the thermodynamic functions.  These activity (or
 * generalized) concentrations are used by kinetics manager classes to compute
 * the forward and reverse rates of elementary reactions. The activity
 * concentration,\f$  C^a_k \f$,is given by the following expression.
 *
 * \f[
 *      C^a_k = C^s_k  X_k  = \frac{P}{R T} X_k
 * \f]
 *
 * The standard concentration for species *k* is independent of *k* and equal to
 *
 * \f[
 *     C^s_k =  C^s = \frac{P}{R T}
 * \f]
 *
 * For example, a bulk-phase binary gas reaction between species j and k,
 * producing a new gas species l would have the following equation for its rate
 * of progress variable, \f$ R^1 \f$, which has units of kmol m-3 s-1.
 *
 * \f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1 (C^s a_j) (C^s a_k)
 * \f]
 * where
 * \f[
 *    C_j^a = C^s a_j \mbox{\quad and \quad} C_k^a = C^s a_k
 * \f]
 *
 * \f$ C_j^a \f$ is the activity concentration of species j, and \f$ C_k^a \f$
 * is the activity concentration of species k. \f$ C^s \f$ is the standard
 * concentration. \f$ a_j \f$ is the activity of species j which is equal to the
 * mole fraction of j.
 *
 * The reverse rate constant can then be obtained from the law of microscopic
 * reversibility and the equilibrium expression for the system.
 *
 * \f[
 *       \frac{a_j a_k}{ a_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 * \f]
 *
 * \f$ K_a^{o,1} \f$ is the dimensionless form of the equilibrium constant,
 * associated with the pressure dependent standard states \f$ \mu^o_l(T,P) \f$
 * and their associated activities, \f$ a_l \f$, repeated here:
 *
 * \f[
 *      \mu_l(T,P) = \mu^o_l(T, P) + R T \log(a_l)
 * \f]
 *
 * We can switch over to expressing the equilibrium constant in terms of the
 * reference state chemical potentials
 *
 * \f[
 *     K_a^{o,1} = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{P}
 * \f]
 *
 * The concentration equilibrium constant, \f$ K_c \f$, may be obtained by
 * changing over to activity concentrations. When this is done:
 *
 * \f[
 *     \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
 *         \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{RT}
 * \f]
 *
 * %Kinetics managers will calculate the concentration equilibrium constant, \f$
 * K_c \f$, using the second and third part of the above expression as a
 * definition for the concentration equilibrium constant.
 *
 * For completeness, the pressure equilibrium constant may be obtained as well
 *
 * \f[
 *     \frac{P_j P_k}{ P_l P_{ref}} = K_p^1 = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} )
 * \f]
 *
 * \f$ K_p \f$ is the simplest form of the equilibrium constant for ideal gases.
 * However, it isn't necessarily the simplest form of the equilibrium constant
 * for other types of phases; \f$ K_c \f$ is used instead because it is
 * completely general.
 *
 * The reverse rate of progress may be written down as
 * \f[
 *    R^{-1} = k^{-1} C_l^a =  k^{-1} (C^o a_l)
 * \f]
 *
 * where we can use the concept of microscopic reversibility to write the
 * reverse rate constant in terms of the forward rate constant and the
 * concentration equilibrium constant, \f$ K_c \f$.
 *
 * \f[
 *     k^{-1} =  k^1 K^1_c
 * \f]
 *
 * \f$k^{-1} \f$ has units of s-1.
 *
 * @ingroup thermoprops
 */
class RedlichKisterVPSSTP : public GibbsExcessVPSSTP
{
public:
    //! Constructor
    /*!
     * This doesn't do much more than initialize constants with default values.
     */
    RedlichKisterVPSSTP();

    //! Construct a RedlichKisterVPSSTP object from an input file
    /*!
     * @param inputFile Name of the input file containing the phase definition
     * @param id        name (ID) of the phase in the input file. If empty, the
     *                  first phase definition in the input file will be used.
     */
    RedlichKisterVPSSTP(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize a RedlichKisterVPSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    RedlichKisterVPSSTP(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "RedlichKister";
    }

    //! @name  Molar Thermodynamic Properties
    //! @{

    virtual doublereal enthalpy_mole() const;
    virtual doublereal entropy_mole() const;
    virtual doublereal cp_mole() const;
    virtual doublereal cv_mole() const;

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and pressure.
     * @{
     */

    virtual void getLnActivityCoefficients(doublereal* lnac) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    virtual void getChemPotentials(doublereal* mu) const;

    //! Returns an array of partial molar enthalpies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the molality-based
     * activity coefficient wrt temperature
     *
     *  \f[
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *  \f]
     *
     * @param hbar  Vector of returned partial molar enthalpies
     *              (length m_kk, units = J/kmol)
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the activity coefficient
     * wrt temperature
     *
     *  \f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *  \f]
     *
     * @param sbar  Vector of returned partial molar entropies
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Returns an array of partial molar entropies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the activity coefficient
     * wrt temperature
     *
     *  \f[
     *   ???????????????
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *   ???????????????
     *  \f]
     *
     * @param cpbar  Vector of returned partial molar heat capacities
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Get the array of temperature second derivatives of the log activity
    //! coefficients
    /*!
     *  units = 1/Kelvin
     *
     * @param d2lnActCoeffdT2  Output vector of temperature 2nd derivatives of
     *                         the log Activity Coefficients. length = m_kk
     */
    virtual void getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const;

    virtual void getdlnActCoeffdT(doublereal* dlnActCoeffdT) const;

    /// @}
    /// @name Initialization
    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see importPhase().

    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Add a binary species interaction with the specified parameters
    /*!
     * @param speciesA         name of the first species
     * @param speciesB         name of the second species
     * @param excess_enthalpy  coefficients of the excess enthalpy polynomial
     * @param n_enthalpy       number of excess enthalpy polynomial coefficients
     * @param excess_entropy   coefficients of the excess entropy polynomial
     * @param n_entropy        number of excess entropy polynomial coefficients
     */
    void addBinaryInteraction(const std::string& speciesA, const std::string& speciesB,
        const double* excess_enthalpy, size_t n_enthalpy,
        const double* excess_entropy, size_t n_entropy);

    //! @}
    //! @name  Derivatives of Thermodynamic Variables needed for Applications
    //! @{

    virtual void getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds, doublereal* dlnActCoeffds) const;
    virtual void getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const;
    virtual void getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const;
    virtual void getdlnActCoeffdlnN(const size_t ld, doublereal* const dlnActCoeffdlnN);

    //@}

private:
    //! Process an XML node called "binaryNeutralSpeciesParameters"
    /*!
     * This node contains all of the parameters necessary to describe the
     * Redlich-Kister model for a particular binary interaction. This function
     * reads the XML file and writes the coefficients it finds to an internal
     * data structures.
     *
     * @param xmlBinarySpecies  Reference to the XML_Node named
     *     "binaryNeutralSpeciesParameters" containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node& xmlBinarySpecies);

    //! Initialize lengths of local variables after all species have been
    //! identified.
    void initLengths();

    //! Update the activity coefficients
    /*!
     * This function will be called to update the internally stored natural
     * logarithm of the activity coefficients
     */
    void s_update_lnActCoeff() const;

    //! Update the derivative of the log of the activity coefficients wrt T
    /*!
     * This function will be called to update the internally stored derivative
     * of the natural logarithm of the activity coefficients wrt temperature.
     */
    void s_update_dlnActCoeff_dT() const;

    //! Internal routine that calculates the derivative of the activity
    //! coefficients wrt the mole fractions.
    /*!
     * This routine calculates the the derivative of the activity coefficients
     * wrt to mole fraction with all other mole fractions held constant. This is
     * strictly not permitted. However, if the resulting matrix is multiplied by
     * a permissible deltaX vector then everything is ok.
     *
     * This is the natural way to handle concentration derivatives in this
     * routine.
     */
    void s_update_dlnActCoeff_dX_() const;

    //! Internal routine that calculates the total derivative of the activity
    //! coefficients with respect to the log of the mole fractions.
    /*!
     * This function will be called to update the internally stored vector of
     * the total derivatives (i.e. not assuming other mole fractions are
     * constant) of the natural logarithm of the activity coefficients with
     * respect to the log of the mole fraction.
     */
    void s_update_dlnActCoeff_dlnX_diag() const;

protected:
    //! number of binary interaction expressions
    size_t numBinaryInteractions_;

    //! vector of species indices representing species A in the interaction
    /*!
     *  Each Redlich-Kister excess Gibbs free energy term involves two species,
     *  A and B. This vector identifies species A.
     */
    std::vector<size_t> m_pSpecies_A_ij;

    //! vector of species indices representing species B in the interaction
    /*!
     *  Each Redlich-Kister excess Gibbs free energy term involves two species,
     *  A and B. This vector identifies species B.
     */
    std::vector<size_t> m_pSpecies_B_ij;

    //! Vector of the length of the polynomial for the interaction.
    std::vector<size_t> m_N_ij;

    //! Enthalpy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable std::vector< vector_fp> m_HE_m_ij;

    //! Entropy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable std::vector< vector_fp> m_SE_m_ij;

    //! form of the RedlichKister interaction expression. Currently there is
    //! only one form.
    int formRedlichKister_;

    //! form of the temperature dependence of the Redlich-Kister interaction
    //! expression. Currently there is only one form -> constant wrt
    //! temperature.
    int formTempModel_;

    //! Two dimensional array of derivatives of activity coefficients wrt mole
    //! fractions
    mutable Array2D dlnActCoeff_dX_;
};

}

#endif
