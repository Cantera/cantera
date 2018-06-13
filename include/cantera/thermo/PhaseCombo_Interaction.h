/**
 *  @file  PhaseCombo_Interaction.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ the Margules Gibbs free energy formulation and eliminates the ideal mixing term.
 *  (see \ref thermoprops
 *    and class \link Cantera::PhaseCombo_Interaction PhaseCombo_Interaction\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_PHASECOMBO_INTERACTION_H
#define CT_PHASECOMBO_INTERACTION_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera
{

//! PhaseCombo_Interaction is a derived class of GibbsExcessVPSSTP that employs
//! the Margules approximation for the excess Gibbs free energy while
//! eliminating the entropy of mixing term.
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
 * PhaseCombo_Interaction derives from class GibbsExcessVPSSTP which is derived
 * from VPStandardStateTP, and overloads the virtual methods defined there with
 * ones that use expressions appropriate for the Margules Excess Gibbs free
 * energy approximation. The reader should refer to the MargulesVPSSTP class for
 * information on that class. This class in addition adds a term to the activity
 * coefficient that eliminates the ideal solution mixing term within the
 * chemical potential. This is a very radical thing to do, but it is supported
 * by experimental evidence under some conditions.
 *
 * The independent unknowns are pressure, temperature, and mass fraction.
 *
 * This class is introduced to represent specific conditions observed in thermal
 * batteries. HOwever, it may be physically motivated to represent conditions
 * where there may be a mixture of compounds that are not "mixed" at the
 * molecular level. Therefore, there is no mixing term.
 *
 * The lack of a mixing term has profound effects. First, the mole fraction of a
 * species can now be identically zero due to thermodynamic considerations. The
 * phase behaves more like a series of phases. That's why we named it
 * PhaseCombo.
 *
 * ## Specification of Species Standard State Properties
 *
 * All species are defined to have standard states that depend upon both the
 * temperature and the pressure. The Margules approximation assumes symmetric
 * standard states, where all of the standard state assume that the species are
 * in pure component states at the temperature and pressure of the solution.  I
 * don't think it prevents, however, some species from being dilute in the
 * solution.
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * The molar excess Gibbs free energy is given by the following formula which is
 * a sum over interactions *i*. Each of the interactions are binary interactions
 * involving two of the species in the phase, denoted, *Ai* and *Bi*. This is
 * the generalization of the Margules formulation for a phase that has more than
 * 2 species. The second term in the excess Gibbs free energy is a negation of
 * the ideal solution's mixing term.
 *
 * \f[
 *     G^E = \sum_i \left(  H_{Ei} - T S_{Ei} \right) -  \sum_i \left( n_i  R T \ln{X_i} \right)
 * \f]
 * \f[
 *    H^E_i = n X_{Ai} X_{Bi} \left( h_{o,i} +  h_{1,i} X_{Bi} \right)
 * \f]
 * \f[
 *    S^E_i = n X_{Ai} X_{Bi} \left( s_{o,i} +  s_{1,i} X_{Bi} \right)
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
 *
 * \f[
 *      R T \ln( \gamma_k )= \sum_i \left( \left( \delta_{Ai,k} X_{Bi} + \delta_{Bi,k} X_{Ai}  - X_{Ai} X_{Bi} \right)
 *       \left( g^E_{o,i} +  g^E_{1,i} X_{Bi} \right) +
 *       \left( \delta_{Bi,k} - X_{Bi} \right)      X_{Ai} X_{Bi}  g^E_{1,i} \right) - RT \ln{X_k}
 * \f]
 *
 * where \f$  g^E_{o,i} =  h_{o,i} - T s_{o,i} \f$ and
 * \f$ g^E_{1,i} =  h_{1,i} - T s_{1,i} \f$ and where \f$ X_k \f$ is the mole
 * fraction of species *k*.
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
 * The partial molar volume for  species *k* is
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
 *
 * where
 *
 * \f[
 *      C_j^a = C^s a_j \mbox{\quad and \quad} C_k^a = C^s a_k
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
 *     \frac{a_j a_k}{ a_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 * \f]
 *
 * \f$ K_a^{o,1} \f$ is the dimensionless form of the equilibrium constant,
 * associated with the pressure dependent standard states \f$ \mu^o_l(T,P) \f$
 * and their associated activities, \f$ a_l \f$, repeated here:
 *
 * \f[
 *     \mu_l(T,P) = \mu^o_l(T, P) + R T \log(a_l)
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
 * Kinetics managers will calculate the concentration equilibrium constant,
 * \f$  K_c \f$, using the second and third part of the above expression as a
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
 *     R^{-1} = k^{-1} C_l^a =  k^{-1} (C^o a_l)
 * \f]
 *
 * where we can use the concept of microscopic reversibility to write the
 * reverse rate constant in terms of the forward reate constant and the
 * concentration equilibrium constant, \f$ K_c \f$.
 *
 * \f[
 *    k^{-1} =  k^1 K^1_c
 * \f]
 *
 * \f$k^{-1} \f$ has units of s-1.
 *
 * ## Instantiation of the Class
 *
 * The constructor for this phase is located in the default ThermoFactory for
 * %Cantera. A new PhaseCombo_Interaction object may be created by the following
 * code snippet:
 *
 * @code
 *    XML_Node *xc = get_XML_File("LiFeS_X_combo.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "LiFeS_X");
 *    ThermoPhase *l_tp = newPhase(*xs);
 *    PhaseCombo_Interaction *LiFeS_X_solid = dynamic_cast <PhaseCombo_Interaction *>(l_tp);
 * @endcode
 *
 * or by the following code
 *
 * @code
 *    std::string id = "LiFeS_X";
 *    Cantera::ThermoPhase *LiFeS_X_Phase = Cantera::newPhase("LiFeS_X_combo.xml", id);
 *    PhaseCombo_Interaction *LiFeS_X_solid = dynamic_cast <PhaseCombo_Interaction *>(l_tp);
 * @endcode
 *
 * or by the following constructor:
 *
 * @code
 *    XML_Node *xc = get_XML_File("LiFeS_X_combo.xml");
 *    XML_Node * const xs = xc->findNameID("phase", "LiFeS_X");
 *    PhaseCombo_Interaction *LiFeS_X_solid = new PhaseCombo_Interaction(*xs);
 * @endcode
 *
 * ## XML Example
 *
 * An example of an XML Element named phase setting up a PhaseCombo_Interaction
 * object named LiFeS_X  is given below.
 *
 * @code
 * <phase dim="3" id="LiFeS_X">
 *   <elementArray datasrc="elements.xml">
 *      Li Fe S
 *   </elementArray>
 *   <speciesArray datasrc="#species_LiFeS">
 *     LiTFe1S2(S)  Li2Fe1S2(S)
 *   </speciesArray>
 *   <thermo model="PhaseCombo_Interaction">
 *     <activityCoefficients model="Margules" TempModel="constant">
 *       <binaryNeutralSpeciesParameters speciesA="LiTFe1S2(S)" speciesB="Li2Fe1S2(S)">
 *         <excessEnthalpy model="poly_Xb" terms="2" units="kJ/mol">
 *             84.67069219, -269.1959421
 *         </excessEnthalpy>
 *         <excessEntropy  model="poly_Xb" terms="2" units="J/mol/K">
 *             100.7511565, -361.4222659
 *         </excessEntropy>
 *         <excessVolume_Enthalpy model="poly_Xb" terms="2" units="ml/mol">
 *              0, 0
 *         </excessVolume_Enthalpy>
 *         <excessVolume_Entropy  model="poly_Xb" terms="2" units="ml/mol/K">
 *              0, 0
 *         </excessVolume_Entropy>
 *       </binaryNeutralSpeciesParameters>
 *     </activityCoefficients>
 *   </thermo>
 *   <transport model="none"/>
 *   <kinetics model="none"/>
 * </phase>
 * @endcode
 *
 * The model attribute "PhaseCombo_Interaction" of the thermo XML element
 * identifies the phase as being of the type handled by the
 * PhaseCombo_Interaction object.
 *
 * @ingroup thermoprops
 */
class PhaseCombo_Interaction : public GibbsExcessVPSSTP
{
public:
    //! Constructor
    PhaseCombo_Interaction();

    //! Construct and initialize a PhaseCombo_Interaction ThermoPhase object
    //! directly from an XML input file
    /*!
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    PhaseCombo_Interaction(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize a PhaseCombo_Interaction ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    PhaseCombo_Interaction(XML_Node& phaseRef, const std::string& id = "");

    //! @name  Utilities
    //! @{

    virtual std::string type() const {
        return "PhaseCombo_Interaction";
    }

    //! @}
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
     * The activity \f$a_k\f$ of a species in solution is related to the
     * chemical potential by \f[ \mu_k = \mu_k^0(T) + \hat R T \log a_k. \f] The
     * quantity \f$\mu_k^0(T,P)\f$ is the chemical potential at unit activity,
     * which depends only on temperature and pressure.
     * @{
     */

    virtual void getActivityCoefficients(doublereal* ac) const;

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
     * \f[
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     * \f]
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
     * \f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     * \f]
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

    //! Return an array of partial molar volumes for the species in the mixture.
    //! Units: m^3/kmol.
    /*!
     * Frequently, for this class of thermodynamics representations, the excess
     * Volume due to mixing is zero. Here, we set it as a default. It may be
     * overridden in derived classes.
     *
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
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
    /// @{

    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

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
     * Margules model for a particular binary interaction. This function reads
     * the XML file and writes the coefficients it finds to an internal data
     * structures.
     *
     * @param xmlBinarySpecies  Reference to the XML_Node named
     *     "binaryNeutralSpeciesParameters" containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node& xmlBinarySpecies);

    //! Resize internal arrays within the object that depend upon the number of
    //! binary Margules interaction terms
    /*!
     *  @param num Number of binary Margules interaction terms
     */
    void resizeNumInteractions(const size_t num);

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

    //! Update the derivative of the log of the activity coefficients wrt
    //! log(mole fraction)
    /*!
     * This function will be called to update the internally stored derivative
     * of the natural logarithm of the activity coefficients wrt logarithm of
     * the mole fractions.
     */
    void s_update_dlnActCoeff_dlnX_diag() const;

    //! Update the derivative of the log of the activity coefficients wrt
    //! log(moles) - diagonal only
    /*!
     * This function will be called to update the internally stored diagonal
     * entries for the derivative of the natural logarithm of the activity
     * coefficients wrt logarithm of the moles.
     */
    void s_update_dlnActCoeff_dlnN_diag() const;

    //! Update the derivative of the log of the activity coefficients wrt
    //! log(moles_m)
    /*!
     * This function will be called to update the internally stored derivative
     * of the natural logarithm of the activity coefficients wrt logarithm of
     * the mole number of species
     */
    void s_update_dlnActCoeff_dlnN() const;

protected:
    //! number of binary interaction expressions
    size_t numBinaryInteractions_;

    //! Enthalpy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_HE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_HE_c_ij;

    //! Enthalpy term for the quaternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_HE_d_ij;

    //! Entropy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_SE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_SE_c_ij;

    //! Entropy term for the quaternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_SE_d_ij;

    //! Enthalpy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_VHE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_VHE_c_ij;

    //! Enthalpy term for the quaternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_VHE_d_ij;

    //! Entropy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_VSE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_VSE_c_ij;

    //! Entropy term for the quaternary mole fraction interaction of the excess
    //! Gibbs free energy expression
    mutable vector_fp m_VSE_d_ij;

    //! vector of species indices representing species A in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and
     *  B. This vector identifies species A.
     */
    std::vector<size_t> m_pSpecies_A_ij;

    //! vector of species indices representing species B in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and
     *  B. This vector identifies species B.
     */
    std::vector<size_t> m_pSpecies_B_ij;

    //! form of the Margules interaction expression
    /*!
     *  Currently there is only one form.
     */
    int formMargules_;

    //! form of the temperature dependence of the Margules interaction expression
    /*!
     *  Currently there is only one form -> constant wrt temperature.
     */
    int formTempModel_;
};

}

#endif
