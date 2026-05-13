/**
 *  @file  MargulesVPSSTP.h (see @ref thermoprops and class @link
 *      Cantera::MargulesVPSSTP MargulesVPSSTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MARGULESVPSSTP_H
#define CT_MARGULESVPSSTP_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera
{

//! MargulesVPSSTP is a derived class of GibbsExcessVPSSTP that employs the
//! Margules approximation for the excess Gibbs free energy
/*!
 * MargulesVPSSTP derives from class GibbsExcessVPSSTP which is derived from
 * VPStandardStateTP, and overloads the virtual methods defined there with ones
 * that use expressions appropriate for the Margules Excess Gibbs free energy
 * approximation.
 *
 * The independent unknowns are pressure, temperature, and mass fraction.
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
 * 2 species.
 *
 * @f[
 *     G^E = \sum_i \left(  H_{Ei} - T S_{Ei} + (P - P_{ref}) V^E_i \right)
 * @f]
 * @f[
 *    H^E_i = n X_{Ai} X_{Bi} \left( h_{o,i} +  h_{1,i} X_{Bi} \right)
 * @f]
 * @f[
 *    S^E_i = n X_{Ai} X_{Bi} \left( s_{o,i} +  s_{1,i} X_{Bi} \right)
 * @f]
 * @f[
 *    V^E_i = n X_{Ai} X_{Bi} \left( v_{o,i} +  v_{1,i} X_{Bi} \right)
 * @f]
 *
 * where n is the total moles in the solution, @f$ P_{ref} @f$ is the
 * standard-state reference pressure (1 atm), and @f$ v_{o,i} = v^{HE}_{o,i} -
 * T \cdot v^{SE}_{o,i} @f$, @f$ v_{1,i} = v^{HE}_{1,i} - T \cdot
 * v^{SE}_{1,i} @f$ where @f$ v^{HE} @f$ and @f$ v^{SE} @f$ are the
 * `excess-volume-enthalpy` and `excess-volume-entropy` interaction parameters.
 * The @f$ (P - P_{\rm ref}) V^E_i @f$ term is what makes the Maxwell relation
 * @f$ (\partial S/\partial P)_T = -(\partial V/\partial T)_P @f$ consistent
 * with the non-zero, temperature-dependent excess volume.
 *
 * The activity of a species defined in the phase is given by an excess Gibbs
 * free energy formulation.
 *
 * @f[
 *      a_k = \gamma_k  X_k
 * @f]
 *
 * where
 *
 * @f[
 *      R T \ln( \gamma_k )= \frac{d(n G^E)}{d(n_k)}\Bigg|_{n_i}
 * @f]
 *
 * Taking the derivatives results in the following expression
 *
 * @f[
 *      R T \ln( \gamma_k )= \sum_i \left( \left( \delta_{Ai,k} X_{Bi} + \delta_{Bi,k} X_{Ai}  - X_{Ai} X_{Bi} \right)
 *       \left( g^E_{o,i} +  g^E_{1,i} X_{Bi} \right) +
 *       \left( \delta_{Bi,k} - X_{Bi} \right)      X_{Ai} X_{Bi}  g^E_{1,i} \right)
 * @f]
 * where
 * @f$  g^E_{o,i} =  h_{o,i} - T s_{o,i} + (P - P_{\rm ref})(v^{HE}_{o,i} - T \cdot v^{SE}_{o,i}) @f$
 * and
 * @f$ g^E_{1,i} =  h_{1,i} - T s_{1,i} + (P - P_{\rm ref})(v^{HE}_{1,i} - T \cdot v^{SE}_{1,i}) @f$
 * and where @f$ X_k @f$ is the mole fraction of species *k*.
 *
 * This object inherits from the class VPStandardStateTP. Therefore, the
 * specification and calculation of all standard state and reference state
 * values are handled at that level. Various functional forms for the standard
 * state are permissible. The chemical potential for species *k* is equal to
 *
 * @f[
 *      \mu_k(T,P) = \mu^o_k(T, P) + R T \ln(\gamma_k X_k)
 * @f]
 *
 * The partial molar entropy for species *k* is given by
 *
 * @f[
 *       \tilde{s}_k(T,P) =  s^o_k(T,P)  - R \ln( \gamma_k X_k )
 *              - R T \frac{d \ln(\gamma_k) }{dT}
 * @f]
 *
 * The partial molar enthalpy for species *k* is given by
 *
 * @f[
 *      \tilde{h}_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
 * @f]
 *
 * The partial molar volume for species *k* is
 *
 * @f[
 *        \tilde V_k(T,P)  = V^o_k(T,P)  + R T \frac{d \ln(\gamma_k) }{dP}
 * @f]
 *
 * The partial molar Heat Capacity for species *k* is
 *
 * @f[
 *      \tilde{C}_{p,k}(T,P) = C^o_{p,k}(T,P)   - 2 R T \frac{d \ln( \gamma_k )}{dT}
 *              - R T^2 \frac{d^2 \ln(\gamma_k) }{{dT}^2}
 * @f]
 *
 * @ingroup thermoprops
 */
class MargulesVPSSTP : public GibbsExcessVPSSTP
{
public:
    //! Construct a MargulesVPSSTP object from an input file
    /*!
     * @param inputFile Name of the input file containing the phase definition.
     *                  If blank, an empty phase will be created.
     * @param id        name (ID) of the phase in the input file. If empty, the
     *                  first phase definition in the input file will be used.
     */
    explicit MargulesVPSSTP(const string& inputFile="", const string& id="");

    string type() const override {
        return "Margules";
    }

    //! @name  Molar Thermodynamic Properties
    //! @{

    //! Molar heat capacity at constant volume.
    //!
    //! For a Margules phase with pressure-dependent standard state volumes, the
    //! thermodynamic identity @f$ c_v = c_p - T v \beta^2 / \kappa_T @f$ is used.
    //! If @f$ \kappa_T = 0 @f$ (all species have pressure-independent standard state
    //! volumes), then @f$ c_v = c_p @f$.
    double cv_mole() const override;

    //! Isothermal compressibility of the mixture.
    /*!
     * Computed as
     * @f[
     *   \kappa_T = -\frac{1}{v} \sum_k X_k \frac{\partial V^\circ_k}{\partial P}\Bigg|_T
     * @f]
     * The Margules excess volume has no pressure dependence, so only the standard
     * state PDSS derivatives contribute.
     */
    double isothermalCompressibility() const override;

    //! Thermal expansion coefficient of the mixture.
    /*!
     * Computed as
     * @f[
     *   \beta = \frac{1}{v} \left(
     *       \sum_k X_k \frac{\partial V^\circ_k}{\partial T}\Bigg|_P
     *       + \frac{\partial V_\mathrm{excess}}{\partial T}\Bigg|_{P,X}
     *   \right)
     * @f]
     * where the excess volume T-derivative is computed analytically from the
     * Margules excess-volume-entropy interaction parameters.
     */
    double thermalExpansionCoeff() const override;

    //! @}
    //! @name Activities, Standard States, and Activity Concentrations
    //!
    //! The activity @f$ a_k @f$ of a species in solution is related to the
    //! chemical potential by @f[ \mu_k = \mu_k^0(T) + \hat R T \ln a_k. @f] The
    //! quantity @f$ \mu_k^0(T,P) @f$ is the chemical potential at unit activity,
    //! which depends only on temperature and pressure.
    //! @{

    void getLnActivityCoefficients(span<double> lnac) const override;

    //! @}
    //! @name  Partial Molar Properties of the Solution
    //! @{

    void getChemPotentials(span<double> mu) const override;

    //! Returns an array of partial molar enthalpies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the molality-based
     * activity coefficient wrt temperature
     *
     * @f[
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     * @f]
     *
     * @param hbar  Vector of returned partial molar enthalpies
     *              (length m_kk, units = J/kmol)
     */
    void getPartialMolarEnthalpies(span<double> hbar) const override;

    //! Returns an array of partial molar entropies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the activity coefficient
     * wrt temperature
     *
     * @f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     * @f]
     *
     * @param sbar  Vector of returned partial molar entropies
     *              (length m_kk, units = J/kmol/K)
     */
    void getPartialMolarEntropies(span<double> sbar) const override;

    //! Returns an array of partial molar entropies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the standard
     * state enthalpies modified by the derivative of the activity coefficient
     * wrt temperature
     *
     *  @f[
     *   ???????????????
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *   ???????????????
     *  @f]
     *
     * @param cpbar  Vector of returned partial molar heat capacities
     *              (length m_kk, units = J/kmol/K)
     */
    void getPartialMolarCp(span<double> cpbar) const override;

    void getPartialMolarVolumes(span<double> vbar) const override;

    //! Get the array of temperature second derivatives of the log activity
    //! coefficients
    /*!
     *  units = 1/Kelvin
     *
     * @param d2lnActCoeffdT2  Output vector of temperature 2nd derivatives of
     *                         the log Activity Coefficients. length = m_kk
     */
    void getd2lnActCoeffdT2(span<double> d2lnActCoeffdT2) const;

    void getdlnActCoeffdT(span<double> dlnActCoeffdT) const override;

    //! @}
    //! @name Initialization
    //!
    //! The following methods are used in the process of constructing the phase
    //! and setting its parameters from a specification in an input file. They
    //! are not normally used in application programs. To see how they are used,
    //! see importPhase()
    //! @{

    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;

    //! Add a binary species interaction with the specified parameters
    /*!
     * @param speciesA   name of the first species
     * @param speciesB   name of the second species
     * @param h0         first excess enthalpy coefficient [J/kmol]
     * @param h1         second excess enthalpy coefficient [J/kmol]
     * @param s0         first excess entropy coefficient [J/kmol/K]
     * @param s1         second excess entropy coefficient [J/kmol/K]
     * @param vh0        first enthalpy coefficient for excess volume [m^3/kmol]
     * @param vh1        second enthalpy coefficient for excess volume [m^3/kmol]
     * @param vs0        first entropy coefficient for excess volume [m^3/kmol/K]
     * @param vs1        second entropy coefficient for excess volume [m^3/kmol/K]
     */
    void addBinaryInteraction(const string& speciesA,
        const string& speciesB, double h0, double h1, double s0, double s1,
        double vh0, double vh1, double vs0, double vs1);

    //! @}
    //! @name  Derivatives of Thermodynamic Variables needed for Applications
    //! @{

    void getdlnActCoeffds(const double dTds, span<const double> dXds,
                          span<double> dlnActCoeffds) const override;
    void getdlnActCoeffdlnX_diag(span<double> dlnActCoeffdlnX_diag) const override;
    void getdlnActCoeffdlnN_diag(span<double> dlnActCoeffdlnN_diag) const override;
    void getdlnActCoeffdlnN(const size_t ld,
                            span<double> const dlnActCoeffdlnN) override;

    //! @}

private:
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
    size_t numBinaryInteractions_ = 0;

    //! Enthalpy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_HE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_HE_c_ij;

    //! Entropy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_SE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_SE_c_ij;

    //! Enthalpy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_VHE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_VHE_c_ij;

    //! Entropy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_VSE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector<double> m_VSE_c_ij;

    //! vector of species indices representing species A in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and
     *  B. This vector identifies species A.
     */
    vector<size_t> m_pSpecies_A_ij;

    //! vector of species indices representing species B in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and
     *  B. This vector identifies species B.
     */
    vector<size_t> m_pSpecies_B_ij;

    //! form of the Margules interaction expression
    /*!
     *  Currently there is only one form.
     */
    int formMargules_ = 0;

    //! form of the temperature dependence of the Margules interaction expression
    /*!
     *  Currently there is only one form -> constant wrt temperature.
     */
    int formTempModel_ = 0;
};

}

#endif
