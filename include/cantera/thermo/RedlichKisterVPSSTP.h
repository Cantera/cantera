/**
 *  @file  RedlichKisterVPSSTP.h (see @ref thermoprops and class @link
 *      Cantera::RedlichKisterVPSSTP RedlichKisterVPSSTP@endlink).
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
 * @f[
 *     G^E = \sum_{i} G^E_{i}
 * @f]
 *
 * where
 *
 * @f[
 *    G^E_{i} =   n X_{Ai} X_{Bi} \sum_m \left( A^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 * @f]
 *
 * where n is the total moles in the solution and where we can break down the Gibbs free
 * energy contributions into enthalpy and entropy contributions by defining
 * @f$ A^i_m = H^i_m - T S^i_m @f$ :
 *
 * @f[
 *    H^E_i = n X_{Ai} X_{Bi} \sum_m \left( H^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 * @f]
 *
 * @f[
 *    S^E_i = n X_{Ai} X_{Bi} \sum_m \left( S^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 * @f]
 *
 * The activity of a species defined in the phase is given by an excess Gibbs free
 * energy formulation:
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
 * @f[
 *      R T \ln( \gamma_k )=  \sum_i \delta_{Ai,k} (1 - X_{Ai}) X_{Bi} \sum_m \left( A^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 *                           + \sum_i \delta_{Ai,k} X_{Ai} X_{Bi} \sum_m \left(  A^{i}_0 +  A^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^{m-1} (1 - X_{Ai} + X_{Bi}) \right)
 * @f]
 *
 * Evaluating thermodynamic properties requires the following derivatives of
 * @f$ \ln(\gamma_k) @f$:
 *
 * @f[
 *    \frac{d \ln( \gamma_k )}{dT} = - \frac{1}{RT^2} \left( \sum_i \delta_{Ai,k} (1 - X_{Ai}) X_{Bi} \sum_m \left( H^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^m \right)
 *        + \sum_i \delta_{Ai,k} X_{Ai} X_{Bi} \sum_m \left(  H^{i}_0 +  H^{i}_m {\left( X_{Ai} -  X_{Bi} \right)}^{m-1} (1 - X_{Ai} + X_{Bi}) \right) \right)
 * @f]
 *
 * and
 *
 * @f[
 *    \frac{d^2 \ln( \gamma_k )}{dT^2} = -\frac{2}{T} \frac{d \ln( \gamma_k )}{dT}
 * @f]
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
 * The partial molar entropy for species *k* is given by the following relation,
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
 *              - R T^2 \frac{d^2 \ln(\gamma_k) }{{dT}^2} = C^o_{p,k}(T,P)
 * @f]
 *
 * @ingroup thermoprops
 */
class RedlichKisterVPSSTP : public GibbsExcessVPSSTP
{
public:
    //! Construct a RedlichKisterVPSSTP object from an input file
    /*!
     * @param inputFile Name of the input file containing the phase definition.
     *                  If blank, an empty phase will be created.
     * @param id        name (ID) of the phase in the input file. If empty, the
     *                  first phase definition in the input file will be used.
     */
    explicit RedlichKisterVPSSTP(const string& inputFile="", const string& id="");

    string type() const override {
        return "Redlich-Kister";
    }

    //! @name  Molar Thermodynamic Properties
    //! @{

    double cv_mole() const override;

    //! @}
    //! @name Activities, Standard States, and Activity Concentrations
    //!
    //! The activity @f$ a_k @f$ of a species in solution is
    //! related to the chemical potential by @f[ \mu_k = \mu_k^0(T)
    //! + \hat R T \ln a_k. @f] The quantity @f$ \mu_k^0(T,P) @f$ is
    //! the chemical potential at unit activity, which depends only
    //! on temperature and pressure.
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
     *  @f[
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *  @f]
     *
     * @param hbar  Vector of returned partial molar enthalpies
     *              (length m_kk, units = J/kmol)
     */
    void getPartialMolarEnthalpies(span<double> hbar) const override;

    //! Returns an array of partial molar entropies for the species in the
    //! mixture.
    /*!
     * For this phase, the partial molar entropies are equal to the standard
     * state entropies modified by the derivative of the activity coefficient
     * with respect to temperature:
     *
     *  @f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *  @f]
     *
     * @param sbar  Vector of returned partial molar entropies
     *              (length m_kk, units = J/kmol/K)
     */
    void getPartialMolarEntropies(span<double> sbar) const override;

    //! Returns an array of partial molar heat capacities for the species in the
    //! mixture.
    /*!
     * Units (J/kmol/K)
     *
     * For this phase, the first and second temperature derivative terms of
     * the activity coefficients cancel, so the partial molar heat capacities
     * are equal to the standard state heat capacities:
     *
     * @f[
     *      \tilde{C}_{p,k}(T,P) = C^o_{p,k}(T,P)
     * @f]
     *
     * @param cpbar  Vector of returned partial molar heat capacities
     *              (length m_kk, units = J/kmol/K)
     */
    void getPartialMolarCp(span<double> cpbar) const override;

    void getPartialMolarVolumes(span<double> vbar) const override;
    //! @}

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

    //! @name Initialization
    //!
    //! The following methods are used in the process of constructing
    //! the phase and setting its parameters from a specification in an
    //! input file. They are not normally used in application programs.
    //! To see how they are used, see importPhase().

    void initThermo() override;
    void getParameters(AnyMap& phaseNode) const override;

    //! Add a binary species interaction with the specified parameters
    /*!
     * @param speciesA         name of the first species
     * @param speciesB         name of the second species
     * @param excess_enthalpy  coefficients of the excess enthalpy polynomial
     * @param excess_entropy   coefficients of the excess entropy polynomial
     */
    void addBinaryInteraction(const string& speciesA, const string& speciesB,
        span<const double> excess_enthalpy, span<const double> excess_entropy);

    //! @name  Derivatives of Thermodynamic Variables needed for Applications
    //! @{

    void getdlnActCoeffds(const double dTds, span<const double> dXds,
                          span<double> dlnActCoeffds) const override;
    void getdlnActCoeffdlnX_diag(span<double> dlnActCoeffdlnX_diag) const override;
    void getdlnActCoeffdlnN_diag(span<double> dlnActCoeffdlnN_diag) const override;
    void getdlnActCoeffdlnN(const size_t ld, span<double> const dlnActCoeffdlnN) override;
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
     * the total derivatives (that is, not assuming other mole fractions are
     * constant) of the natural logarithm of the activity coefficients with
     * respect to the log of the mole fraction.
     */
    void s_update_dlnActCoeff_dlnX_diag() const;

protected:
    //! vector of species indices representing species A in the interaction
    /*!
     *  Each Redlich-Kister excess Gibbs free energy term involves two species,
     *  A and B. This vector identifies species A.
     */
    vector<size_t> m_pSpecies_A_ij;

    //! vector of species indices representing species B in the interaction
    /*!
     *  Each Redlich-Kister excess Gibbs free energy term involves two species,
     *  A and B. This vector identifies species B.
     */
    vector<size_t> m_pSpecies_B_ij;

    //! Enthalpy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    vector<vector<double>> m_HE_m_ij;

    //! Entropy term for the binary mole fraction interaction of the excess
    //! Gibbs free energy expression
    vector<vector<double>> m_SE_m_ij;

    //! Two dimensional array of derivatives of activity coefficients wrt mole
    //! fractions
    mutable Array2D dlnActCoeff_dX_;
};

}

#endif
