/**
 *  @file IdealMolalSoln.h
 *   ThermoPhase object for the ideal molal equation of
 * state (see @ref thermoprops
 * and class @link Cantera::IdealMolalSoln IdealMolalSoln@endlink).
 *
 * Header file for a derived class of ThermoPhase that handles variable pressure
 * standard state methods for calculating thermodynamic properties that are
 * further based upon activities on the molality scale. The Ideal molal solution
 * assumes that all molality-based activity coefficients are equal to one. This
 * turns out to be highly nonlinear in the limit of the solvent mole fraction
 * going to zero.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALMOLALSOLN_H
#define CT_IDEALMOLALSOLN_H

#include "MolalityVPSSTP.h"

namespace Cantera
{

/**
 * This phase is based upon the mixing-rule assumption that all molality-based
 * activity coefficients are equal to one.
 *
 * This is a full instantiation of a ThermoPhase object. The assumption is that
 * the molality-based activity coefficient is equal to one. This also implies
 * that the osmotic coefficient is equal to one.
 *
 * Note, this does not mean that the solution is an ideal solution. In fact,
 * there is a singularity in the formulation as the solvent concentration goes
 * to zero.
 *
 * The mechanical equation of state is currently assumed to be that of an
 * incompressible solution. This may change in the future. Each species has its
 * own molar volume. The molar volume is a constant.
 *
 * Class IdealMolalSoln represents a condensed phase. The phase and the pure
 * species phases which comprise the standard states of the species are assumed
 * to have zero volume expansivity and zero isothermal compressibility. Each
 * species does, however, have constant but distinct partial molar volumes equal
 * to their pure species molar volumes. The class derives from class
 * ThermoPhase, and overloads the virtual methods defined there with ones that
 * use expressions appropriate for incompressible mixtures.
 *
 * The standard concentrations can have three different forms.
 * See setStandardConcentrationModel().
 *
 * @f$ V^0_0 @f$ is the solvent standard molar volume. @f$ m^{\Delta} @f$ is a
 * constant equal to a molality of @f$ 1.0 \quad\mbox{gm kmol}^{-1} @f$.
 *
 * The current default is to have mformGC = 2.
 *
 * The value and form of the activity concentration will affect reaction rate
 * constants involving species in this phase.
 *
 * An example phase definition is given in the
 * <a href="../../sphinx/html/yaml/phases.html#ideal-molal-solution">
 * YAML API Reference</a>.
 *
 * @ingroup thermoprops
 */
class IdealMolalSoln : public MolalityVPSSTP
{
public:
    //! Constructor for phase initialization
    /*!
     * This constructor will initialize a phase, by reading the required
     * information from an input file.
     *
     *  @param inputFile   Name of the Input file that contains information
     *      about the phase. If blank, an empty phase will be created.
     *  @param id          id of the phase within the input file
     */
    explicit IdealMolalSoln(const string& inputFile="", const string& id="");

    string type() const override {
        return "ideal-molal-solution";
    }

    bool isIdeal() const override {
        return true;
    }

    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    //! Molar enthalpy of the solution. Units: J/kmol.
    /*!
     * Returns the amount of enthalpy per mole of solution. For an ideal molal
     * solution,
     * @f[
     * \bar{h}(T, P, X_k) = \sum_k X_k \bar{h}_k(T)
     * @f]
     * The formula is written in terms of the partial molar enthalpies.
     * @f$ \bar{h}_k(T, p, m_k) @f$.
     * See the partial molar enthalpy function, getPartialMolarEnthalpies(),
     * for details.
     *
     * Units: J/kmol
     */
    double enthalpy_mole() const override;

    //! Molar internal energy of the solution: Units: J/kmol.
    /*!
     * Returns the amount of internal energy per mole of solution. For an ideal
     * molal solution,
     * @f[
     * \bar{u}(T, P, X_k) = \sum_k X_k \bar{u}_k(T)
     * @f]
     * The formula is written in terms of the partial molar internal energy.
     * @f$ \bar{u}_k(T, p, m_k) @f$.
     */
    double intEnergy_mole() const override;

    //! Molar entropy of the solution. Units: J/kmol/K.
    /*!
     * Returns the amount of entropy per mole of solution. For an ideal molal
     * solution,
     * @f[
     * \bar{s}(T, P, X_k) = \sum_k X_k \bar{s}_k(T)
     * @f]
     * The formula is written in terms of the partial molar entropies.
     * @f$ \bar{s}_k(T, p, m_k) @f$.
     * See the partial molar entropies function, getPartialMolarEntropies(),
     * for details.
     *
     * Units: J/kmol/K.
     */
    double entropy_mole() const override;

    //! Molar Gibbs function for the solution: Units J/kmol.
    /*!
     * Returns the Gibbs free energy of the solution per mole of the solution.
     *
     * @f[
     * \bar{g}(T, P, X_k) = \sum_k X_k \mu_k(T)
     * @f]
     *
     * Units: J/kmol
     */
    double gibbs_mole() const override;

    //! Molar heat capacity of the solution at constant pressure. Units: J/kmol/K.
    /*!
     * @f[
     * \bar{c}_p(T, P, X_k) = \sum_k X_k \bar{c}_{p,k}(T)
     * @f]
     *
     * Units: J/kmol/K
     */
    double cp_mole() const override;

    //! @}
    //! @name Mechanical Equation of State Properties
    //!
    //! In this equation of state implementation, the density is a function only
    //! of the mole fractions. Therefore, it can't be an independent variable.
    //! Instead, the pressure is used as the independent variable. Functions
    //! which try to set the thermodynamic state by calling setDensity() will
    //! cause an exception to be thrown.
    //! @{

protected:
    void calcDensity() override;

public:
    //! The isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * @f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * @f]
     *
     * It's equal to zero for this model, since the molar volume doesn't change
     * with pressure or temperature.
     */
    double isothermalCompressibility() const override;

    //! The thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     *
     * @f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * @f]
     *
     * It's equal to zero for this model, since the molar volume doesn't change
     * with pressure or temperature.
     */
    double thermalExpansionCoeff() const override;

    //! @}
    //! @name Activities and Activity Concentrations
    //!
    //! The activity @f$ a_k @f$ of a species in solution is related to the
    //! chemical potential by @f[ \mu_k = \mu_k^0(T) + \hat R T \ln a_k. @f] The
    //! quantity @f$ \mu_k^0(T) @f$ is the chemical potential at unit activity,
    //! which depends only on temperature and the pressure.
    //! @{

    Units standardConcentrationUnits() const override;
    void getActivityConcentrations(double* c) const override;
    double standardConcentration(size_t k=0) const override;

    /**
     * Get the array of non-dimensional activities at the current solution
     * temperature, pressure, and solution concentration.
     *
     * (note solvent is on molar scale)
     *
     * @param ac      Output activity coefficients. Length: m_kk.
     */
    void getActivities(double* ac) const override;

    /**
     * Get the array of non-dimensional molality-based activity coefficients at
     * the current solution temperature, pressure, and solution concentration.
     *
     * (note solvent is on molar scale. The solvent molar
     *  based activity coefficient is returned).
     *
     * @param acMolality      Output Molality-based activity coefficients.
     *                        Length: m_kk.
     */
    void getMolalityActivityCoefficients(double* acMolality) const override;

    //! @}
    //! @name  Partial Molar Properties of the Solution
    //! @{

    //!Get the species chemical potentials: Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solution.
     *
     * @f[
     *    \mu_k = \mu^{o}_k(T,P) + R T \ln(\frac{m_k}{m^\Delta})
     * @f]
     * @f[
     *    \mu_w = \mu^{o}_w(T,P) +
     *            R T ((X_w - 1.0) / X_w)
     * @f]
     *
     * @f$ w @f$ refers to the solvent species.
     * @f$ X_w @f$ is the mole fraction of the solvent.
     * @f$ m_k @f$ is the molality of the kth solute.
     * @f$ m^\Delta @f$ is 1 gmol solute per kg solvent.
     *
     * Units: J/kmol.
     *
     * @param mu     Output vector of species chemical potentials. Length: m_kk.
     */
    void getChemPotentials(double* mu) const override;

    //! Returns an array of partial molar enthalpies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol). For this phase, the partial molar enthalpies are equal to
     * the species standard state enthalpies.
     *  @f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * @f]
     * The reference-state pure-species enthalpies, @f$ \hat h^{ref}_k(T) @f$,
     * at the reference pressure,@f$ P_{ref} @f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param hbar   Output vector of partial molar enthalpies.
     *               Length: m_kk.
     */
    void getPartialMolarEnthalpies(double* hbar) const override;

    //! Returns an array of partial molar internal energies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol). For this phase, the partial molar internal energies are equal to
     * the species standard state internal energies (which are equal to the reference
     * state internal energies)
     *  @f[
     * \bar u_k(T,P) = \hat u^{ref}_k(T)
     * @f]
     * @param hbar   Output vector of partial molar internal energies, length #m_kk
     */
    void getPartialMolarIntEnergies(double* hbar) const override;

    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol.
    /*!
     * Maxwell's equations provide an insight in how to calculate this
     * (p.215 Smith and Van Ness)
     * @f[
     *      \frac{d(\mu_k)}{dT} = -\bar{s}_i
     * @f]
     * For this phase, the partial molar entropies are equal to the standard
     * state species entropies plus the ideal molal solution contribution.
     *
     * @f[
     *   \bar{s}_k(T,P) =  s^0_k(T) - R \ln( \frac{m_k}{m^{\triangle}} )
     * @f]
     * @f[
     *   \bar{s}_w(T,P) =  s^0_w(T) - R ((X_w - 1.0) / X_w)
     * @f]
     *
     * The subscript, w, refers to the solvent species. @f$ X_w @f$ is the mole
     * fraction of solvent. The reference-state pure-species entropies,@f$
     * s^0_k(T) @f$, at the reference pressure, @f$ P_{ref} @f$, are computed by
     * the species thermodynamic property manager. They are polynomial functions
     * of temperature.
     * @see MultiSpeciesThermo
     *
     * @param sbar Output vector of partial molar entropies.
     *             Length: m_kk.
     */
    void getPartialMolarEntropies(double* sbar) const override;

    // partial molar volumes of the species Units: m^3 kmol-1.
    /*!
     * For this solution, the partial molar volumes are equal to the constant
     * species molar volumes.
     *
     * Units: m^3 kmol-1.
     *  @param vbar Output vector of partial molar volumes.
     */
    void getPartialMolarVolumes(double* vbar) const override;

    //! Partial molar heat capacity of the solution:. UnitsL J/kmol/K
    /*!
     * The kth partial molar heat capacity is equal to the temperature
     * derivative of the partial molar enthalpy of the kth species in the
     * solution at constant P and composition (p. 220 Smith and Van Ness).
     * @f[
     *    \bar{Cp}_k(T,P) =  {Cp}^0_k(T)
     * @f]
     *
     * For this solution, this is equal to the reference state heat capacities.
     *
     * Units: J/kmol/K
     *
     * @param cpbar  Output vector of partial molar heat capacities.
     *               Length: m_kk.
     */
    void getPartialMolarCp(double* cpbar) const override;

    //! @}

    // -------------- Utilities -------------------------------

    bool addSpecies(shared_ptr<Species> spec) override;

    void initThermo() override;

    void getParameters(AnyMap& phaseNode) const override;

    //! Set the standard concentration model.
    /*!
     * Must be one of 'unity', 'species-molar-volume', or 'solvent-molar-volume'.
     * The default is 'solvent-molar-volume'.
     *
     * | model                | ActivityConc                     | StandardConc       |
     * | -------------------- | -------------------------------- | ------------------ |
     * | unity                | @f$ {m_k}/ { m^{\Delta}} @f$      | @f$ 1.0        @f$ |
     * | species-molar-volume | @f$  m_k / (m^{\Delta} V_k) @f$   | @f$ 1.0 / V_k  @f$ |
     * | solvent-molar-volume | @f$  m_k / (m^{\Delta} V^0_0) @f$ | @f$ 1.0 / V^0_0 @f$ |
     */
    void setStandardConcentrationModel(const string& model);

    //! Set cutoff model. Must be one of 'none', 'poly', or 'polyExp'.
    void setCutoffModel(const string& model);

    //! Report the molar volume of species k
    /*!
     * units - @f$ m^3 kmol^{-1} @f$
     *
     * @param k Species index.
     */
    double speciesMolarVolume(int k) const;

    /**
     * Fill in a return vector containing the species molar volumes
     * units - @f$ m^3 kmol^{-1} @f$
     *
     * @param smv Output vector of species molar volumes.
     */
    void getSpeciesMolarVolumes(double* smv) const;

protected:
    //! Species molar volume @f$ m^3 kmol^{-1} @f$
    vector<double> m_speciesMolarVolume;

    /**
     * The standard concentrations can have one of three different forms:
     * 0 = 'unity', 1 = 'species-molar-volume', 2 = 'solvent-molar-volume'. See
     * setStandardConcentrationModel().
     */
    int m_formGC = 2;

    //! Cutoff type
    int IMS_typeCutoff_ = 0;

private:
    //! vector of size m_kk, used as a temporary holding area.
    mutable vector<double> m_tmpV;

    //! Logarithm of the molal activity coefficients
    /*!
     *   Normally these are all one. However, stability schemes will change that
     */
    mutable vector<double> IMS_lnActCoeffMolal_;
public:
    //! value of the solute mole fraction that centers the cutoff polynomials
    //! for the cutoff =1 process;
    double IMS_X_o_cutoff_;

    //! gamma_o value for the cutoff process at the zero solvent point
    double IMS_gamma_o_min_;

    //! gamma_k minimum for the cutoff process at the zero solvent point
    double IMS_gamma_k_min_;

    //! Parameter in the polyExp cutoff treatment. This is the slope of the f
    //! function at the zero solvent point. Default value is 0.6
    double IMS_slopefCut_;

    //! Parameter in the polyExp cutoff treatment. This is the slope of the g
    //! function at the zero solvent point. Default value is 0.0
    double IMS_slopegCut_;

    //! @name Parameters in the polyExp cutoff having to do with rate of exp decay
    //! @{
    double IMS_cCut_;
    double IMS_dfCut_ = 0.0;
    double IMS_efCut_ = 0.0;
    double IMS_afCut_ = 0.0;
    double IMS_bfCut_ = 0.0;
    double IMS_dgCut_ = 0.0;
    double IMS_egCut_ = 0.0;
    double IMS_agCut_ = 0.0;
    double IMS_bgCut_ = 0.0;
    //! @}

private:
    //! This function will be called to update the internally stored
    //! natural logarithm of the molality activity coefficients
    /*!
     * Normally the solutes are all zero. However, sometimes they are not,
     * due to stability schemes.
     *
     *    gamma_k_molar =  gamma_k_molal / Xmol_solvent
     *
     *    gamma_o_molar = gamma_o_molal
     */
    void s_updateIMS_lnMolalityActCoeff() const;

    //! Calculate parameters for cutoff treatments of activity coefficients
    /*!
     * Some cutoff treatments for the activity coefficients actually require
     * some calculations to create a consistent treatment.
     *
     * This routine is called during the setup to calculate these parameters
     */
    void calcIMSCutoffParams_();
};

}

#endif
