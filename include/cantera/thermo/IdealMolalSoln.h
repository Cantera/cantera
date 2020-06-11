/**
 *  @file IdealMolalSoln.h
 *   ThermoPhase object for the ideal molal equation of
 * state (see \ref thermoprops
 * and class \link Cantera::IdealMolalSoln IdealMolalSoln\endlink).
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
 * \f$ V^0_0 \f$ is the solvent standard molar volume. \f$ m^{\Delta} \f$ is a
 * constant equal to a molality of \f$ 1.0 \quad\mbox{gm kmol}^{-1} \f$.
 *
 * The current default is to have mformGC = 2.
 *
 * The value and form of the activity concentration will affect reaction rate
 * constants involving species in this phase.
 *
 *      <thermo model="IdealMolalSoln">
 *         <standardConc model="solvent_volume" />
 *         <solvent> H2O(l) </solvent>
 *         <activityCoefficients model="IdealMolalSoln" >
 *             <idealMolalSolnCutoff model="polyExp">
 *                 <gamma_O_limit> 1.0E-5  </gamma_O_limit>
 *                 <gamma_k_limit> 1.0E-5  <gamma_k_limit>
 *                 <X_o_cutoff>    0.20    </X_o_cutoff>
 *                 <C_0_param>     0.05    </C_0_param>
 *                 <slope_f_limit> 0.6     </slope_f_limit>
 *                 <slope_g_limit> 0.0     </slope_g_limit>
 *             </idealMolalSolnCutoff>
 *          </activityCoefficients>
 *      </thermo>
 *
 * @ingroup thermoprops
 */
class IdealMolalSoln : public MolalityVPSSTP
{
public:
    /// Constructor
    IdealMolalSoln();

    //! Constructor for phase initialization
    /*!
     * This constructor will initialize a phase, by reading the required
     * information from an input file.
     *
     *  @param inputFile   Name of the Input file that contains information
     *      about the phase
     *  @param id          id of the phase within the input file
     */
    IdealMolalSoln(const std::string& inputFile, const std::string& id = "");

    //! Constructor for phase initialization
    /*!
     * This constructor will initialize a phase, by reading the required
     * information from XML_Node tree.
     *
     *  @param phaseRef    reference for an XML_Node tree that contains
     *                     the information necessary to initialize the phase.
     *  @param id          id of the phase within the input file
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    IdealMolalSoln(XML_Node& phaseRef, const std::string& id = "");

    virtual std::string type() const {
        return "IdealMolalSoln";
    }

    //! @}
    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    //! Molar enthalpy of the solution. Units: J/kmol.
    /*!
     * Returns the amount of enthalpy per mole of solution. For an ideal molal
     * solution,
     * \f[
     * \bar{h}(T, P, X_k) = \sum_k X_k \bar{h}_k(T)
     * \f]
     * The formula is written in terms of the partial molar enthalpies.
     * \f$ \bar{h}_k(T, p, m_k) \f$.
     * See the partial molar enthalpy function, getPartialMolarEnthalpies(),
     * for details.
     *
     * Units: J/kmol
     */
    virtual doublereal enthalpy_mole() const;

    //! Molar internal energy of the solution: Units: J/kmol.
    /*!
     * Returns the amount of internal energy per mole of solution. For an ideal
     * molal solution,
     * \f[
     * \bar{u}(T, P, X_k) = \sum_k X_k \bar{u}_k(T)
     * \f]
     * The formula is written in terms of the partial molar internal energy.
     * \f$ \bar{u}_k(T, p, m_k) \f$.
     */
    virtual doublereal intEnergy_mole() const;

    //! Molar entropy of the solution. Units: J/kmol/K.
    /*!
     * Returns the amount of entropy per mole of solution. For an ideal molal
     * solution,
     * \f[
     * \bar{s}(T, P, X_k) = \sum_k X_k \bar{s}_k(T)
     * \f]
     * The formula is written in terms of the partial molar entropies.
     * \f$ \bar{s}_k(T, p, m_k) \f$.
     * See the partial molar entropies function, getPartialMolarEntropies(),
     * for details.
     *
     * Units: J/kmol/K.
     */
    virtual doublereal entropy_mole() const;

    //! Molar Gibbs function for the solution: Units J/kmol.
    /*!
     * Returns the Gibbs free energy of the solution per mole of the solution.
     *
     * \f[
     * \bar{g}(T, P, X_k) = \sum_k X_k \mu_k(T)
     * \f]
     *
     * Units: J/kmol
     */
    virtual doublereal gibbs_mole() const;

    //! Molar heat capacity of the solution at constant pressure. Units: J/kmol/K.
    /*!
     * \f[
     * \bar{c}_p(T, P, X_k) = \sum_k X_k \bar{c}_{p,k}(T)
     * \f]
     *
     * Units: J/kmol/K
     */
    virtual doublereal cp_mole() const;

    //@}
    /** @name Mechanical Equation of State Properties
     *
     * In this equation of state implementation, the density is a function only
     * of the mole fractions. Therefore, it can't be an independent variable.
     * Instead, the pressure is used as the independent variable. Functions
     * which try to set the thermodynamic state by calling setDensity() will
     * cause an exception to be thrown.
     */
    //@{

protected:
    /**
     * Calculate the density of the mixture using the partial molar volumes and
     * mole fractions as input
     *
     * The formula for this is
     *
     * \f[
     * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
     * \f]
     *
     * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are the molecular
     * weights, and \f$V_k\f$ are the pure species molar volumes.
     *
     * Note, the basis behind this formula is that in an ideal solution the
     * partial molar volumes are equal to the pure species molar volumes. We
     * have additionally specified in this class that the pure species molar
     * volumes are independent of temperature and pressure.
     */
    void calcDensity();

public:
    /**
     * Overridden setDensity() function is necessary because the density is not
     * an independent variable.
     *
     * This function will now throw an error condition
     *
     * @internal May have to adjust the strategy here to make the eos for these
     *     materials slightly compressible, in order to create a condition where
     *     the density is a function of the pressure.
     *
     * @param rho   Input Density
     * @deprecated Functionality merged with base function after Cantera 2.5.
     *             (superseded by isCompressible check in Phase::setDensity)
     */
    virtual void setDensity(const doublereal rho);

    /**
     * Overridden setMolarDensity() function is necessary because the density
     * is not an independent variable.
     *
     * This function will now throw an error condition.
     *
     * @param rho   Input Density
     * @deprecated Functionality merged with base function after Cantera 2.5.
     *             (superseded by isCompressible check in Phase::setDensity)
     */
    virtual void setMolarDensity(const doublereal rho);

    //! The isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *
     * It's equal to zero for this model, since the molar volume doesn't change
     * with pressure or temperature.
     */
    virtual doublereal isothermalCompressibility() const;

    //! The thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     *
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     *
     * It's equal to zero for this model, since the molar volume doesn't change
     * with pressure or temperature.
     */
    virtual doublereal thermalExpansionCoeff() const;

    /**
     * @}
     * @name Activities and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is related to the
     * chemical potential by \f[ \mu_k = \mu_k^0(T) + \hat R T \log a_k. \f] The
     * quantity \f$\mu_k^0(T)\f$ is the chemical potential at unit activity,
     * which depends only on temperature and the pressure.
     * @{
     */

    virtual Units standardConcentrationUnits() const;
    virtual void getActivityConcentrations(doublereal* c) const;
    virtual doublereal standardConcentration(size_t k=0) const;

    /*!
     * Get the array of non-dimensional activities at the current solution
     * temperature, pressure, and solution concentration.
     *
     * (note solvent is on molar scale)
     *
     * @param ac      Output activity coefficients. Length: m_kk.
     */
    virtual void getActivities(doublereal* ac) const;

    /*!
     * Get the array of non-dimensional molality-based activity coefficients at
     * the current solution temperature, pressure, and solution concentration.
     *
     * (note solvent is on molar scale. The solvent molar
     *  based activity coefficient is returned).
     *
     * @param acMolality      Output Molality-based activity coefficients.
     *                        Length: m_kk.
     */
    virtual void getMolalityActivityCoefficients(doublereal* acMolality) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

    //!Get the species chemical potentials: Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the species in
     * solution.
     *
     * \f[
     *    \mu_k = \mu^{o}_k(T,P) + R T \ln(\frac{m_k}{m^\Delta})
     * \f]
     * \f[
     *    \mu_w = \mu^{o}_w(T,P) +
     *            R T ((X_w - 1.0) / X_w)
     * \f]
     *
     * \f$ w \f$ refers to the solvent species.
     * \f$ X_w \f$ is the mole fraction of the solvent.
     * \f$ m_k \f$ is the molality of the kth solute.
     * \f$ m^\Delta \f$ is 1 gmol solute per kg solvent.
     *
     * Units: J/kmol.
     *
     * @param mu     Output vector of species chemical potentials. Length: m_kk.
     */
    virtual void getChemPotentials(doublereal* mu) const;

    //! Returns an array of partial molar enthalpies for the species in the
    //! mixture.
    /*!
     * Units (J/kmol). For this phase, the partial molar enthalpies are equal to
     * the species standard state enthalpies.
     *  \f[
     * \bar h_k(T,P) = \hat h^{ref}_k(T) + (P - P_{ref}) \hat V^0_k
     * \f]
     * The reference-state pure-species enthalpies, \f$ \hat h^{ref}_k(T) \f$,
     * at the reference pressure,\f$ P_{ref} \f$, are computed by the species
     * thermodynamic property manager. They are polynomial functions of
     * temperature.
     * @see MultiSpeciesThermo
     *
     * @param hbar   Output vector of partial molar enthalpies.
     *               Length: m_kk.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies of the species in the
    //! solution. Units: J/kmol.
    /*!
     * Maxwell's equations provide an insight in how to calculate this
     * (p.215 Smith and Van Ness)
     * \f[
     *      \frac{d(\mu_k)}{dT} = -\bar{s}_i
     * \f]
     * For this phase, the partial molar entropies are equal to the standard
     * state species entropies plus the ideal molal solution contribution.
     *
     * \f[
     *   \bar{s}_k(T,P) =  s^0_k(T) - R \ln( \frac{m_k}{m^{\triangle}} )
     * \f]
     * \f[
     *   \bar{s}_w(T,P) =  s^0_w(T) - R ((X_w - 1.0) / X_w)
     * \f]
     *
     * The subscript, w, refers to the solvent species. \f$ X_w \f$ is the mole
     * fraction of solvent. The reference-state pure-species entropies,\f$
     * s^0_k(T) \f$, at the reference pressure, \f$ P_{ref} \f$, are computed by
     * the species thermodynamic property manager. They are polynomial functions
     * of temperature.
     * @see MultiSpeciesThermo
     *
     * @param sbar Output vector of partial molar entropies.
     *             Length: m_kk.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    // partial molar volumes of the species Units: m^3 kmol-1.
    /*!
     * For this solution, the partial molar volumes are equal to the constant
     * species molar volumes.
     *
     * Units: m^3 kmol-1.
     *  @param vbar Output vector of partial molar volumes.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //! Partial molar heat capacity of the solution:. UnitsL J/kmol/K
    /*!
     * The kth partial molar heat capacity is equal to the temperature
     * derivative of the partial molar enthalpy of the kth species in the
     * solution at constant P and composition (p. 220 Smith and Van Ness).
     * \f[
     *    \bar{Cp}_k(T,P) =  {Cp}^0_k(T)
     * \f]
     *
     * For this solution, this is equal to the reference state heat capacities.
     *
     * Units: J/kmol/K
     *
     * @param cpbar  Output vector of partial molar heat capacities.
     *               Length: m_kk.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //@}

    // -------------- Utilities -------------------------------

    virtual bool addSpecies(shared_ptr<Species> spec);

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id="");

    virtual void initThermo();

    //! Set the standard concentration model.
    /*!
     * Must be one of 'unity', 'molar_volume', or 'solvent_volume'.
     * The default is 'solvent_volume'.
     *
     * | model          | ActivityConc                     | StandardConc       |
     * | -------------- | -------------------------------- | ------------------ |
     * | unity          | \f$ {m_k}/ { m^{\Delta}}\f$      | \f$ 1.0        \f$ |
     * | molar_volume   | \f$  m_k / (m^{\Delta} V_k)\f$   | \f$ 1.0 / V_k  \f$ |
     * | solvent_volume | \f$  m_k / (m^{\Delta} V^0_0)\f$ | \f$ 1.0 / V^0_0\f$ |
     */
    void setStandardConcentrationModel(const std::string& model);

    //! Set cutoff model. Must be one of 'none', 'poly', or 'polyExp'.
    void setCutoffModel(const std::string& model);

    //! Report the molar volume of species k
    /*!
     * units - \f$ m^3 kmol^{-1} \f$
     *
     * @param k Species index.
     */
    double speciesMolarVolume(int k) const;

    /*!
     * Fill in a return vector containing the species molar volumes
     * units - \f$ m^3 kmol^{-1} \f$
     *
     * @param smv Output vector of species molar volumes.
     */
    void getSpeciesMolarVolumes(double* smv) const;
    //@}

protected:
    //! Species molar volume \f$ m^3 kmol^{-1} \f$
    vector_fp m_speciesMolarVolume;

    /**
     * The standard concentrations can have one of three different forms:
     * 0 = 'unity', 1 = 'molar_volume', 2 = 'solvent_volume'. See
     * setStandardConcentrationModel().
     */
    int m_formGC;

    //! Cutoff type
    int IMS_typeCutoff_;

private:
    //! vector of size m_kk, used as a temporary holding area.
    mutable vector_fp m_tmpV;

    //! Logarithm of the molal activity coefficients
    /*!
     *   Normally these are all one. However, stability schemes will change that
     */
    mutable vector_fp IMS_lnActCoeffMolal_;
public:
    //! value of the solute mole fraction that centers the cutoff polynomials
    //! for the cutoff =1 process;
    doublereal IMS_X_o_cutoff_;

    //! gamma_o value for the cutoff process at the zero solvent point
    doublereal IMS_gamma_o_min_;

    //! gamma_k minimum for the cutoff process at the zero solvent point
    doublereal IMS_gamma_k_min_;

    //! Parameter in the polyExp cutoff treatment. This is the slope of the f
    //! function at the zero solvent point. Default value is 0.6
    doublereal IMS_slopefCut_;

    //! Parameter in the polyExp cutoff treatment. This is the slope of the g
    //! function at the zero solvent point. Default value is 0.0
    doublereal IMS_slopegCut_;

    //! @name Parameters in the polyExp cutoff treatment having to do with rate
    //!     of exp decay
    //! @{
    doublereal IMS_cCut_;
    doublereal IMS_dfCut_;
    doublereal IMS_efCut_;
    doublereal IMS_afCut_;
    doublereal IMS_bfCut_;
    doublereal IMS_dgCut_;
    doublereal IMS_egCut_;
    doublereal IMS_agCut_;
    doublereal IMS_bgCut_;
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
