/**
 * @file PlasmaPhase.h
 * Header file for class PlasmaPhase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAPHASE_H
#define CT_PLASMAPHASE_H

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class Reaction;
class ElectronCollisionPlasmaRate;

//! Base class for handling plasma properties, specifically focusing on the
//! electron energy distribution.
/*!
 * This class provides functionality to manage the the electron energy distribution
 * using two primary methods for defining the electron distribution and electron
 * temperature.
 *
 * The first method utilizes setElectronTemperature(), which sets the electron
 * temperature and calculates the electron energy distribution assuming an
 * isotropic-velocity model. Note that all units in PlasmaPhase are in SI, except
 * for electron energy, which is measured in volts.
 *
 * The generalized electron energy distribution for an isotropic-velocity
 * distribution (as described by Gudmundsson @cite gudmundsson2001 and Khalilpour
 * and Foroutan @cite khalilpour2020)
 * is given by:
 *   @f[
 *          f(\epsilon) = c_1 \frac{\sqrt{\epsilon}}{\epsilon_m^{3/2}}
 *          \exp \left(-c_2 \left(\frac{\epsilon}{\epsilon_m}\right)^x \right),
 *   @f]
 * where @f$ x = 1 @f$ corresponds to a Maxwellian distribution and
 * @f$ x = 2 @f$ corresponds to a Druyvesteyn distribution.
 * Here, @f$ \epsilon_m = \frac{3}{2} T_e @f$ [V] represents the
 * mean electron energy.
 *
 * The total probability distribution integrates to one:
 *   @f[
 *           \int_0^{\infty} f(\epsilon) d\epsilon = 1.
 *   @f]
 * According to Hagelaar and Pitchford @cite hagelaar2005, the electron energy
 * probability function can be defined as
 * @f$ F(\epsilon) = \frac{f(\epsilon)}{\sqrt{\epsilon}} @f$ with units of
 * [V@f$^{-3/2}@f$]. The generalized form of the electron energy probability
 * function for isotropic-velocity distributions is:
 *   @f[
 *          F(\epsilon) = c_1 \frac{1}{\epsilon_m^{3/2}}
 *          \exp\left(-c_2 \left(\frac{\epsilon}{\epsilon_m}\right)^x\right),
 *   @f]
 * and this form is used to model the isotropic electron energy distribution
 * in PlasmaPhase.
 *
 * The second method allows for manual definition of the electron energy
 * distribution using setDiscretizedElectronEnergyDist(). In this approach,
 * the electron temperature is derived from the mean electron energy,
 * @f$ \epsilon_m @f$, which can be calculated as follows @cite hagelaar2005 :
 *   @f[
 *          \epsilon_m = \int_0^{\infty} \epsilon^{3/2} F(\epsilon) d\epsilon.
 *   @f]
 * This integral can be approximated using the trapezoidal rule,
 *   @f[
 *          \epsilon_m = \sum_i \left(\epsilon_{i+1}^{5/2} - \epsilon_i^{5/2}\right)
 *                       \frac{F(\epsilon_{i+1}) + F(\epsilon_i)}{2},
 *   @f]
 * where @f$ i @f$ is the index of discrete energy levels, or Simpson's rule.
 *
 * ## Thermodynamic consistency of the two-temperature model
 *
 * When the electron temperature @f$ T_\text{e} @f$ differs from the heavy-species
 * temperature @f$ T @f$, this phase does not satisfy all of Cantera's standard-state
 * consistency relations simultaneously: one of them is deliberately left unsatisfied.
 * This is a consequence of two design commitments that this implementation makes, both
 * of which are individually reasonable but which cannot be reconciled once
 * @f$ T_\text{e} \neq T @f$:
 *
 * 1. **Kinetics uses true concentrations.** The generalized (activity) concentrations
 *    that drive the law of mass action are kept equal to the true molar concentrations,
 *    @f$ C^a_k = C_k @f$, so that reaction rates depend on actual number densities.
 *    Together with the standard concentration @f$ C^0_k = P/RT @f$ (the gas
 *    temperature, for all species; see #standardConcentration()), this requires the
 *    reported activities to be
 *    @f[ a_k = \frac{C^a_k}{C^0_k} = X_k \frac{T}{\overline{T}}, @f]
 *    where @f$ \overline{T} @f$ is the mole-fraction-weighted mean temperature (see
 *    #meanTemperature()). This is what #getActivities() returns.
 *
 * 2. **Chemical potentials keep the ideal-gas form.** The chemical potentials retain
 *    the usual ideal-gas composition dependence,
 *    @f$ \mu_k = \mu^0_k(T_k) + R T_k \ln X_k @f$,
 *    with each species evaluated at its own temperature @f$ T_k @f$ (@f$ T @f$ for
 *    heavy species, @f$ T_\text{e} @f$ for electrons; see #getChemPotentials() and
 *    #getStandardChemPotentials()). This is equivalent to an activity of @f$ X_k @f$.
 *
 * When @f$ T_\text{e} \neq T @f$, these two demands imply different values for the
 * activities -- @f$ X_k\,T/\overline{T} @f$ versus @f$ X_k @f$ -- which differ by the
 * factor @f$ T/\overline{T} @f$. #getActivities() reports the first;
 * #getChemPotentials() embodies the second; and the identity
 * @f$ \mu_k = \mu^0_k + R T_k \ln a_k @f$ is therefore the relation left unsatisfied.
 * This is why the consistency test `chem_potentials_to_activities` is a known failure
 * for two-temperature states.
 *
 * The practical consequence of this inconsistency is limited. It formally perturbs any
 * *chemical* equilibrium computed for the phase, but the concept of chemical
 * equilibrium in the absence of *thermal* equilibrium (@f$ T_\text{e} \neq T @f$) is
 * itself ill-posed, so the perturbed quantity has little physical meaning in exactly
 * the regime where the inconsistency appears. The trade-offs and potential approaches
 * for a fully-consistent implementation are documented in
 * https://github.com/Cantera/enhancements/issues/258#issuecomment-4857093158.
 *
 * @warning  This class is an experimental part of %Cantera and may be
 *           changed or removed without notice.
 * @todo Implement electron Boltzmann equation solver to solve EEDF.
 *       https://github.com/Cantera/enhancements/issues/127
 * @ingroup thermoprops
 */
class PlasmaPhase: public IdealGasPhase
{
public:
    //! Construct and initialize a PlasmaPhase object directly from an input file.
    //! The constructor initializes the electron energy distribution to be a
    //! Maxwellian distribution (#m_isotropicShapeFactor = 1.0).
    //! The initial electron energy grid is set to a linear space which starts
    //! at 0.01 eV and ends at 1 eV with 1000 points.
    /*!
     * @param  inputFile Name of the input file containing the phase definition
     *                   to set up the object. If blank, an empty phase will be
     *                   created.
     * @param  id        ID of the phase in the input file. Defaults to the
     *                   empty string.
     */
    explicit PlasmaPhase(const string& inputFile="", const string& id="");

    ~PlasmaPhase();

    string type() const override {
        return "plasma";
    }

    void initThermo() override;

    //! @name Overridden from IdealGasPhase or ThermoPhase
    //! @{
    bool addSpecies(shared_ptr<Species> spec) override;
    virtual void setSolution(std::weak_ptr<Solution> soln) override;
    void getParameters(AnyMap& phaseNode) const override;
    void setParameters(const AnyMap& phaseNode,
                       const AnyMap& rootNode=AnyMap()) override;
    //! @}
    //! @name Electron Species Information
    //! @{

    //! Electron Species Index
    size_t electronSpeciesIndex() const {
        return m_electronSpeciesIndex;
    }

    //! Electron species name
    string electronSpeciesName() const {
        return speciesName(m_electronSpeciesIndex);
    }

    //! @}
    //! @name Electron Energy Distribution Functions
    //! @{

    //! Set electron energy levels.
    //! @param  levels The vector of electron energy levels [eV].
    //!                Length: #m_nPoints.
    void setElectronEnergyLevels(span<const double> levels);

    //! Get electron energy levels.
    //! @param  levels The vector of electron energy levels [eV]. Length: #m_nPoints
    void getElectronEnergyLevels(span<double> levels) const {
        Eigen::Map<Eigen::ArrayXd>(levels.data(), levels.size()) = m_electronEnergyLevels;
    }

    //! Set discretized electron energy distribution.
    //! @param  levels The vector of electron energy levels [eV].
    //!                Length: #m_nPoints.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    void setDiscretizedElectronEnergyDist(span<const double> levels,
                                          span<const double> distrb);

    //! Get electron energy distribution.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    void getElectronEnergyDistribution(span<double> distrb) const {
        Eigen::Map<Eigen::ArrayXd>(distrb.data(), distrb.size()) = m_electronEnergyDist;
    }

    //! Set the shape factor of isotropic electron energy distribution.
    //! Note that @f$ x = 1 @f$ and @f$ x = 2 @f$ correspond to the
    //! Maxwellian and Druyvesteyn distribution, respectively.
    //! @param  x The shape factor
    void setIsotropicShapeFactor(double x);

    //! The shape factor of isotropic electron energy distribution.
    double isotropicShapeFactor() const {
        return m_isotropicShapeFactor;
    }

    //! Electron Temperature [K].
    //!     @return The electron temperature of the phase.
    double electronTemperature() const override {
        return m_electronTemp;
    }
    //! Set the internally stored electron temperature of the phase [K].
    //! @param  Te Electron temperature in Kelvin.
    void setElectronTemperature(double Te) override;

    //! Mean electron energy [eV].
    double meanElectronEnergy() const {
        return 3.0 / 2.0 * electronTemperature() * Boltzmann / ElectronCharge;
    }

    //! Set mean electron energy [eV].
    //! This method also sets electron temperature accordingly.
    void setMeanElectronEnergy(double energy);

    //! Get electron energy distribution type.
    string electronEnergyDistributionType() const {
        return m_distributionType;
    }

    //! Set electron energy distribution type.
    void setElectronEnergyDistributionType(const string& type);

    //! Numerical quadrature method. Method: #m_quadratureMethod
    string quadratureMethod() const {
        return m_quadratureMethod;
    }

    //! Set numerical quadrature method for integrating electron
    //! energy distribution function. Method: #m_quadratureMethod
    void setQuadratureMethod(const string& method) {
        m_quadratureMethod = method;
    }

    //! Set flag of automatically normalize electron energy distribution.
    //! Flag: #m_do_normalizeElectronEnergyDist
    void enableNormalizeElectronEnergyDist(bool enable) {
        m_do_normalizeElectronEnergyDist = enable;
    }

    //! Flag of automatically normalize electron energy distribution.
    //! Flag: #m_do_normalizeElectronEnergyDist
    bool normalizeElectronEnergyDistEnabled() const {
        return m_do_normalizeElectronEnergyDist;
    }

    //! Number of electron levels.
    size_t nElectronEnergyLevels() const {
        return m_nPoints;
    }

    //! Number of electron collision cross sections.
    size_t nCollisions() const {
        return m_collisions.size();
    }

    //! Get the Reaction object associated with electron collision *i*.
    //! @since New in %Cantera 3.2.
    const shared_ptr<Reaction> collision(size_t i) const {
        return m_collisions[i];
    }

    //! Get the ElectronCollisionPlasmaRate object associated with electron collision
    //! *i*.
    //! @since New in %Cantera 3.2.
    const shared_ptr<ElectronCollisionPlasmaRate> collisionRate(size_t i) const {
        return m_collisionRates[i];
    }

    //! Update the electron energy distribution.
    void updateElectronEnergyDistribution();

    //! Return the distribution number #m_distNum.
    int distributionNumber() const {
        return m_distNum;
    }

    //! Return the electron energy level number #m_levelNum.
    int levelNumber() const {
        return m_levelNum;
    }

    //! Get the indicies for inelastic electron collisions.
    //! @since New in %Cantera 3.2.
    const vector<size_t>& kInelastic() const {
        return m_kInelastic;
    }

    //! Get the indices for elastic electron collisions.
    //! @since New in %Cantera 3.2.
    const vector<size_t>& kElastic() const {
        return m_kElastic;
    }

    //! Return the target of a specific process.
    //! @since New in %Cantera 3.2.
    size_t targetIndex(size_t i) const {
        return m_targetSpeciesIndices[i];
    }

    //! Get the frequency of the applied electric field [Hz].
    //! @since New in %Cantera 3.2.
    double electricFieldFrequency() const {
        return m_electricFieldFrequency;
    }

    //! Get the applied electric field strength [V/m].
    double electricField() const {
        return m_electricField;
    }

    //! Set the absolute electric field strength [V/m].
    void setElectricField(double E) {
        m_electricField = E;
    }

    //! Calculate the degree of ionization
    //double ionDegree() const {
    //    double ne = concentration(m_electronSpeciesIndex); // [kmol/m³]
    //    double n_total = molarDensity();                   // [kmol/m³]
    //    return ne / n_total;
    //}

    //! Get the reduced electric field strength [V·m²].
    double reducedElectricField() const {
        return m_electricField / (molarDensity() * Avogadro);
    }

    //! Set reduced electric field given in [V·m²].
    void setReducedElectricField(double EN) {
        m_electricField = EN * molarDensity() * Avogadro; // [V/m]
    }

    /**
     * The electron mobility (m²/V/s)
     *   @f[
     *     \mu = \nu_d / E,
     *   @f]
     * where @f$ \nu_d @f$ is the drift velocity (m²/s), and @f$ E @f$ is the electric
     * field strength (V/m).
     */
    double electronMobility() const;

    /**
     * The elastic power loss [J/s/m³]
     *   @f[
     *     P_k = N_A N_A C_e e \sum_k C_k K_k,
     *   @f]
     * where @f$ C_k @f$ and @f$ C_e @f$ are the concentration (kmol/m³) of the
     * target species and electrons, respectively. @f$ K_k @f$ is the elastic
     * electron energy loss coefficient (eV-m³/s).
     */
    double elasticPowerLoss();

    /**
     * The joule heating power (W/m³)
     *   @f[
     *     q_J = \sigma * E^2,
     *   @f]
     * where @f$ \sigma @f$ is the conductivity (S/m), defined by:
     *   @f[
     *     \sigma = e * n_e * \mu_e
     *   @f]
     * and @f$ E @f$ is the electric field strength (V/m).
     */
    double jouleHeatingPower() const;

    void beginEquilibrate() override;

    void endEquilibrate() override;

    double intrinsicHeating() override;

    //! Set parameters for the electron energy distribution.
    /*!
    * This method configures the electron energy distribution using the same
    * mapping as the YAML `electron-energy-distribution` entry for a plasma phase.
    *
    * The mapping must contain the key `type`, which selects the electron energy
    * distribution model. Supported values are:
    *
    * - `isotropic`
    * - `discretized`
    * - `Boltzmann-two-term`
    *
    * For `isotropic`, the mapping must contain:
    *
    * - `shape-factor`
    * - `mean-electron-energy`
    *
    * The optional key `energy-levels` may be used to specify the energy grid.
    *
    * For `discretized`, the mapping must contain:
    *
    * - `energy-levels`
    * - `distribution`
    *
    * The optional key `normalize` controls whether the distribution is normalized.
    *
    * For `Boltzmann-two-term`, the mapping may specify either:
    *
    * - `energy-levels`, defining an explicit fixed grid; or
    * - `initial-max-energy-level` and `grid-cell-count`, defining a generated grid.
    *
    * Generated grids may also use:
    *
    * - `energy-level-spacing`, with supported values `linear`, `quadratic`,
    *   and `geometric`. If omitted, `linear` is used.
    * - `geometric-grid-ratio`, used only for `geometric` grids.
    * - `energy-grid-adaptation`, a mapping controlling automatic grid adaptation (see
    *   the [`plasma` input documentation](../reference/yaml/phases.html) for further
    *   details on energy-grid adaptation).
    *
    * The optional key `reduced-field-threshold-before-Maxwellian` applies to
    * `Boltzmann-two-term` for both explicit and generated grids: below this threshold,
    * a gas-temperature Maxwellian is imposed.
    *
    * This method updates the electron energy levels and electron energy
    * distribution stored by `PlasmaPhase`. Changing the EEDF model or grid
    * invalidates the previously stored EEDF.
    *
    * @param eedf  Mapping using the same keys as the YAML
    *              `electron-energy-distribution` entry.
    *
    * @since New in Cantera 4.0
    */
    void setElectronEnergyDistributionParameters(const AnyMap& eedf);


    //! @}
    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    //! Return the Molar enthalpy. Units: J/kmol.
    /*!
     * For an ideal gas mixture with electrons at a different temperature
     * than the heavy species, the molar enthalpy is calculated as:
     * @f[
     * \begin{align}
     *  \hat h(T, T_\text{e})
     *  &= \sum_{k} X_k h^0_k(T_k) \\
     *  &= \sum_{k \neq e} X_k h^0_k(T) + X_\text{e} h_\text{e}^0(T_\text{e})
     * \end{align}
     * @f]
     * where heavy-species properties are evaluated at @f$ T @f$, electron properties
     * at @f$ T_\text{e} @f$, and @f$ P @f$ is the total pressure of the mixture.
     * For an ideal gas, enthalpy is independent of pressure.
     * The standard-state pure-species enthalpies @f$ h^0_k(T_k) @f$ are
     * computed by the species thermodynamic property manager.
     */
    double enthalpy_mole() const override;

    //! Return the molar entropy. Units: J/kmol/K.
    /*!
     * For an ideal gas mixture with electrons at a different temperature
     * than the heavy species, the molar entropy is calculated as:
     * @f[
     * \hat s(T,T_\text{e},P) = \sum_k X_k \hat s^0_k(T_k) - \hat R \ln \frac{P}{P^0}
     * @f]
     * where heavy-species properties are evaluated at @f$ T @f$, electron properties
     * at @f$ T_\text{e} @f$, and @f$ P @f$ is the total pressure of the mixture.
     */
    double entropy_mole() const override;

    //! Return the molar Gibbs free energy. Units: J/kmol.
    /*!
     * For an ideal gas mixture with electrons at a different temperature
     * than the heavy species, the molar Gibbs free energy is calculated as:
     * @f[
     * \begin{align}
     *  \hat g(T,T_\text{e},P)
     *  &= \sum_{k} X_k \tilde{g}_k(T_k, P)
     *   = \sum_{k} X_k \mu_k(T_k, P) \\
     *  &= \sum_{k} X_k \left[\mu^0_k(T_k)
     *     + R T_k \ln \left( \frac{X_k P}{P^0} \right) \right] \\
     * \end{align}
     * @f]
     * where heavy-species properties are evaluated at @f$ T @f$, electron properties
     * at @f$ T_\text{e} @f$, and @f$ P @f$ is the total pressure of the mixture.
     */
    double gibbs_mole() const override;

    //! Return the molar internal energy. Units: J/kmol.
    /*!
     * For an ideal gas mixture with electrons at a different temperature
     * than the heavy species, the molar internal energy is calculated as:
     * @f[
     * \hat u(T,T_\text{e},P)
     *  = \sum_{k} X_k \tilde{u}_k(T_k)
     *  = \sum_{k} X_k (\tilde{h}_k(T_k) - R T_k)
     * @f]
     * where heavy-species properties are evaluated at @f$ T @f$, electron properties
     * at @f$ T_\text{e} @f$, and @f$ P @f$ is the total pressure of the mixture.
     * For an ideal gas, internal energy is independent of pressure.
     */
    double intEnergy_mole() const override;

    //! @}
    //! @name Mechanical Equation of State
    //! @{

    //! Return the mean temperature of the plasma phase. Units: K.
    /*!
     * In a plasma phase, the electron temperature can differ from the
     * heavy-species (gas) temperature. Therefore, the mean temperature is
     * defined as a mole-fraction-weighted average of the electron and
     * heavy-species temperatures:
     * @f[
     * \begin{align}
     *   \overline{T} &= \sum_{k \neq e} X_k T + X_\text{e} T_\text{e} \\
     *                &= (1 - X_\text{e}) T + X_\text{e} T_\text{e} \\
     *                &= T + X_\text{e} (T_\text{e} - T)
     * \end{align}
     * @f]
     * See the #pressure() method for usage of the mean temperature in
     * calculating the pressure of the plasma phase.
     */
    double meanTemperature() const;

    //! Return the pressure of the plasma phase. Units: Pa.
    /*!
     * The pressure of the plasma phase is calculated using the mean temperature,
     * which is a mole-fraction-weighted average of the electron and heavy-species
     * temperatures.
     * @f[
     * \begin{align}
     *      P &= \sum_k n_k k_\text{B} T_k \\
     *        &= \sum_{k \neq e} n_k k_\text{B} T
     *           + n_\text{e} k_\text{B} T_\text{e} \\
     *        &= (n_\text{total} - n_\text{e}) k_\text{B} T
     *           + n_\text{e} k_\text{B} T_\text{e} \\
     *        &= n_\text{total} (1 - X_\text{e}) k_\text{B} T
     *           + n_\text{total} X_\text{e} k_\text{B} T_\text{e} \\
     *        &= n_\text{total} k_\text{B} \overline{T} \\
     * \end{align}
     * @f]
     * where @f$ \overline{T} @f$ is the mean temperature of the plasma phase,
     * defined in the #meanTemperature() method.
     * Here, @f$ n_k @f$ is the number density of species *k* [in 1/m³],
     * @f$ n_\text{total} @f$ is the total number density of the phase [in 1/m³],
     * @f$ X_\text{e} @f$ is the mole fraction of electrons,
     * @f$ T @f$ is the heavy-species (gas) temperature [K],
     * @f$ T_\text{e} @f$ is the electron temperature [K],
     * and @f$ k_\text{B} @f$ is the Boltzmann constant [J/K].
     *
     * The number density times Boltzmann constant can be expressed as
     * @f$ n_\text{total} k_\text{B} = C_\text{total} R @f$, where
     * @f$ C_\text{total} @f$ is the total molar concentration of the phase [kmol/m³]
     * and @f$ R @f$ is the gas constant [J/kmol/K], so that:
     * @f[
     *     P = C_\text{total} R \overline{T}.
     * @f]
     *
     * Here, the pressure is set via the density, via:
     * @f[
     *     P = R \overline{T} \frac{\rho}{\overline{W}}
     * @f]
     * where @f$ \overline{W} @f$ is the mean molecular weight.
     */
    double pressure() const override;

    //! Set the pressure at constant temperature and composition. Units: Pa.
    /*!
     * This method is implemented by setting the mass density to
     * @f[
     * \rho = \frac{P \overline{W}}{\hat R \overline{T} }.
     * @f]
     *
     * where @f$ \overline{W} @f$ is the mean molecular weight,
     * and  @f$ \overline{T} @f$ is the mean temperature (see #meanTemperature()).
     *
     * @param p Pressure [Pa].
     */
    void setPressure(double p) override {
        setDensity(p * meanMolecularWeight() / (GasConstant * meanTemperature()));
    }

    //! Return the Gas Constant multiplied by the current electron temperature [J/kmol].
    double RTe() const {
        return electronTemperature() * GasConstant;
    }

    //! Return the electron pressure. Units: Pa.
    /*!
     * The partial pressure of electrons is calculated using the electron temperature:
     * @f[
     *    P_\text{e} = n_\text{e} k_\text{B} T_\text{e} = C_\text{e} R T_\text{e},
     * @f]
     * where @f$ n_\text{e} @f$ is the particular density of electrons [in 1/m³],
     * @f$ k_\text{B} @f$ is the Boltzmann constant [J/K],
     * @f$ C_\text{e} @f$ is the molar concentration of electrons [in kmol/m³],
     * @f$ R @f$ is the gas constant [J/kmol/K], and
     * @f$ T_\text{e} @f$ is the electron temperature [K].
     */
    virtual double electronPressure() const {
        return GasConstant * concentration(m_electronSpeciesIndex) *
               electronTemperature();
    }

    //! Raise NotImplementedError, as there is an ambiguity on the temperature to use.
    double thermalExpansionCoeff() const override {
        throw NotImplementedError("PlasmaPhase::thermalExpansionCoeff");
    }

    //! Raise NotImplementedError, as there is an ambiguity on the temperature to use.
    double soundSpeed() const override {
        throw NotImplementedError("PlasmaPhase::soundSpeed");
    }

    //! @}
    //! @name Chemical Potentials and Activities
    //! @{

    //! Returns the standard concentration @f$ C^0_k @f$, which is used to
    //! normalize the generalized concentration. Units: m³/kmol.
    /*!
     * The standard concentration is chosen to be equal to @f[ C^0_k = \frac{P}{R T} @f]
     * where @f$ P @f$ is the total pressure of the mixture,
     * @f$ R @f$ is the gas constant, and @f$ T @f$ is the gas temperature.
     *
     * @note The gas temperature @f$ T @f$ (not the electron temperature or the
     *     mean temperature) is used for *all* species, including electrons.
     *     Consequently @f$ C^0_\text{e} = P/RT \neq 1/V^0_\text{e} = P/R T_\text{e} @f$
     *     when @f$ T_\text{e} \neq T @f$. This is intentional: the reciprocal
     *     relation @f$ C^0_k = 1/V^0_k @f$ is a coincidence of the
     *     single-temperature ideal gas and is not a thermodynamic requirement.
     *     See the class documentation.
     *
     * @param k Parameter indicating the species. The default
     *          is to assume this refers to species 0.
     */
    double standardConcentration(size_t k=0) const override;

    //! Get the array of non-dimensional activities at the current solution
    //! temperature, pressure, and solution concentration.
    /*!
     * Since @f$ a_k = C^a_k / C^0_k @f$, where @f$ C^a_k = C_k @f$ is the true molar
     * concentration, and @f$ C^0_k @f$ is the standard concentration defined above,
     * the activities are equal to @f$ a_k = X_k T / \overline{T} @f$.
     *
     * @note These activities are defined for consistency with the kinetics
     *     (law-of-mass-action) concentrations and, when @f$ T_\text{e} \neq T @f$,
     *     differ by a factor of @f$ T/\overline{T} @f$ from the activities implied by
     *     the chemical-potential model (@f$ a_k = X_k @f$). See the class documentation
     *     for a discussion of this intentional inconsistency.
     *
     * @param a   Output vector of activities. Length: m_kk.
     */
    void getActivities(span<double> a) const override;

    //! Get the array of non-dimensional activity coefficients at the current
    //! solution temperature, pressure, and solution concentration.
    /*!
     * The activity coefficients are defined by @f$ \gamma_k = a_k / X_k @f$.
     * With the current definition of the activities, the activity coefficients are
     * equal to @f$ \gamma_k = T / \overline{T} @f$.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    void getActivityCoefficients(span<double> ac) const override;

    //! @}
    //! @name Partial Molar Properties of the Solution
    //! @{

    //! Return the chemical potentials of the species in the solution. Units: J/kmol.
    /*!
     * The chemical potential of species *k* is calculated as:
     * @f[
     * \begin{align}
     *  \mu_k(T_k, X_k, P)
     *      &= \mu_k^o(T_k, P) + R T_k \ln(X_k) \\
     *      &= \mu^0_k(T_k)(T_k)
     *         + R T_k \ln\left(\frac{P}{P^0}\right)
     *         + R T_k \ln(X_k) \\
     *      &= h^0_k(T_k)
     *         - T_k s^0_k(T_k)
     *         + R T_k \ln\left(\frac{P X_k}{P^0}\right) \\
     * \end{align}
     * @f]
     * Note that it is also equal to the partial molar Gibbs free energy.
     */
    void getChemPotentials(span<double> mu) const override;

    //! Return the partial molar enthalpies of the species in the solution. Units: J/kmol.
    /*!
     * The partial molar enthalpy of species *k* is
     * @f$ \tilde{h}_k(T_k,P) = h^o_k(T,P) = h^{ref}_k(T_k) @f$,
     * where heavy-species properties are evaluated at the gas temperature @f$ T @f$,
     * and electron properties are evaluated at the electron temperature @f$ T_e @f$.
     */
    void getPartialMolarEnthalpies(span<double> hbar) const override;

    //! Return the partial molar entropies of the species in the solution. Units: J/kmol/K.
    /*!
     * The partial molar enthalpy of species *k* is:
     * @f[
     * \tilde{s}_k(T_k, P) = s^0_k(T_k) - R \ln \frac{X_k P}{P^0}
     * @f]
     * where heavy-species properties are evaluated at the gas temperature,
     * and electron properties are evaluated at the electron temperature.
     * Here, @f$ P @f$ is the total pressure of the mixture, as defined in #pressure().
     */
    void getPartialMolarEntropies(span<double> sbar) const override;

    //! Return the partial molar internal energies of the species in the solution. Units: J/kmol.
    /*!
     * The partial molar internal energy of species *k* is calculated as:
     * @f[
     * \tilde{u}_k(T_k,P) = \tilde{h}_k(T_k,P) - R T_k = h^{ref}_k(T_k) - R T_k
     * @f]
     * where @f$ \tilde{h}_k @f$ is the partial molar enthalpy of species *k*, and
     * @f$ T_k @f$ is the temperature at which the partial molar enthalpy is evaluated.
     * For heavy species, the partial molar enthalpy is evaluated at the gas temperature,
     * while for electrons, it is evaluated at the electron temperature.
     */
    void getPartialMolarIntEnergies(span<double> ubar) const override;

    //! Return the partial molar volumes of the species in the solution. Units: m³/kmol.
    /*!
     * For a multitemperature system,
     * @f[
     *   v_k = \frac{R T_k}{P}
     * @f]
     * where @f$ T_k @f$ is the temperature at which the
     * partial molar enthalpy is evaluated.
     */
    void getPartialMolarVolumes(span<double> vbar) const override;

    //! @}
    //! @name  Properties of the Standard State of the Species in the Solution
    //! @{

    //! Return the standard chemical potentials of the species. Units: J/kmol.
    /*!
     * The standard chemical potentials (or standard state Gibbs free energy)
     * of species *k* is calculated as:
     * @f[
     * \begin{align}
     *  \mu^0_k(T_k, X_k, P)
     *      &= mu^0_k(T_k)(T_k)
     *         + R T_k \ln\left(\frac{P}{P^0}\right) \\
     *      &= h^0_k(T_k) - T_k s^0_k(T_k)
     *         + R T_k \ln\left(\frac{P X_k}{P^0}\right) \\
     * \end{align}
     * @f]
     * Here, @f$ P @f$ is the total pressure of the mixture, as defined in #pressure().
     */
    void getStandardChemPotentials(span<double> muStar) const override;

    //! Return the standard molar volumes of the species. Units: m³/kmol.
    /*!
     * For a multitemperature system,
     * @f[
     *   v_k = \frac{R T_k}{P},
     * @f]
     * where @f$ T_k @f$ is the temperature at which
     * the standard molar volume is evaluated.
     */
    void getStandardVolumes(span<double> vol) const override;

    //! @}
    //! @name Thermodynamic Values for the Species Reference States
    //! @{

    //! Return the reference chemical potentials of the species. Units: J/kmol.
    /*!
     * The reference chemical potentials (or reference state Gibbs free energy)
     * of species *k* is calculated as:
     * @f[
     *  \mu^0_k(T_k) = h^0_k(T_k) - T_k s^0_k(T_k)
     * @f]
     */
    void getGibbs_ref(span<double> g) const override;

    //! Return the molar volumes of the species reference states. Units: m³/kmol.
    /*!
     * The molar volumes of the species reference states of species *k*
     * is calculated as:
     * @f[
     *  v^o_k(T_k) = \frac{R T_k}{P^0}
     * @f]
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.

     */
    void getStandardVolumes_ref(span<double> vol) const override;

    //! @}
    //! @name Setting the State
    //!
    //! For a plasma phase, setting the state requires specifying both
    //! the heavy-species (gas) temperature and the electron temperature.
    //! @{

    //! Set the state using an AnyMap containing any combination of properties
    //! supported by the thermodynamic model
    /*!
     * Accepted keys are:
     * * `X` (mole fractions)
     * * `Y` (mass fractions)
     * * `T` or `Tg` or `gas-temperature` [K]
     * * `Te` or `electron-temperature` [K]
     * * `P` or `pressure` [Pa]
     * * `H` or `enthalpy` [J/kg]
     * * `U` or `internal-energy` [J/kg]
     * * `S` or `entropy` [J/kg/K]
     * * `V` or `specific-volume` [m^3/kg]
     * * `D` or `density` [kg/m^3]
     *
     * Composition can be specified as either an AnyMap of species names to
     * values or as a composition string. All other values can be given as
     * floating point values in Cantera's default units, or as strings with the
     * units specified, which will be converted using the Units class.
     *
     * If no thermodynamic property pair is given, or only one of temperature or
     * pressure is given, then 298.15 K and 101325 Pa will be used as necessary
     * to fully set the state.
     *
     * Set the electron temperature first, and call ThermoPhase::setState.
     */
    void setState(const AnyMap& state) override;


    //! @}

protected:
    //! Update the species reference state thermodynamic functions
    /*!
     *  This method is called each time a thermodynamic property is requested,
     *  to check whether the internal species properties within the object
     *  need to be updated (like getPartialMolarCp, getEnthalpy_RT, getEntropy_R,
     *  getGibbs_RT, getIntEnergy_RT, getCp_R, and getEnthalpy_RT_ref and alike).
     *  Currently, this updates the species thermo polynomial values for the current
     *  value of the gas temperature for heavy species, and of the electron
     *  temperature for the electron species. A check is made to see if gas
     *  or electron temperatures have changed since the last evaluation.
     *  This object does not contain any persistent data that depends on the
     *  concentration, that needs to be updated. The state object modifies its
     *  concentration dependent information at the time the setMoleFractions()
     *  (or equivalent) call is made.
     */
    void updateThermo() const override;

    //! When electron energy distribution changed, plasma properties such as
    //! electron-collision reaction rates need to be re-evaluated.
    void electronEnergyDistributionChanged();

    //! When electron energy level changed, plasma properties such as
    //! electron-collision reaction rates need to be re-evaluate.
    //! In addition, the cross-sections need to be interpolated at
    //! the new level.
    void electronEnergyLevelChanged();

    //! Check the electron energy levels
    /*!
     *  The values of electron energy levels need to be positive and
     *  monotonically increasing.
     */
    void checkElectronEnergyLevels() const;

    //! Check the electron energy distribution
    /*!
     *  This method check the electron energy distribution for the criteria
     *  below.
     *
     *  1. The values of electron energy distribution cannot be negative.
     *
     *  2. If the last value of electron energy distribution is larger
     *  than 0.01, it will raise a warning to suggest using a higher electron
     *  energy levels.
     */
    void checkElectronEnergyDistribution() const;

    //! Set isotropic electron energy distribution
    void setIsotropicElectronEnergyDistribution();

    //! Update electron temperature (K) From energy distribution.
    //! #m_electronTemp
    void updateElectronTemperatureFromEnergyDist();

    //! Electron energy distribution norm
    void normalizeElectronEnergyDistribution();

    //! Update interpolated cross section of a collision
    bool updateInterpolatedCrossSection(size_t k);

    //! Update electron energy distribution difference
    void updateElectronEnergyDistDifference();

    // Electron energy order in the exponential term
    double m_isotropicShapeFactor = 1.0;

    //! Number of points of electron energy levels
    size_t m_nPoints = 1001;

    //! electron energy levels [ev]. Length: #m_nPoints
    Eigen::ArrayXd m_electronEnergyLevels;

    //! Normalized electron energy distribution vector [-]
    //! Length: #m_nPoints
    Eigen::ArrayXd m_electronEnergyDist;

    //! Index of electron species.
    size_t m_electronSpeciesIndex = npos;

    //! Electron temperature [K].
    double m_electronTemp;

    //! Electron energy distribution type. Can be "isotropic", "discretized" or "Boltzmann-two-term".
    string m_distributionType = "isotropic";

    //! Numerical quadrature method for electron energy distribution.
    string m_quadratureMethod = "simpson";

    //! Flag of normalizing electron energy distribution.
    bool m_do_normalizeElectronEnergyDist = true;

    //! Indices of inelastic collisions in m_crossSections.
    vector<size_t> m_kInelastic;

    //! Indices of elastic collisions in m_crossSections.
    vector<size_t> m_kElastic;

    //! electric field [V/m].
    double m_electricField = 0.0;

    //! electric field freq [Hz].
    double m_electricFieldFrequency = 0.0;

    //! Cross section data. m_crossSections[i][j], where i is the specific process,
    //! j is the index of vector. Unit: [m^2]
    vector<vector<double>> m_crossSections;

    //! Electron energy levels corresponding to the cross section data. m_energyLevels[i][j],
    //! where i is the specific process, j is the index of vector. Unit: [eV]
    vector<vector<double>> m_energyLevels;

    //! ionization degree for the electron-electron collisions (tmp is the previous one)
    //double m_ionDegree = 0.0;

    //! Electron energy distribution Difference dF/dε (V^-5/2)
    Eigen::ArrayXd m_electronEnergyDistDiff;

    //! Elastic electron energy loss coefficients (eV m3/s)
    /*! The elastic electron energy loss coefficient for species k is,
     *   @f[
     *     K_k = \frac{2 m_e}{m_k} \sqrt{\frac{2 e}{m_e}} \int_0^{\infty} \sigma_k
     *           \epsilon^2 \left( F_0 + \frac{k_\text{B} T}{e}
     *           \frac{\partial F_0}{\partial \epsilon} \right) d \epsilon,
     *   @f]
     * where @f$ m_e @f$ [kg] is the electron mass, @f$ \epsilon @f$ [V] is the
     * electron energy, @f$ \sigma_k @f$ [m2] is the reaction collision cross section,
     * @f$ F_0 @f$ [V^(-3/2)] is the normalized electron energy distribution function.
     */
    vector<double> m_elasticElectronEnergyLossCoefficients;

    //! Updates the elastic electron energy loss coefficient for collision index i
    /*! Calculates the elastic energy loss coefficient using the current electron
        energy distribution and cross sections.
    */
    void updateElasticElectronEnergyLossCoefficient(size_t i);

    //! Update elastic electron energy loss coefficients
    /*! Used by elasticPowerLoss() and other plasma property calculations that
        depends on #m_elasticElectronEnergyLossCoefficients. This function calls
        updateInterpolatedCrossSection() before calling
        updateElasticElectronEnergyLossCoefficient()
    */
    void updateElasticElectronEnergyLossCoefficients();

private:

    //! Solver used to calculate the EEDF based on electron collision rates
    unique_ptr<EEDFTwoTermApproximation> m_eedfSolver = nullptr;

    //! Electron energy distribution change variable. Whenever
    //! #m_electronEnergyDist changes, this int is incremented.
    int m_distNum = -1;

    //! Electron energy level change variable. Whenever
    //! #m_electronEnergyLevels changes, this int is incremented.
    int m_levelNum = -1;

    //! The list of shared pointers of plasma collision reactions
    vector<shared_ptr<Reaction>> m_collisions;

    //! The list of shared pointers of collision rates
    vector<shared_ptr<ElectronCollisionPlasmaRate>> m_collisionRates;

    //! The collision-target species indices of #m_collisions
    vector<size_t> m_targetSpeciesIndices;

    //! The list of whether the interpolated cross sections is ready
    vector<bool> m_interp_cs_ready;

    //! Set collisions. This function sets the list of collisions and
    //! the list of target species using #addCollision.
    void setCollisions();

    //! Add a collision and record the target species
    void addCollision(shared_ptr<Reaction> collision);

    //! Saved electron temperature during an equilibrium solve
    double m_electronTempEquil = 0.0;

    //! Lock flag (default off)
    bool m_inEquilibrate = false;

    //! Work array
    mutable std::vector<double> m_work;
};

}

#endif
