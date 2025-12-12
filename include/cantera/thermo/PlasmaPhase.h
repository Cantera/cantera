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
 * @warning  This class is an experimental part of %Cantera and may be
 *           changed or removed without notice.
 * @todo Implement electron Boltzmann equation solver to solve EEDF.
 *       https://github.com/Cantera/enhancements/issues/127
 * @ingroup thermoprops
 */
class PlasmaPhase: public IdealGasPhase
{
public:
    //! Construct and initialize a PlasmaPhase object
    //! directly from an input file. The constructor initializes the electron
    //! energy distribution to be Druyvesteyn distribution (m_x = 2.0). The initial
    //! electron energy grid is set to a linear space which starts at 0.01 eV and ends
    //! at 1 eV with 1000 points.
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

    // ================================================================= //
    // ================================================================= //
    //! @name Overridden from IdealGasPhase or ThermoPhase
    //! @{
    bool addSpecies(shared_ptr<Species> spec) override;
    virtual void setSolution(std::weak_ptr<Solution> soln) override;
    void getParameters(AnyMap& phaseNode) const override;
    void setParameters(const AnyMap& phaseNode,
                       const AnyMap& rootNode=AnyMap()) override;
    //! @}
    // ================================================================= //
    // ================================================================= //
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
    // ================================================================= //
    // ================================================================= //
    //! @name Electron Energy Distribution Functions
    //! @{

    //! Set electron energy levels.
    //! @param  levels The vector of electron energy levels [eV].
    //!                Length: #m_nPoints.
    //! @param  length The length of the @c levels.
    void setElectronEnergyLevels(const double* levels, size_t length);

    //! Get electron energy levels.
    //! @param  levels The vector of electron energy levels [eV]. Length: #m_nPoints
    void getElectronEnergyLevels(double* levels) const {
        Eigen::Map<Eigen::ArrayXd>(levels, m_nPoints) = m_electronEnergyLevels;
    }

    //! Set discretized electron energy distribution.
    //! @param  levels The vector of electron energy levels [eV].
    //!                Length: #m_nPoints.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    //! @param  length The length of the vectors, which equals #m_nPoints.
    void setDiscretizedElectronEnergyDist(const double* levels,
                                          const double* distrb,
                                          size_t length);

    //! Get electron energy distribution.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    void getElectronEnergyDistribution(double* distrb) const {
        Eigen::Map<Eigen::ArrayXd>(distrb, m_nPoints) = m_electronEnergyDist;
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
     * The elastic power loss [J/s/m³]
     *   @f[
     *     P_k = N_A N_A C_e e \sum_k C_k K_k,
     *   @f]
     * where @f$ C_k @f$ and @f$ C_e @f$ are the concentration (kmol/m³) of the
     * target species and electrons, respectively. @f$ K_k @f$ is the elastic
     * electron energy loss coefficient (eV-m³/s).
     */
    double elasticPowerLoss();

    //! @}
    // ================================================================= //
    // ================================================================= //
    //! @name Molar Thermodynamic Properties of the Solution
    //! @{

    //! Return the Molar enthalpy. Units: J/kmol.
    /*!
     * For an ideal gas mixture with additional electron,
     * @f[
     * \hat h(T) = \sum_{k \neq k_e} X_k \hat h^0_k(T) + X_{k_e} \hat h^0_{k_e}(T_e),
     * @f]
     * and is a function only of temperature. The standard-state pure-species
     * enthalpies @f$ \hat h^0_k(T) @f$ are computed by the species
     * thermodynamic property manager.
     *
     * @see MultiSpeciesThermo
     */
    double enthalpy_mole() const override;

    double entropy_mole() const override {
        throw NotImplementedError("PlasmaPhase::entropy_mole");
    }

    double gibbs_mole() const override {
        throw NotImplementedError("PlasmaPhase::gibbs_mole");
    }

    double intEnergy_mole() const override {
        throw NotImplementedError("PlasmaPhase::intEnergy_mole");
    }

    double cp_mole() const override;
    // double cp_mass() const; // Already defined in ThermoPhase
    // double cv_mole() const; // Already defined in IdealGasPhase

    //! Return the constant-pressure specific heat capacity of electrons [J/kmol/K].
    /*!
     * For electrons, the specific heat capacity at constant pressure does not
     * depend on temperature. It is constant, and given by:
     * @f[
     * \hat c_{p,e} = \frac{5}{2} R
     * @f]
     * This comes from statistical mechanics, where only translational degrees of
     * freedom contribute to the specific heat capacity (no electronic, vibrational, or
     * rotational contributions). The factor of 5/2 arises from the 3/2R contribution
     * from translational motion, plus an additional R from the pressure-volume work
     * done at constant pressure.
     */
    double cp_mole_e() const;
    //! Return the constant-pressure specific heat capacity of electrons [J/kg/K].
    double cp_mass_e() const {
        return cp_mole_e() / m_Me;
    };
    //! Return the constant-volume specific heat capacity of electrons [J/kmol/K].
    double cv_mole_e() const {
        return cp_mole_e() - GasConstant;
    };
    //! Return the constant-volume specific heat capacity of electrons [J/kg/K].
    double cv_mass_e() const {
        return cp_mass_e() - GasConstant / m_Me;
    };

    //! Return the constant-pressure specific heat capacity of heavy species [J/kmol/K].
    /*!
     * The mean heavy-species specific heat capacity at constant pressure is
     * calculated as:
     * @f[
     * \hat c_{p,h} = \sum_{k \neq k_e} X_k \hat c^0_{p,k}(T_g)
     * @f]
     * where @f$ \hat c^0_{p,k}(T_g) @f$ is the standard-state specific heat capacity
     * at constant pressure for species @f$ k @f$ and @f$ T_g @f$ is the gas temperature.
     */
    double cp_mole_h() const;
    //! Return the constant-pressure specific heat capacity of heavy species [J/kg/K].
    double cp_mass_h() const;
    //! Return the constant-volume specific heat capacity of heavy species [J/kmol/K].
    double cv_mole_h() const;
    //! Return the constant-volume specific heat capacity of heavy species [J/kg/K].
    double cv_mass_h() const;

    //! @}
    // ================================================================= //
    // ================================================================= //
    //! @name Mechanical Equation of State
    //! @{

    //! Return the mean temperature of the plasma phase [K].
    /*!
     * In a plasma phase, the electron temperature can differ from the
     * heavy-species (gas) temperature. Therefore, the mean temperature is
     * defined as a mole-fraction-weighted average of the electron and
     * heavy-species temperatures:
     * @f[
     *      \overline{T} = \sum_{k \neq k_e} X_k T_g + X_{k_e} T_e
     *             = (1 - X_{k_e}) T_g + X_{k_e} T_e
     *             = T_g + X_{k_e} (T_e - T_g)
     * @f]
     * See the `pressure()` method for usage of the mean temperature in
     * calculating the pressure of the plasma phase.
     */
    double meanTemperature() const;

    //! Return the pressure of the plasma phase [Pa].
    /*!
     * The pressure of the plasma phase is calculated using the mean temperature,
     * which is a mole-fraction-weighted average of the electron and heavy-species
     * temperatures.
     * @f[
     *      P = \sum_k n_k k_B T_k
     *        = \sum_{k \neq e} n_k k_B T_g + n_{e} k_B T_e
     *        = (n_{total} - n_{e}) k_B T_g + n_{e} k_B T_e
     *        = n_{total} (1 - X_e) k_B T_g + n_{total} X_e k_B T_e
     *        = n_{total} k_B \overline{T}
     * @f]
     * where @f$ \overline{T} @f$ is the mean temperature of the plasma phase,
     * defined in the `meanTemperature()` method.
     * Here, @f$ n_k @f$ is the number density of species @f$ k @f$ [in 1/m³],
     * @f$ n_{total} @f$ is the total number density of the phase [in 1/m³],
     * @f$ X_e @f$ is the mole fraction of electrons, @f$ T_g @f$ is the
     * heavy-species (gas) temperature [K], @f$ T_e @f$ is the electron temperature [K],
     * and @f$ k_B @f$ is the Boltzmann constant [J/K].
     * 
     * The number density times Boltzmann constant can be expressed as
     * @f$ n_{total} k_B = C_{total} R @f$, where @f$ C_{total} @f$ is the
     * total molar concentration of the phase [in kmol/m³] and
     * @f$ R @f$ is the gas constant [J/kmol/K], so that:
     * @f[
     *     P = C_{total} R \overline{T}.
     * @f]
     */
    double pressure() const override;

    //! Set the pressure at constant temperature and composition.
    /*!
     * Units: Pa.
     * This method is implemented by setting the mass density to
     * @f[
     * \rho = \frac{P \overline{W}}{\hat R \overline{T} }.
     * @f]
     *
     * @param p Pressure (Pa)
     */
    void setPressure(double p) override {
        setDensity(p * meanMolecularWeight() / (GasConstant * meanTemperature()));
    }

    //! Return the Gas Constant multiplied by the current electron temperature [J/kmol].
    double RTe() const {
        return electronTemperature() * GasConstant;
    }

    //! Return the electron pressure [Pa]. 
    /*
     * @f[P = n_{k_e} R T_e @f]
     */
    virtual double electronPressure() const {
        return GasConstant * concentration(m_electronSpeciesIndex) *
               electronTemperature();
    }

    //! @}
    // ================================================================= //
    // ================================================================= //
    //! @name Partial Molar Properties of the Solution
    //! @{

    void getChemPotentials(double* mu) const override;
    void getPartialMolarEnthalpies(double* hbar) const override;
    // void getPartialMolarEntropies(double* sbar) const override;
    void getPartialMolarIntEnergies(double* ubar) const override;
    // void getPartialMolarCp(double* cpbar) const override;
    // void getPartialMolarVolumes(double* vbar) const override;

    //! @}
    // ================================================================= //
    // ================================================================= //
    //! @name  Properties of the Standard State of the Species in the Solution
    //! @{

    //! Return the standard chemical potentials of the species [J/kmol].
    /*!
     * For heavy species, this is identical to the IdealGasPhase
     * implementation. For electrons, the standard chemical potential
     * is evaluated at the electron temperature:
     * @f[
     *  \mu^0_{e}(T_e) = g_k^0(T_e) + RT_e \ln \left(\frac{P}{P^0}\right).
     * @f]
     */
    void getStandardChemPotentials(double* muStar) const override;

    // void getEnthalpy_RT(double* hrt) const override;
    // void getEntropy_R(double* sr) const override;
    // void getGibbs_RT(double* grt) const override;
    // void getIntEnergy_RT(double* urt) const override;
    // void getCp_R(double* cpr) const override;
    // void getStandardVolumes(double* vol) const override;

    //! @}
    // ================================================================= //
    // ================================================================= //
    //! @name Thermodynamic Values for the Species Reference States
    //! @{

    // void getEnthalpy_RT_ref(double* hrt) const override;
    // void getGibbs_RT_ref(double* grt) const override;

    //! Return the reference-state Gibbs free energys of the species [J/kmol].
    /*!
     * For heavy species, this is identical to the IdealGasPhase
     * implementation. For electrons, the reference-state Gibbs free energy
     * is evaluated at the electron temperature:
     * @f[
     *  \hat{g}^0_{e}(T_e) = \hat{h}^0_{e}(T_e) - T_e \hat{s}^0_{e}(T_e).
     * @f]
     */
    void getGibbs_ref(double* g) const override;

    // void getEntropy_R_ref(double* er) const override;
    // void getIntEnergy_RT_ref(double* urt) const override;
    // void getCp_R_ref(double* cprt) const override;
    void getStandardVolumes_ref(double* vol) const override;

    //! @}



protected:
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

    //! Index of electron species
    size_t m_electronSpeciesIndex = npos;

    //! Electron temperature [K]
    double m_electronTemp;

    //! Electron energy distribution type
    string m_distributionType = "isotropic";

    //! Numerical quadrature method for electron energy distribution
    string m_quadratureMethod = "simpson";

    //! Flag of normalizing electron energy distribution
    bool m_do_normalizeElectronEnergyDist = true;

    //! Indices of inelastic collisions in m_crossSections
    vector<size_t> m_kInelastic;

    //! Indices of elastic collisions in m_crossSections
    vector<size_t> m_kElastic;

    //! electric field [V/m]
    double m_electricField = 0.0;

    //! electric field freq [Hz]
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
     *           \epsilon^2 \left( F_0 + \frac{k_B T}{e}
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

    // Electron molar mass [kg/kmol].
    const double m_Me = ElectronMass * Avogadro;
};

}

#endif
