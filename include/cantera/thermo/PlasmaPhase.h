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

    //! Set electron energy levels.
    //! @param  levels The vector of electron energy levels (eV).
    //!                Length: #m_nPoints.
    void setElectronEnergyLevels(span<const double> levels);

    //! Get electron energy levels.
    //! @param  levels The vector of electron energy levels (eV). Length: #m_nPoints
    void getElectronEnergyLevels(span<double> levels) const {
        Eigen::Map<Eigen::ArrayXd>(levels.data(), levels.size()) = m_electronEnergyLevels;
    }

    //! Set discretized electron energy distribution.
    //! @param  levels The vector of electron energy levels (eV).
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

    //! The shape factor of isotropic electron energy distribution
    double isotropicShapeFactor() const {
        return m_isotropicShapeFactor;
    }

    //! Set the internally stored electron temperature of the phase (K).
    //! @param  Te Electron temperature in Kelvin
    void setElectronTemperature(double Te) override;

    //! Set mean electron energy [eV]. This method also sets electron temperature
    //! accordingly.
    void setMeanElectronEnergy(double energy);

    //! Get electron energy distribution type
    string electronEnergyDistributionType() const {
        return m_distributionType;
    }

    //! Set electron energy distribution type
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

    //! Mean electron energy [eV]
    double meanElectronEnergy() const {
        return 3.0 / 2.0 * electronTemperature() * Boltzmann / ElectronCharge;
    }

    //! Set flag of automatically normalize electron energy distribution
    //! Flag: #m_do_normalizeElectronEnergyDist
    void enableNormalizeElectronEnergyDist(bool enable) {
        m_do_normalizeElectronEnergyDist = enable;
    }

    //! Flag of automatically normalize electron energy distribution.
    //! Flag: #m_do_normalizeElectronEnergyDist
    bool normalizeElectronEnergyDistEnabled() const {
        return m_do_normalizeElectronEnergyDist;
    }

    bool addSpecies(shared_ptr<Species> spec) override;

    //! Electron Temperature (K)
    //!     @return The electron temperature of the phase
    double electronTemperature() const override {
        return m_electronTemp;
    }

    //! Return the Gas Constant multiplied by the current electron temperature
    /*!
     *  The units are Joules kmol-1
     */
    double RTe() const {
        return electronTemperature() * GasConstant;
    }

    /**
     * Electron pressure. Units: Pa.
     * @f[P = n_{k_e} R T_e @f]
     */
    virtual double electronPressure() const {
        return GasConstant * concentration(m_electronSpeciesIndex) *
               electronTemperature();
    }

    //! Number of electron levels
    size_t nElectronEnergyLevels() const {
        return m_nPoints;
    }

    //! Number of electron collision cross sections
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

    //! Electron Species Index
    size_t electronSpeciesIndex() const {
        return m_electronSpeciesIndex;
    }

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

    //! Return the molar entropy. Units: J/kmol/K.
    /*!
    * For an ideal gas mixture with an additional electron species,
    * @f[
    * \hat s(T,T_e,p,X)
    * = \sum_{k\neq k_e} X_k\left[\hat s_k^0(T) - R\ln\left(\frac{X_k p}{p^0}\right)\right]
    * + X_{k_e}\left[\hat s_{k_e}^0(T_e) - R\ln\left(\frac{X_{k_e} p_e}{p^0}\right)\right],
    * @f]
    * where heavy-species properties are evaluated at @f$T@f$, electron properties at
    * @f$T_e@f$, and the electron mixing term uses the electron pressure
    * @f$p_e = n_{k_e} R T_e@f$.
    *
    * @see MultiSpeciesThermo
    */
    double entropy_mole() const override;

    //! Return the molar Gibbs free energy. Units: J/kmol.
    /*!
    * For an ideal gas mixture with an additional electron species,
    * @f[
    *   \hat g(T, T_e, p, X) = \sum_k X_k \mu_k(T_k, p_k, X),
    * @f]
    * where heavy species use the gas temperature @f$T@f$ and bulk pressure,
    * while the electron chemical potential uses @f$T_e@f$ and
    * @f$p_e = n_{k_e} R T_e@f$.
    *
    * @see MultiSpeciesThermo
    */
    double gibbs_mole() const override;

    //! Return the molar internal energy. Units: J/kmol.
    /*!
    * For an ideal gas mixture with an additional electron species,
    * @f[
    *   \hat u(T,T_e) = \sum_{k \neq k_e} X_k \hat u^0_k(T) + X_{k_e} \hat u^0_{k_e}(T_e),
    * @f]
    * where @f$\hat u^0_k(T) = \hat h^0_k(T) - R T@f$ for heavy species and
    * @f$\hat u^0_{k_e}(T_e) = \hat h^0_{k_e}(T_e) - R T_e@f$ for electrons.
    *
    * @see MultiSpeciesThermo
    */
    double intEnergy_mole() const override;

    void getEntropy_R(span<double> sr) const override;

    void getGibbs_RT(span<double> grt) const override;

    void getGibbs_ref(span<double> g) const override;

    void getStandardVolumes_ref(span<double> vol) const override;

    void getChemPotentials(span<double> mu) const override;

    void getStandardChemPotentials(span<double> muStar) const override;

    void getPartialMolarEnthalpies(span<double> hbar) const override;

    void getPartialMolarEntropies(span<double> sbar) const override;

    void getPartialMolarIntEnergies(span<double> ubar) const override;

    void getParameters(AnyMap& phaseNode) const override;

    void setParameters(const AnyMap& phaseNode,
                       const AnyMap& rootNode=AnyMap()) override;

    //! Update the electron energy distribution.
    void updateElectronEnergyDistribution();

    //! Electron species name
    string electronSpeciesName() const {
        return speciesName(m_electronSpeciesIndex);
    }

    //! Return the distribution Number #m_distNum
    int distributionNumber() const {
        return m_distNum;
    }

    //! Return the electron energy level Number #m_levelNum
    int levelNumber() const {
        return m_levelNum;
    }

    //! Get the indicies for inelastic electron collisions
    //! @since New in %Cantera 3.2.
    const vector<size_t>& kInelastic() const {
        return m_kInelastic;
    }

    //! Get the indices for elastic electron collisions
    //! @since New in %Cantera 3.2.
    const vector<size_t>& kElastic() const {
        return m_kElastic;
    }

    //! target of a specific process
    //! @since New in %Cantera 3.2.
    size_t targetIndex(size_t i) const {
        return m_targetSpeciesIndices[i];
    }

    //! Get the frequency of the applied electric field [Hz]
    //! @since New in %Cantera 3.2.
    double electricFieldFrequency() const {
        return m_electricFieldFrequency;
    }

    //! Get the applied electric field strength [V/m]
    double electricField() const {
        return m_electricField;
    }

    //! Set the absolute electric field strength [V/m]
    void setElectricField(double E) {
        m_electricField = E;
    }

    //! Calculate the degree of ionization
    //double ionDegree() const {
    //    double ne = concentration(m_electronSpeciesIndex); // [kmol/m³]
    //    double n_total = molarDensity();                   // [kmol/m³]
    //    return ne / n_total;
    //}

    //! Get the reduced electric field strength [V·m²]
    double reducedElectricField() const {
        return m_electricField / (molarDensity() * Avogadro);
    }

    //! Set reduced electric field given in [V·m²]
    void setReducedElectricField(double EN) {
        m_electricField = EN * molarDensity() * Avogadro; // [V/m]
    }

    virtual void setSolution(std::weak_ptr<Solution> soln) override;

    /**
     * The elastic power loss (J/s/m³)
     *   @f[
     *     P_k = N_A N_A C_e e \sum_k C_k K_k,
     *   @f]
     * where @f$ C_k @f$ and @f$ C_e @f$ are the concentration (kmol/m³) of the
     * target species and electrons, respectively. @f$ K_k @f$ is the elastic
     * electron energy loss coefficient (eV-m³/s).
     */
    double elasticPowerLoss();

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

    //! Saved electron temperature during an equilibrium solve
    double m_electronTempEquil = 0.0;

    //! Lock flag (default off)
    bool m_inEquilibrate = false;

    //! Work array
    mutable std::vector<double> m_work;
};

}

#endif
