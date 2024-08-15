/**
 * @file PlasmaPhase.h
 * Header file for class PlasmaPhase.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_PLASMAPHASE_H
#define CT_PLASMAPHASE_H

#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

/**
 * Base class for a phase with plasma properties. This class manages the
 * plasma properties such as electron energy distribution function (EEDF).
 * There are two ways to define the electron distribution and electron
 * temperature. The first method uses setElectronTemperature() to set
 * the electron temperature which is used to calculate the electron energy
 * distribution with isotropic-velocity model. The generalized electron
 * energy distribution for isotropic-velocity distribution can be
 * expressed as [1,2],
 *   @f[
 *          f(\epsilon) = c_1 \frac{\sqrt{\epsilon}}{\epsilon_m^{3/2}}
 *          \exp(-c_2 (\frac{\epsilon}{\epsilon_m})^x),
 *   @f]
 * where @f$ x = 1 @f$ and @f$ x = 2 @f$ correspond to the Maxwellian and
 * Druyvesteyn (default) electron energy distribution, respectively.
 * @f$ \epsilon_m = 3/2 T_e @f$ [eV] (mean electron energy). The second
 * method uses setDiscretizedElectronEnergyDist() to manually set electron
 * energy distribution and calculate electron temperature from mean electron
 * energy, which is calculated as [3],
 *   @f[
 *          \epsilon_m = \int_0^{\infty} \epsilon^{3/2} f(\epsilon) d\epsilon,
 *   @f]
 * which can be calculated using trapezoidal rule,
 *   @f[
 *          \epsilon_m = \sum_i (\epsilon^{5/2}_{i+1} - \epsilon^{5/2}_i)
 *                       (f(\epsilon_{i+1}) + f(\epsilon_i)) / 2,
 *   @f]
 * where @f$ i @f$ is the index of energy levels.
 *
 * For references, see Gudmundsson @cite gudmundsson2001; Khalilpour and Foroutan
 * @cite khalilpour2020; Hagelaar and Pitchford @cite hagelaar2005, and BOLOS
 * @cite BOLOS.
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

    string type() const override {
        return "plasma";
    }

    void initThermo() override;

    //! Set electron energy levels.
    //! @param  levels The vector of electron energy levels (eV).
    //!                Length: #m_nPoints.
    //! @param  length The length of the @c levels.
    void setElectronEnergyLevels(const double* levels, size_t length);

    //! Get electron energy levels.
    //! @param  levels The vector of electron energy levels (eV). Length: #m_nPoints
    void getElectronEnergyLevels(double* levels) const {
        Eigen::Map<Eigen::ArrayXd>(levels, m_nPoints) = m_electronEnergyLevels;
    }

    //! Set discretized electron energy distribution.
    //! @param  levels The vector of electron energy levels (eV).
    //!                Length: #m_nPoints.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    //! @param  length The length of the vectors, which equals #m_nPoints.
    void setDiscretizedElectronEnergyDist(const double* levels,
                                          const double* distrb,
                                          size_t length);

    //! Set discretized electron energy distribution.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    //! @param  length The length of the vectors, which equals #m_nPoints.
    void setDiscretizedElectronEnergyDist(const double* distrb,
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
     * \see MultiSpeciesThermo
     */
    double enthalpy_mole() const override;

    double cp_mole() const override {
        throw NotImplementedError("PlasmaPhase::cp_mole");
    }

    double entropy_mole() const override {
        throw NotImplementedError("PlasmaPhase::entropy_mole");
    }

    double gibbs_mole() const override {
        throw NotImplementedError("PlasmaPhase::gibbs_mole");
    }

    double intEnergy_mole() const override {
        throw NotImplementedError("PlasmaPhase::intEnergy_mole");
    }

    void getEntropy_R(double* sr) const override;

    void getGibbs_RT(double* grt) const override;

    void getGibbs_ref(double* g) const override;

    void getStandardVolumes_ref(double* vol) const override;

    void getChemPotentials(double* mu) const override;

    void getStandardChemPotentials(double* muStar) const override;

    void getPartialMolarEnthalpies(double* hbar) const override;

    void getPartialMolarEntropies(double* sbar) const override;

    void getPartialMolarIntEnergies(double* ubar) const override;

    void getParameters(AnyMap& phaseNode) const override;

    void setParameters(const AnyMap& phaseNode,
                       const AnyMap& rootNode=AnyMap()) override;

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

    //! Update electron energy distribution.
    void updateElectronEnergyDistribution();

    //! Set isotropic electron energy distribution
    void setIsotropicElectronEnergyDistribution();

    //! Update electron temperature (K) From energy distribution.
    //! #m_electronTemp
    void updateElectronTemperatureFromEnergyDist();

    //! Electron energy distribution norm
    void normalizeElectronEnergyDistribution();

    // Electron energy order in the exponential term
    double m_isotropicShapeFactor = 2.0;

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

private:
    //! Electron energy distribution change variable. Whenever
    //! #m_electronEnergyDist changes, this int is incremented.
    int m_distNum = -1;

    //! Electron energy level change variable. Whenever
    //! #m_electronEnergyLevels changes, this int is incremented.
    int m_levelNum = -1;
};

}

#endif
