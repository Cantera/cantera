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
 *   \f[
 *          f(\epsilon) = c_1 \frac{\sqrt{\epsilon}}{\epsilon_m^{3/2}}
 *          \exp(-c_2 (\frac{\epsilon}{\epsilon_m})^x),
 *   \f]
 * where \f$ x = 1 \f$ and \f$ x = 2 \f$ correspond to the Maxwellian and
 * Druyvesteyn (default) electron energy distribution, respectively.
 * \f$ \epsilon_m = 3/2 T_e \f$ [eV] (mean electron energy). The second
 * method uses setDiscretizedElectronEnergyDist() to manually set electron
 * energy distribution and calculate electron temperature from mean electron
 * energy, which is calculated as [3],
 *   \f[
 *          \epsilon_m = \int_0^{\infty} \epsilon^{3/2} f(\epsilon) d\epsilon,
 *   \f]
 * which can be calculated using trapezoidal rule,
 *   \f[
 *          \epsilon_m = \sum_i (\epsilon^{5/2}_{i+1} - \epsilon^{5/2}_i)
 *                       (f(\epsilon_{i+1}) + f(\epsilon_i)) / 2,
 *   \f]
 * where \f$ i \f$ is the index of energy levels.
 *
 * References:
 *
 * [1] J. T. Gudmundsson. On the effect of the electron energy distribution on the
 * plasma parameters of an argon discharge: a global (volume-averaged) model study.
 * Plasma Sources Science and Technology, 10.1 (2001): 76.
 * doi: https://doi.org/10.1088/0963-0252/10/1/310
 *
 * [2] H. Khalilpour and G. Foroutan. The effects of electron energy distribution
 * function on the plasma sheath structure in the presence of charged nanoparticles
 * Journal of Plasma Physics 86.2 (2020).
 * doi: https://doi.org/10.1017/S0022377820000161
 *
 * [3] G. J. M. Hagelaar and L. C. Pitchford
 * "Solving the Boltzmann equation to obtain electron transport
 * coefficients and rate coefficients for fluid models."
 * Plasma Sources Science and Technology 14.4 (2005): 722.
 * doi: https://doi.org/10.1088/0963-0252/14/4/011
 *
 * [4] A. Luque, "BOLOS: An open source solver for the Boltzmann equation,"
 * https://github.com/aluque/bolos.
 *
 * @warning  This class is an experimental part of %Cantera and may be
 *           changed or removed without notice.
 * @todo Implement electron Boltzmann equation solver to solve EEDF.
 *       https://github.com/Cantera/enhancements/issues/127
 * @ingroup phase
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
    explicit PlasmaPhase(const std::string& inputFile="",
                         const std::string& id="");

    virtual std::string type() const {
        return "plasma";
    }

    virtual void initThermo();

    //! Set electron energy levels.
    //! @param  levels The vector of electron energy levels (eV).
    //!                Length: #m_nPoints.
    //! @param  length The length of the \p levels.
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

    //! Get electron energy distribution.
    //! @param  distrb The vector of electron energy distribution.
    //!                Length: #m_nPoints.
    void getElectronEnergyDistribution(double* distrb) const {
        Eigen::Map<Eigen::ArrayXd>(distrb, m_nPoints) = m_electronEnergyDist;
    }

    //! Set the shape factor of isotropic electron energy distribution.
    //! Note that \f$ x = 1 \f$ and \f$ x = 2 \f$ correspond to the
    //! Maxwellian and Druyvesteyn distribution, respectively.
    //! @param  x The shape factor
    void setIsotropicShapeFactor(double x);

    //! The shape factor of isotropic electron energy distribution
    double isotropicShapeFactor() const {
        return m_isotropicShapeFactor;
    }

    //! Set the internally stored electron temperature of the phase (K).
    //! @param  Te Electron temperature in Kelvin
    virtual void setElectronTemperature(double Te);

    //! Set mean electron energy [eV]. This method also sets electron temperature
    //! accordingly.
    void setMeanElectronEnergy(double energy);

    //! Get electron energy distribution type
    std::string electronEnergyDistributionType() const {
        return m_distributionType;
    }

    //! Set electron energy distribution type
    void setElectronEnergyDistributionType(const std::string& type);

    //! Numerical quadrature method. Method: #m_quadratureMethod
    std::string quadratureMethod() const {
        return m_quadratureMethod;
    }

    //! Set numerical quadrature method for intergating electron
    //! energy distribution function. Method: #m_quadratureMethod
    void setQuadratureMethod(const std::string& method) {
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

    virtual bool addSpecies(shared_ptr<Species> spec);

    //! Electron Temperature (K)
    //!     @return The electron temperature of the phase
    virtual double electronTemperature() const {
        return m_electronTemp;
    }

    //! Number of electron levels
    size_t nElectronEnergyLevels() const {
        return m_nPoints;
    }

    virtual void getParameters(AnyMap& phaseNode) const;

    virtual void setParameters(const AnyMap& phaseNode,
                               const AnyMap& rootNode=AnyMap());

protected:
    virtual void updateThermo() const;

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
    double m_isotropicShapeFactor;

    //! Number of points of electron energy levels
    size_t m_nPoints;

    //! electron energy levels [ev]. Length: #m_nPoints
    Eigen::ArrayXd m_electronEnergyLevels;

    //! Normalized electron energy distribution vector [-]
    //! Length: #m_nPoints
    Eigen::ArrayXd m_electronEnergyDist;

    //! Index of electron species
    size_t m_electronSpeciesIndex;

    //! Electron temperature [K]
    double m_electronTemp;

    //! Electron energy distribution type
    std::string m_distributionType;

    //! Numerical quadrature method for electron energy distribution
    std::string m_quadratureMethod;

    //! Flag of normalizing electron energy distribution
    bool m_do_normalizeElectronEnergyDist;
};

}

#endif
