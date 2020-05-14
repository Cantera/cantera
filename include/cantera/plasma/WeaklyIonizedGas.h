/**
 * @file WeaklyIonizedGas.h
 * Header file for class WeaklyIonizedGas.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_WEAKLYIONIZEDGAS_H
#define CT_WEAKLYIONIZEDGAS_H

#include "cantera/plasma/PlasmaPhase.h"
#include <Eigen/Sparse>

namespace Cantera
{
typedef Eigen::SparseMatrix<double> SparseMat;

/**
 * This class calculates the electron energy distribution function (EEDF) in a weakly
 * ionized gas by modeling collisions between electrons and other species represented
 * by the class ElectronCrossSection. EEDF is used to calculate reaction rate coefficient
 * for plasma reaction and electron temperature for electron-temperature reaction
 * used in kinetics, and diffusivity/mobility of electron in transport.
 * Only electron-neutral collisions are considered for calculating the
 * electron energy distribution function (EEDF). Equation of EEDF becomes,
 * \f[
 *   \frac{d}{d \epsilon}\left(\tilde{W} F_0 - \tilde{D} \frac{d F_0}{d \epsilon}\right)
 *   = \tilde{S}
 * \f]
 * where
 * \f[
 *     \tilde{W} = -\gamma\epsilon^2\sigma_{\epsilon},
 * \f]
 * \f[
 *     \tilde{D} = \frac{\gamma}{3} \left(\frac{E}{N} \right)^2 \frac{\epsilon}{\tilde{\sigma}_m} +
 *                 \frac{\gamma k_B T}{e} \epsilon^2 \sigma_{\epsilon},
 * \f]
 * \f[
 *     \tilde{S} = \sum_{k=inelastic} \tilde{C}_{0,k} + G
 * \f]
 * where \f$ \gamma = (\frac{2 e}{m})^{1/2} \f$, \f$ \epsilon \f$ is the electron energy,
 * \f$ \sigma_{\epsilon} \f$ is the total elastic collision cross section,
 * \f$ E \f$ is electric field strength, \f$ N \f$ is gas number density,
 * \f$ \tilde{\sigma}_m \f$ is the effective total momentum-transfer cross section,
 * \f$ k_B \f$ is Boltzmann constant, \f$ e \f$ is the elementary charge,
 * \f$ \tilde{C}_{0,k} \f$ represents the rate of change in EEDF due to collisions,
 * and \f$ F_0 \f$ is the normalized EEDF.
 * 
 * The collision terms are,
 * \f[
 *     \tilde{C}_{0,k=elastic} = \gamma X_k \frac{2 m}{M_k} \frac{d}{d \epsilon}
 *         \left[\epsilon^2 \sigma_k \left( F_0 + \frac{k_B T}{e} \frac{d F_0}{\epsilon} \right) \right],
 * \f]
 * \f[
 *     \tilde{C}_{0,k=inelastic} = -\gamma X_k \epsilon \sigma_k F_0
 *                                 \big|^{\epsilon=\epsilon}_{\epsilon=\epsilon + u_k}
 * \f]
 * \f[
 *     \tilde{C}_{0,k=ionization} = -\gamma X_k \epsilon \sigma_k F_0
 *                                  \big|^{\epsilon=\epsilon}_{\epsilon=2\epsilon + u_k}
 * \f]
 * \f[
 *     \tilde{C}_{0,k=attachment} = -\gamma X_k \epsilon \sigma_k F_0
 * \f]
 * For exponential temporal growth model,
 * \f[
 *     G = \left[ \int_0^\infty \left(\sum_{k=ionization} X_k \sigma_k - \sum_{k=attachment} X_k \sigma_k \right)
 *         \epsilon F_0 d \epsilon \right] \epsilon^{1/2} F_0
 * \f]
 * where \f$ X_k \f$ is the mole fraction of the target species, and \f$ \sigma_k \f$ is the cross section.
 *
 * Reference:
 * [1] G. J. M. Hagelaar and L. C. Pitchford
 * "Solving the Boltzmann equation to obtain electron transport
 * coefficients and rate coefficients for fluid models."
 * Plasma Sources Science and Technology 14.4 (2005): 722.
 * doi: https://doi.org/10.1088/0963-0252/14/4/011
 * [2] A. Luque, "BOLOS: An open source solver for the Boltzmann equation,"
 * https://github.com/aluque/bolos.
 * [3] D. McElroy, C. Walsh, A. Markwick, M. Cordiner, K. Smith, T. Millar,
 * "The umist database for astrochemistry 2012,"
 * Astronomy & Astrophysics 550 (2013) A36.
 * doi: https://doi.org/10.1051/0004-6361/201220465
 * @ingroup plasma
 */
class WeaklyIonizedGas: public PlasmaPhase
{
public:
    WeaklyIonizedGas();

    virtual std::string type() const {
        return "WeaklyIonizedGas";
    }

    virtual double electronDiffusivity();
    virtual double electronMobility();
    virtual double meanElectronEnergy();
    virtual double powerGain();
    virtual double elasticPowerLoss();
    virtual double inelasticPowerLoss();
    virtual double totalCollisionFreq();
    virtual double reverseRateCoefficient(size_t k);
    virtual double rateCoefficient(size_t k);

    //! The real part of the mobility. This is used in power gain for case of AC.
    double realMobility();

    //! electron temperature
    //! If the reduced electric field is set, electron temperature is calculated
    //! from EEDF.
    virtual double electronTemperature();

protected:
    //! Calculate distribution function by solving Boltzmann equation
    //! with two-term approximate method.
    void calculateDistributionFunction();

    //! Calculate total cross section. The total cross section is defined as the summation
    //! of weighted (by mole fractions) cross sections of all coliision processes.
    void calculateTotalCrossSection();

    //! Calculate total elastic cross section \f$ \sigma_{\epsilon} \f$.
    /**
     * \f[
     *     \sigma_{\epsilon} = \sum_{k=elastic} \frac{2 m}{M_k} X_k \sigma_k,
     * \f]
     * where \f$ m \f$ is the mass of electron, \f$ M_k \f$ is the mass of the target species,
     * \f$ X_k \f$ is the mole fraction of the target species, \f$ \sigma_k \f$ is the cross section.
     */
    void calculateTotalElasticCrossSection();

    //! The integral in [a, b] of \f$x u(x) \exp[g (x_0 - x)]\f$
    //! assuming that u is linear with u(a) = u0 and u(b) = u1
    double integralPQ(double a, double b, double u0, double u1,
                       double g, double x0);

    //! Vector g is used by matrix_PQ.
    /**
     * \f[
     * g_i = \frac{1}{\epsilon_{i+1} - \epsilon_{i-1}} \ln(\frac{F_{0, i+1}}{F_{0, i-1}})
     * \f]
     */
    vector_fp vector_g(Eigen::VectorXd& f0);

    //! The matrix of scattering-out.
    /**
     * \f[
     * P_i = \sum_{k=inelastic} \gamma X_k \int_{\epsilon_i - 1/2}^{\epsilon_i + 1/2}
     * \epsilon \sigma_k exp[(\epsilon_i - \epsilon)g_i] d \epsilon
     * \f]
     */
    SparseMat matrix_P(vector_fp& g, size_t k);

    //! The matrix of scattering-in
    /**
     * \f[
     * Q_{i,j} = \sum_{k=inelastic} \gamma X_k \int_{\epsilon_1}^{\epsilon_2}
     * \epsilon \sigma_k exp[(\epsilon_j - \epsilon)g_j] d \epsilon
     * \f]
     */
    //! where the interval \f$[\epsilon_1, \epsilon_2]\f$ is the overlap of cell j,
    //! and cell i shifted by the threshold energy u_k:
    /**
     * \f[
     * \epsilon_1 = \min(\max(\epsilon_{i-1/2}+u_k, \epsilon_{j-1/2}),\epsilon_{j+1/2}),
     * \f]
     * \f[
     * \epsilon_2 = \min(\max(\epsilon_{i+1/2}+u_k, \epsilon_{j-1/2}),\epsilon_{j+1/2})
     * \f]
     */
    SparseMat matrix_Q(vector_fp& g, size_t k);

    //! Matrix A (Ax = b) of the equation of EEDF, which is discretized by the exponential scheme
    //! of Scharfetter and Gummel,
    /**
     * \f[
     *     \left[ \tilde{W} F_0 - \tilde{D} \frac{d F_0}{\epsilon} \right]_{i+1/2} =
     *     \frac{\tilde{W}_{i+1/2} F_{0,i}}{1 - \exp[-z_{i+1/2}]} +
     *     \frac{\tilde{W}_{i+1/2} F_{0,i+1}}{1 - \exp[z_{i+1/2}]}
     * \f]
     * where \f$ z_{i+1/2} = \tilde{w}_{i+1/2} / \tilde{D}_{i+1/2} \f$ (Peclet number).
     */
    SparseMat matrix_A(Eigen::VectorXd& f0);

    //! An iteration of solving electron energy distribution function
    Eigen::VectorXd iterate(Eigen::VectorXd& f0, double delta = 1e14);

    //! Iterate until convergence and obtain EEDF
    Eigen::VectorXd converge(Eigen::VectorXd& f0);

    //! Reduced net production frequency. Equation (10) of ref. [1]
    //! divided by N.
    //! @param f0 EEDF
    double netProductionFreq(Eigen::VectorXd& f0);

    //! electron temperature. For internal use only. This function is used to evaluate EEDF by
    //! comparing the resulting electron temeprature to gas temperature.
    double electronTemperature(Eigen::VectorXd f0);

    //! bi-Maxwellian Boltzmann factor. Assume that the excitation
    //! temperature equals to the gas temperature.
    double biMaxwellFraction(size_t k);

    //! Total electron cross section on the cell center of energy grid
    vector_fp m_totalCrossSectionCenter;

    //! Total electron cross section on the cell boundary (i-1/2) of
    //! energy grid
    vector_fp m_totalCrossSectionEdge;

    //! vector of total elastic cross section weighted with mass ratio
    vector_fp m_sigmaElastic;

    //! Set chemionization scattering-in rate
    double m_chemionScatRate;
};

}

#endif
