/**
 *  @file EEDFTwoTermApproximation.h EEDF Two-Term approximation solver.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_EEDF_TWO_TERM_APPROXIMATION_H
#define CT_EEDF_TWO_TERM_APPROXIMATION_H

#include "cantera/base/ct_defs.h"
#include "cantera/numerics/eigen_sparse.h"

namespace Cantera
{

class PlasmaPhase;

/**
 *  EEDF solver options. Used internally by class EEDFTwoTermApproximation.
 */
class TwoTermOpt
{
public:
    TwoTermOpt() = default;

    double m_delta0 = 1e14; //!< Initial value of the iteration parameter
    size_t m_maxn = 200; //!< Maximum number of iterations
    double m_factorM = 4.0; //!< Reduction factor of error
    size_t m_points = 150; //!< Number of points for energy grid
    double m_rtol = 1e-5; //!< Relative tolerance of EEDF for solving Boltzmann equation
    string m_growth = "temporal"; //!< String for the growth model (none, temporal or spatial)
    double m_moleFractionThreshold = 0.01; //!< Threshold for species not considered in the Boltzmann solver but present in the mixture
    string m_firstguess = "maxwell"; //!< String for EEDF first guess
    double m_init_kTe = 2.0; //!< Initial electron mean energy

}; // end of class TwoTermOpt

//! Boltzmann equation solver for the electron energy distribution function based on
//! the two-term approximation.
//!
//! @since New in %Cantera 3.2.
//! @warning This class is an experimental part of %Cantera and may be changed without
//!     notice.
class EEDFTwoTermApproximation
{
public:
    EEDFTwoTermApproximation() = default;

    //! Constructor combined with the initialization function
    /*!
     * This constructor initializes the EEDFTwoTermApproximation object with everything
     * it needs to start solving EEDF.
     *
     * @param s PlasmaPhase object that will be used in the solver calls.
     */
    EEDFTwoTermApproximation(PlasmaPhase& s);

    virtual ~EEDFTwoTermApproximation() = default;

    // compute the EEDF given an electric field
    // CQM The solver will take the species to consider and the set of cross-sections
    // from the PlasmaPhase object.
    // It will write the EEDF and its grid into the PlasmaPhase object.
    // Successful returns are indicated by a return value of 0.
    int calculateDistributionFunction();

    void setLinearGrid(double& kTe_max, size_t& ncell);

    void setGridCache();

    /**
     * Options controlling how the calculation is carried out.
     * @see TwoTermOpt
     */
    TwoTermOpt options;

    vector<double> getGridEdge() const {
        return m_gridEdge;
    }

    vector<double> getEEDFEdge() const {
        return m_f0_edge;
    }

    double getElectronMobility() const {
        return m_electronMobility;
    }

protected:
    //! Pointer to the PlasmaPhase object used to initialize this object.
    /*!
     * This PlasmaPhase object must be compatible with the PlasmaPhase objects
     * input from the compute function. Currently, this means that the 2
     * PlasmaPhases have to have consist of the same species and elements.
     */
    PlasmaPhase* m_phase;

    //! Iterate f0 (EEDF) until convergence
    void converge(Eigen::VectorXd& f0);

    //! An iteration of solving electron energy distribution function
    Eigen::VectorXd iterate(const Eigen::VectorXd& f0, double delta);

    //! The integral in [a, b] of \f$x u(x) \exp[g (x_0 - x)]\f$
    //! assuming that u is linear with u(a) = u0 and u(b) = u1
    double integralPQ(double a, double b, double u0, double u1,
                       double g, double x0);

    //! Vector g is used by matrix_P() and matrix_Q().
    /**
     * \f[
     * g_i = \frac{1}{\epsilon_{i+1} - \epsilon_{i-1}} \ln(\frac{F_{0, i+1}}{F_{0, i-1}})
     * \f]
     */
    vector<double> vector_g(const Eigen::VectorXd& f0);

    //! The matrix of scattering-out.
    /**
     * \f[
     * P_{i,k} = \gamma \int_{\epsilon_i - 1/2}^{\epsilon_i + 1/2}
     * \epsilon \sigma_k exp[(\epsilon_i - \epsilon)g_i] d \epsilon
     * \f]
     */
    Eigen::SparseMatrix<double> matrix_P(const vector<double>& g, size_t k);

    //! The matrix of scattering-in
    /**
     * \f[
     * Q_{i,j,k} = \gamma \int_{\epsilon_1}^{\epsilon_2}
     * \epsilon \sigma_k exp[(\epsilon_j - \epsilon)g_j] d \epsilon
     * \f]
     */
    //! where the interval \f$[\epsilon_1, \epsilon_2]\f$ is the overlap of cell j,
    //! and cell i shifted by the threshold energy:
    /**
     * \f[
     * \epsilon_1 = \min(\max(\epsilon_{i-1/2}+u_k, \epsilon_{j-1/2}),\epsilon_{j+1/2}),
     * \f]
     * \f[
     * \epsilon_2 = \min(\max(\epsilon_{i+1/2}+u_k, \epsilon_{j-1/2}),\epsilon_{j+1/2})
     * \f]
     */
    Eigen::SparseMatrix<double> matrix_Q(const vector<double>& g, size_t k);

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
    Eigen::SparseMatrix<double> matrix_A(const Eigen::VectorXd& f0);

    //! Reduced net production frequency. Equation (10) of ref. [1]
    //! divided by N.
    //! @param f0 EEDF
    double netProductionFreq(const Eigen::VectorXd& f0);

    //! Diffusivity
    double electronDiffusivity(const Eigen::VectorXd& f0);

    //! Mobility
    double electronMobility(const Eigen::VectorXd& f0);

    void initSpeciesIndexCS();

    void checkSpeciesNoCrossSection();

    void updateCS();

    void update_mole_fractions();

    void calculateTotalElasticCrossSection();

    void calculateTotalCrossSection();

    double norm(const Eigen::VectorXd& f, const Eigen::VectorXd& grid);

    double m_electronMobility;

    //! Grid of electron energy (cell center) [eV]
    Eigen::VectorXd m_gridCenter;

    //! Grid of electron energy (cell boundary i-1/2) [eV]
    vector<double> m_gridEdge;

    //! Location of cell j for grid cache
    vector<vector<size_t>> m_j;

    //! Location of cell i for grid cache
    vector<vector<size_t>> m_i;

    //! Cross section at the boundaries of the overlap of cell i and j
    vector<vector<vector<double>>> m_sigma;

    //! The energy boundaries of the overlap of cell i and j
    vector<vector<vector<double>>> m_eps;

    //! The cross section at the center of a cell
    vector<vector<double>> m_sigma_offset;

    //! normalized electron energy distribution function
    Eigen::VectorXd m_f0;

    //! EEDF at grid edges (cell boundaries)
    vector<double> m_f0_edge;

    //! Total electron cross section on the cell center of energy grid
    vector<double> m_totalCrossSectionCenter;

    //! Total electron cross section on the cell boundary (i-1/2) of
    //! energy grid
    vector<double> m_totalCrossSectionEdge;

    //! vector of total elastic cross section weighted with mass ratio
    vector<double> m_sigmaElastic;

    //! list of target species indices in global Cantera numbering (1 index per cs)
    vector<size_t> m_kTargets;

    //! list of target species indices in local X EEDF numbering (1 index per cs)
    vector<size_t> m_klocTargets;

    //! Indices of species which has no cross-section data
    vector<size_t> m_kOthers;

    //! Local to global indices
    vector<size_t> m_k_lg_Targets;

    //! Mole fraction of targets
    vector<double> m_X_targets;

    //! Previous mole fraction of targets used to compute eedf
    vector<double> m_X_targets_prev;

    //! in factor. This is used for calculating the Q matrix of
    //! scattering-in processes.
    vector<int> m_inFactor;

    double m_gamma;

    //! boolean for the electron-electron collisions
    bool m_eeCol = false;

    //! Compute electron-electron collision integrals
    void eeColIntegrals(vector<double>& A1, vector<double>& A2, vector<double>& A3,
                        double& a, size_t nPoints);

    //! flag of having an EEDF
    bool m_has_EEDF;

    //! First call to calculateDistributionFunction
    bool m_first_call;
}; // end of class EEDFTwoTermApproximation

} // end of namespace Cantera

#endif