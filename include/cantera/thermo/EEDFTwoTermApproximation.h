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

//! Boltzmann equation solver for the electron energy distribution function based on
//! the two-term approximation.
/*!
 * This class implements a solver for the electron energy distribution function
 * based on a steady-state solution to the Boltzmann equation using the classical
 * two-term expansion, applicable to weakly ionized plasmas. The numerical approach
 * and theory are primarily derived from the work of Hagelaar and Pitchford
 * @cite hagelaar2005.
 *
 * The two-term approximation assumes that the EEDF can be represented as:
 *   @f[
 *       f(\epsilon, \mu) = f_0(\epsilon) + \mu f_1(\epsilon),
 *   @f]
 * where @f$ \epsilon @f$ is the electron energy and @f$ \mu @f$ is the cosine
 * of the angle between the electron velocity vector and the electric field.
 * The Boltzmann equation is projected onto the zeroth moment over mu to obtain
 * an equation for @f$ f_0(\epsilon) @f$ , the isotropic part of the distribution.
 * The first-order anisotropic term @f$ f_1(\epsilon) @f$ is not solved directly,
 * but is approximated and substituted into the drift and collision terms. This
 * results in a second-order differential equation for @f$ f_0(\epsilon) @f$ alone.
 *
 * The governing equation for @f$ f_0(\epsilon) @f$ is discretized on an energy
 * grid using a finite difference method and solved using a tridiagonal matrix
 * algorithm.
 *
 * @since New in %Cantera 3.2.
 * @warning This class is an experimental part of %Cantera and may be changed without
 *     notice.
 */
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
    EEDFTwoTermApproximation(PlasmaPhase* s);

    virtual ~EEDFTwoTermApproximation() = default;

    //! compute the EEDF given an electric field
    //! CQM The solver will take the species to consider and the set of cross-sections
    //! from the PlasmaPhase object.
    //! It will write the EEDF and its grid into the PlasmaPhase object.
    //! Successful returns are indicated by a return value of 0.
    int calculateDistributionFunction();

    //! Sets a linear energy grid for the EEDF solver, defined by the maximum energy and the number of grid cells.
    //! @since New in %Cantera 4.0
    //! @param kTe_max: maximum grid energy in eV
    //! @param ncell: number of cell to discretize the grid.
    void setLinearGrid(double kTe_max, size_t ncell);
    
    //! Sets a quadratic energy grid for the EEDF solver, defined by the maximum energy and the number of grid cells.
    //! @since New in %Cantera 4.0
    //! @param kTe_max: maximum grid energy in eV
    //! @param ncell: number of cell to discretize the grid.
    void setQuadraticGrid(double kTe_max, size_t ncell);

    //! Sets a geometric energy grid for the EEDF solver, defined by the maximum energy and the number of grid cells.
    //! @since New in %Cantera 4.0
    //! @param kTe_max: maximum grid energy in eV
    //! @param ncell: number of cell to discretize the grid.
    //! @param ratio: the geometric growth ratio.
    void setGeometricGrid(double kTe_max, size_t ncell, double ratio = 1.01);

    //! Sets a custom energy grid for the EEDF solver, defined by the user-provided vector of energy levels.
    void setCustomGrid(span<const double> levels);

    //! Build or rebuild the grid-dependent cache used for scattering matrices.
    void setGridCache();

    //! Set the initial grid parameters used by generated EEDF grids.
    /*!
     * These values are used when creating linear, quadratic, or geometric grids,
     * and are also reused during grid adaptation.
     *
     * @param initialMaxEnergy  Maximum electron energy of the initial grid [eV].
     * @param nGridCells        Number of grid cells. The number of grid edges is
     *                          nGridCells + 1.
     * @param gridType          Type of grid spacing to use when generating or adapting
     *                          the electron energy grid. Supported values are "linear", 
     *                          "quadratic", and "geometric".
     */
    //! @since New in %Cantera 4.0
    void setInitialGridParameters(double initialMaxEnergy, size_t nGridCells, const string& gridType);

    //! Enable or disable automatic grid adaptation for the EEDF solver energy grid.
    //! @since New in %Cantera 4.0
    void enableGridAdaptation(bool enabled);

    //! Set parameters controlling automatic adaptation of the EEDF energy grid.
    /*!
     * Grid adaptation adjusts the maximum grid energy based on the number of
     * decades by which the EEDF decays between the low- and high-energy ends of
     * the grid. The number of grid cells is kept fixed.
     *
     * @param minDecayDecades  Minimum acceptable EEDF tail decay in decades. If
     *                         the decay is smaller, the maximum grid energy is
     *                         increased.
     * @param maxDecayDecades  Maximum acceptable EEDF tail decay in decades. If
     *                         the decay is larger, the maximum grid energy is
     *                         decreased.
     * @param updateFactor     Relative factor used to increase or decrease the
     *                         maximum grid energy during adaptation.
     * @param maxIterations    Maximum number of grid adaptation iterations per
     *                         EEDF solve.
     */
    //! @since New in %Cantera 4.0
    void setGridAdaptationParameters(double minDecayDecades, double maxDecayDecades, double updateFactor, size_t maxIterations);

    //! Return the electron energy grid edges [eV].
    /*!
     * The returned vector contains m_points + 1 values corresponding to cell
     * boundaries.
     */
    //! @since New in %Cantera 4.0
    vector<double> getGridEdge() const {
        return m_gridEdge;
    }

    //! Return the EEDF values interpolated at the electron energy grid edges.
    /*!
     * These values are copied back to PlasmaPhase after a successful
     * Boltzmann-two-term EEDF solve.
     */
    //! @since New in %Cantera 4.0
    vector<double> getEEDFEdge() const {
        return m_f0_edge;
    }

    //! Return the latest value of the computed electron mobility computed from the EEDF
    //! @since New in %Cantera 4.0
    double getElectronMobility() const {
        return m_electronMobility;
    }
    //! Sets the threshold in reduced electric field below which a Maxwellian is imposed instead of computing the EEDF.
    //! @param threshold    The threshold in Td.
    //! @since New in %Cantera 4.0
    void setReducedElectricFieldThresholdForMaxwellian(double threshold){
        m_threshold_to_maxwellian = threshold;
    }

protected:

    //! Formerly options for the EEDF solver

    //! The first step size
    double m_delta0 = 1e14;

    //! Maximum number of iterations
    size_t m_maxn = 200;

    //! The factor for step size change
    double m_factorM = 4.0;

    //! The number of points in the EEDF grid
    size_t m_points = 150;

    //! Error tolerance for convergence
    double m_rtol = 1e-5;

    //! The growth model of EEDF
    std::string m_growth = "temporal";

    //! The threshold for species mole fractions
    double m_moleFractionThreshold = 0.01;

    //! The initial electron temperature [eV]
    double m_init_kTe = 2.0;

    //! Pointer to the PlasmaPhase object used to initialize this object.
    /*!
     * This PlasmaPhase object provides species, element, and cross-section
     * data used by the EEDF solver. It is set during construction and is not
     * modified afterwards. All subsequent calls to compute functions must
     * use the same PlasmaPhase context.
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
    Eigen::SparseMatrix<double> matrix_P(span<const double> g, size_t k);

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
    Eigen::SparseMatrix<double> matrix_Q(span<const double> g, size_t k);

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
    double netProductionFrequency(const Eigen::VectorXd& f0);

    //! Diffusivity
    double electronDiffusivity(const Eigen::VectorXd& f0);

    //! Mobility
    double electronMobility(const Eigen::VectorXd& f0);

    //! Initialize species indices associated with cross-section data
    void initSpeciesIndexCrossSections();

    //! Update the total cross sections based on the current state
    void updateCrossSections();

    //! Update the vector of species mole fractions
    void updateMoleFractions();

    //! Compute the total elastic collision cross section
    void calculateTotalElasticCrossSection();

    //! Compute the total (elastic + inelastic) cross section
    void calculateTotalCrossSection();

    //! Compute the L1 norm of a function f defined over a given energy grid.
    //!
    //! @param f     Vector representing the function values (EEDF)
    //! @param grid  Vector representing the energy grid corresponding to f
    double norm(const Eigen::VectorXd& f, const Eigen::VectorXd& grid);

    //! Electron mobility [m²/V·s]
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

    //! Normalized electron energy distribution function
    Eigen::VectorXd m_f0;

    //! EEDF at grid edges (cell boundaries)
    vector<double> m_f0_edge;

    //! Total electron cross section on the cell center of energy grid
    vector<double> m_totalCrossSectionCenter;

    //! Total electron cross section on the cell boundary (i-1/2) of
    //! energy grid
    vector<double> m_totalCrossSectionEdge;

    //! Vector of total elastic cross section weighted with mass ratio
    vector<double> m_sigmaElastic;

    //! List of target species indices in global Cantera numbering (1 index per cs)
    vector<size_t> m_kTargets;

    //! List of target species indices in local X EEDF numbering (1 index per cs)
    vector<size_t> m_klocTargets;

    //! Indices of species which has no cross-section data
    vector<size_t> m_kOthers;

    //! Local to global indices
    vector<size_t> m_k_lg_Targets;

    //! Mole fraction of targets
    vector<double> m_X_targets;

    //! Previous mole fraction of targets used to compute eedf
    vector<double> m_X_targets_prev;

    //! In factor. This is used for calculating the Q matrix of
    //! scattering-in processes.
    vector<int> m_inFactor;

    //! Defined by the formula: pow(2.0 * ElectronCharge / ElectronMass, 0.5) and comupted during phase initilisation.
    double m_gamma;

    //! Flag of having an EEDF
    bool m_has_EEDF;

    //! First call to calculateDistributionFunction
    bool m_first_call;
    
    //! Energy grid spacing type. Initialised with linear but can also be quadratic or geometric.
    string m_gridType = "linear";

    //! Maximum value of the energy grid in eV. Initialised at 60 eV. [eV]
    double m_kTeMax = 60.0;

    //! Number of cells for the starting energy grid. Initialised at 301.
    size_t m_initialGridCells = 301;

    //! Flag activating or deactivating automatic grid adaptation. Initialised at false.
    bool m_adaptGrid = false;

    //! Minimum amount of decades decay at the tail of the EEDF when grid adaptation is on. Initialised at 10.
    double m_minEedfDecay = 10;

    //! Maximum amount of decades decay at the tail of the EEDF when grid adaptation is on. Initialised at 12.
    double m_maxEedfDecay = 12.0;

    //! Factor by which the EEDF grid maximum energy is increased of shrunk when grid adaptation is on. initialised at 0.1.
    double m_gridUpdateFactor = 0.1;

    //! Maximum number of iterations on the maximum energy accepted for grid adaptation. Initialised at 1000.
    size_t m_maxGridAdaptIterations = 1000;

    //! Updates the grid according to the grid type and the new maximum energy when running grid adaptation. 
    //! @since New in %Cantera 4.0
    void updateGrid(double maxEnergy);

    //! Runs the energy grid adaptation script when this feature is activated.
    //! @since New in %Cantera 4.0
    void adaptEnergyGrid();

    //! Sets a Maxwellian distribution on the 
    void setMaxwellianDistribution(double kTe);

    //! The threshold in reduced electric field below which no EEDF will be computed,
    //! but a Maxwellian at the gas temperature will be imposed instead. 
    //! It is expressed in Townsend, and defaults to 1 Td. [Td]
    double m_threshold_to_maxwellian = 1;

    //! In the case where a geometric grid is employed, this stores the corresponding geometric ratio.
    double m_geometric_ratio = 1.01;
}; // end of class EEDFTwoTermApproximation

} // end of namespace Cantera

#endif