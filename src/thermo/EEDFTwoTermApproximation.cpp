/**
 *  @file EEDFTwoTermApproximation.cpp
 *  EEDF Two-Term approximation solver.  Implementation file for class
 *  EEDFTwoTermApproximation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/thermo/PlasmaPhase.h"
#include <iostream>

namespace Cantera
{

EEDFTwoTermApproximation::EEDFTwoTermApproximation(PlasmaPhase& s)
{
    initialize(s);
}

void EEDFTwoTermApproximation::initialize(PlasmaPhase& s)
{
    // store a pointer to s.
    m_phase = &s;
    m_first_call = true;
    m_has_EEDF = false;
    m_gamma = pow(2.0 * ElectronCharge / ElectronMass, 0.5);
}

void EEDFTwoTermApproximation::setLinearGrid(double& kTe_max, size_t& ncell)
{
    options.m_points = ncell;
    m_gridCenter.resize(options.m_points);
    m_gridEdge.resize(options.m_points + 1);
    m_f0.resize(options.m_points);
    m_f0_edge.resize(options.m_points + 1);
    for (size_t j = 0; j < options.m_points; j++) {
        m_gridCenter[j] = kTe_max * ( j + 0.5 ) / options.m_points;
        m_gridEdge[j] = kTe_max * j / options.m_points;
    }
    m_gridEdge[options.m_points] = kTe_max;
    setGridCache();
}

int EEDFTwoTermApproximation::calculateDistributionFunction()
{
    if (m_first_call) {
        initSpeciesIndexCS();
        m_first_call = false;
    }

    update_mole_fractions();
    checkSpeciesNoCrossSection();
    updateCS();

    if (!m_has_EEDF) {
        writelog("No existing EEDF. Using first guess method: {}\n", options.m_firstguess);
        if (options.m_firstguess == "maxwell") {
            writelog("First guess EEDF maxwell\n");
            for (size_t j = 0; j < options.m_points; j++) {
                m_f0(j) = 2.0 * pow(1.0 / Pi, 0.5) * pow(options.m_init_kTe, -3. / 2.) *
                          exp(-m_gridCenter[j] / options.m_init_kTe);
            }
        } else {
            throw CanteraError("EEDFTwoTermApproximation::calculateDistributionFunction",
                               " unknown EEDF first guess");
        }
    }

    converge(m_f0);

    // write the EEDF at grid edges
    vector<double> f(m_f0.data(), m_f0.data() + m_f0.rows() * m_f0.cols());
    vector<double> x(m_gridCenter.data(), m_gridCenter.data() + m_gridCenter.rows() * m_gridCenter.cols());
    for (size_t i = 0; i < options.m_points + 1; i++) {
        m_f0_edge[i] = linearInterp(m_gridEdge[i], x, f);
    }

    m_has_EEDF = true;

    // update electron mobility
    m_electronMobility = electronMobility(m_f0);

    return 0;

}

void EEDFTwoTermApproximation::converge(Eigen::VectorXd& f0)
{
    double err0 = 0.0;
    double err1 = 0.0;
    double delta = options.m_delta0;

    if (options.m_maxn == 0) {
        throw CanteraError("EEDFTwoTermApproximation::converge",
                           "options.m_maxn is zero; no iterations will occur.");
    }
    if (options.m_points == 0) {
        throw CanteraError("EEDFTwoTermApproximation::converge",
                           "options.m_points is zero; the EEDF grid is empty.");
    }
    if (std::isnan(delta) || delta == 0.0) {
        throw CanteraError("EEDFTwoTermApproximation::converge",
                           "options.m_delta0 is NaN or zero; solver cannot update.");
    }

    for (size_t n = 0; n < options.m_maxn; n++) {
        if (0.0 < err1 && err1 < err0) {
            delta *= log(options.m_factorM) / (log(err0) - log(err1));
        }

        Eigen::VectorXd f0_old = f0;
        f0 = iterate(f0_old, delta);

        err0 = err1;
        Eigen::VectorXd Df0(options.m_points);
        for (size_t i = 0; i < options.m_points; i++) {
            Df0(i) = abs(f0_old(i) - f0(i));
        }
        err1 = norm(Df0, m_gridCenter);

        if ((f0.array() != f0.array()).any()) {
            throw CanteraError("EEDFTwoTermApproximation::converge",
                               "NaN detected in EEDF solution.");
        }
        if ((f0.array().abs() > 1e300).any()) {
            throw CanteraError("EEDFTwoTermApproximation::converge",
                               "Inf detected in EEDF solution.");
        }

        if (err1 < options.m_rtol) {
            break;
        } else if (n == options.m_maxn - 1) {
            throw CanteraError("WeaklyIonizedGas::converge", "Convergence failed");
        }
    }
}

Eigen::VectorXd EEDFTwoTermApproximation::iterate(const Eigen::VectorXd& f0, double delta)
{
    // CQM multiple call to vector_* and matrix_*
    // probably extremely ineficient
    // must be refactored!!

    if ((f0.array() != f0.array()).any()) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
                           "NaN detected in input f0.");
    }

    SparseMat_fp PQ(options.m_points, options.m_points);
    vector_fp g = vector_g(f0);

    for (size_t k : m_phase->kInelastic()) {
        SparseMat_fp Q_k = matrix_Q(g, k);
        SparseMat_fp P_k = matrix_P(g, k);
        PQ += (matrix_Q(g, k) - matrix_P(g, k)) * m_X_targets[m_klocTargets[k]];
    }

    std::vector<std::tuple<int, int, double>> pq_values;
    int count = 0;
    for (int j = 0; j < PQ.outerSize(); ++j) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(PQ, j); it; ++it) {
            if (count < 5) {
                pq_values.push_back({it.row(), it.col(), it.value()});
            }
            count++;
        }
    }

    SparseMat_fp A = matrix_A(f0);
    SparseMat_fp I(options.m_points, options.m_points);
    for (size_t i = 0; i < options.m_points; i++) {
        I.insert(i,i) = 1.0;
    }
    A -= PQ;
    A *= delta;
    A += I;

    if (A.rows() == 0 || A.cols() == 0) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
                           "Matrix A has zero rows/columns.");
    }

    // SparseLU :
    Eigen::SparseLU<SparseMat_fp> solver(A);
    if (solver.info() == Eigen::NumericalIssue) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
            "Error SparseLU solver: NumericalIssue");
    } else if (solver.info() == Eigen::InvalidInput) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
            "Error SparseLU solver: InvalidInput");
    }
    if (solver.info() != Eigen::Success) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
            "Error SparseLU solver", "Decomposition failed");
        return f0;
    }

    // solve f0
    Eigen::VectorXd f1 = solver.solve(f0);
    if(solver.info() != Eigen::Success) {
        throw CanteraError("EEDFTwoTermApproximation::iterate", "Solving failed");
        return f0;
    }

    if ((f1.array() != f1.array()).any()) {
        throw CanteraError("EEDFTwoTermApproximation::iterate",
                           "NaN detected in computed f1.");
    }

    f1 /= norm(f1, m_gridCenter);

    return f1;
}

double EEDFTwoTermApproximation::integralPQ(double a, double b, double u0, double u1,
                                            double g, double x0)
{
    double A1;
    double A2;
    if (g != 0.0) {
        double expm1a = expm1(g * (-a + x0));
        double expm1b = expm1(g * (-b + x0));
        double ag = a * g;
        double ag1 = ag + 1;
        double bg = b * g;
        double bg1 = bg + 1;
        A1 = (expm1a * ag1 + ag - expm1b * bg1 - bg) / (g*g);
        A2 = (expm1a * (2 * ag1 + ag * ag) + ag * (ag + 2) -
              expm1b * (2 * bg1 + bg * bg) - bg * (bg + 2)) / (g*g*g);
    } else {
        A1 = 0.5 * (b*b - a*a);
        A2 = 1.0 / 3.0 * (b*b*b - a*a*a);
    }

    // The interpolation formula of u(x) = c0 + c1 * x
    double c0 = (a * u1 - b * u0) / (a - b);
    double c1 = (u0 - u1) / (a - b);

    return c0 * A1 + c1 * A2;
}

vector_fp EEDFTwoTermApproximation::vector_g(const Eigen::VectorXd& f0)
{
    vector_fp g(options.m_points, 0.0);
    const double f_min = 1e-300;  // Smallest safe floating-point value

    // Handle first point (i = 0)
    double f1 = std::max(f0(1), f_min);
    double f0_ = std::max(f0(0), f_min);
    g[0] = log(f1 / f0_) / (m_gridCenter[1] - m_gridCenter[0]);

    // Handle last point (i = N)
    size_t N = options.m_points - 1;
    double fN = std::max(f0(N), f_min);
    double fNm1 = std::max(f0(N - 1), f_min);
    g[N] = log(fN / fNm1) / (m_gridCenter[N] - m_gridCenter[N - 1]);

    // Handle interior points
    for (size_t i = 1; i < N; ++i) {
        double f_up   = std::max(f0(i + 1), f_min);
        double f_down = std::max(f0(i - 1), f_min);
        g[i] = log(f_up / f_down) / (m_gridCenter[i + 1] - m_gridCenter[i - 1]);
    }

    return g;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_P(const vector_fp& g, size_t k)
{
    vector<Triplet_fp> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        size_t j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double p = m_gamma * r;

        tripletList.push_back(Triplet_fp(j, j, p));
    }
    SparseMat_fp P(options.m_points, options.m_points);
    P.setFromTriplets(tripletList.begin(), tripletList.end());
    return P;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_Q(const vector_fp& g, size_t k)
{
    vector<Triplet_fp> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        size_t i = m_i[k][n];
        size_t j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double q = m_phase->inFactor()[k] * m_gamma * r;

        tripletList.push_back(Triplet_fp(i, j, q));
    }
    SparseMat_fp Q(options.m_points, options.m_points);
    Q.setFromTriplets(tripletList.begin(), tripletList.end());
    return Q;
}

SparseMat_fp EEDFTwoTermApproximation::matrix_A(const Eigen::VectorXd& f0)
{
    vector_fp a0(options.m_points + 1);
    vector_fp a1(options.m_points + 1);
    size_t N = options.m_points - 1;
    // Scharfetter-Gummel scheme
    double nu = netProductionFreq(f0);
    a0[0] = NAN;
    a1[0] = NAN;
    a0[N+1] = NAN;
    a1[N+1] = NAN;

    // Electron-electron collisions declarations
    double a;
    vector_fp A1, A2, A3;
    if (m_eeCol) {
        eeColIntegrals(A1, A2, A3, a, options.m_points);
    }

    double alpha;
    if (options.m_growth == "spatial") {
        double mu = electronMobility(f0);
        double D = electronDiffusivity(f0);
        alpha = (mu * m_phase->E() - sqrt(pow(mu * m_phase->E(), 2) - 4 * D * nu * m_phase->N())) / 2.0 / D / m_phase->N();
    } else {
        alpha = 0.0;
    }

    double sigma_tilde;
    double omega = 2 * Pi * m_phase->F();
    for (size_t j = 1; j < options.m_points; j++) {
        if (options.m_growth == "temporal") {
            sigma_tilde = m_totalCrossSectionEdge[j] + nu / pow(m_gridEdge[j], 0.5) / m_gamma;
        }
        else {
            sigma_tilde = m_totalCrossSectionEdge[j];
        }
        double q = omega / (m_phase->N() * m_gamma * pow(m_gridEdge[j], 0.5));
        double W = -m_gamma * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double F = sigma_tilde * sigma_tilde / (sigma_tilde * sigma_tilde + q * q);
        double DA = m_gamma / 3.0 * pow(m_phase->E() / m_phase->N(), 2.0) * m_gridEdge[j];
        double DB = m_gamma * m_phase->kT() * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double D = DA / sigma_tilde * F + DB;
        if (m_eeCol) {
            W -= 3 * a * m_phase->ionDegree() * A1[j];
            D += 2 * a * m_phase->ionDegree() * (A2[j] + pow(m_gridEdge[j], 1.5) * A3[j]);
        }
        if (options.m_growth == "spatial") {
            W -= m_gamma / 3.0 * 2 * alpha * m_phase->E() / m_phase->N() * m_gridEdge[j] / sigma_tilde;
        }

        double z = W * (m_gridCenter[j] - m_gridCenter[j-1]) / D;
        if (!std::isfinite(z)) {
            throw CanteraError("matrix_A", "Non-finite Peclet number encountered");
        }
        if (std::abs(z) > 500) {
            writelog("Warning: Large Peclet number z = {:.3e} at j = {}. W = {:.3e}, D = {:.3e}, E/N = {:.3e}\n",
                    z, j, W, D, m_phase->E() / m_phase->N());
        }
        a0[j] = W / (1 - std::exp(-z));
        a1[j] = W / (1 - std::exp(z));
    }

    std::vector<Triplet_fp> tripletList;
    // center diagonal
    // zero flux b.c. at energy = 0
    tripletList.push_back(Triplet_fp(0, 0, a0[1]));

    for (size_t j = 1; j < options.m_points - 1; j++) {
        tripletList.push_back(Triplet_fp(j, j, a0[j+1] - a1[j]));
    }

    // upper diagonal
    for (size_t j = 0; j < options.m_points - 1; j++) {
        tripletList.push_back(Triplet_fp(j, j+1, a1[j+1]));
    }

    // lower diagonal
    for (size_t j = 1; j < options.m_points; j++) {
        tripletList.push_back(Triplet_fp(j, j-1, -a0[j]));
    }

    // zero flux b.c.
    tripletList.push_back(Triplet_fp(N, N, -a1[N]));

    SparseMat_fp A(options.m_points, options.m_points);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    //plus G
    SparseMat_fp G(options.m_points, options.m_points);
    if (options.m_growth == "temporal") {
        for (size_t i = 0; i < options.m_points; i++) {
            G.insert(i, i) = 2.0 / 3.0 * (pow(m_gridEdge[i+1], 1.5) - pow(m_gridEdge[i], 1.5)) * nu;
        }
    }
    else if (options.m_growth == "spatial") {
        for (size_t i = 0; i < options.m_points; i++) {
            double sigma_c = 0.5 * (m_totalCrossSectionEdge[i] + m_totalCrossSectionEdge[i + 1]);
            G.insert(i, i) = - alpha * m_gamma / 3 * (alpha * (pow(m_gridEdge[i + 1], 2) - pow(m_gridEdge[i], 2)) / sigma_c / 2
                 - m_phase->E() / m_phase->N() * (m_gridEdge[i + 1] / m_totalCrossSectionEdge[i + 1] - m_gridEdge[i] / m_totalCrossSectionEdge[i]));
        }
    }

    return A + G;
}

double EEDFTwoTermApproximation::netProductionFreq(const Eigen::VectorXd& f0)
{
    double nu = 0.0;
    vector_fp g = vector_g(f0);

    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        if (m_phase->kind(k) == "ionization" ||
            m_phase->kind(k) == "attachment") {
            SparseMat_fp PQ = (matrix_Q(g, k) - matrix_P(g, k)) *
                              m_X_targets[m_klocTargets[k]];
            Eigen::VectorXd s = PQ * f0;
            for (size_t i = 0; i < options.m_points; i++) {
                nu += s[i];
                if (!std::isfinite(s[i])) {
                    writelog("NaN in netProductionFreq at s[{}] for k = {}\n", i, k);
                    break;
                }
            }
        }
    }
    return nu;
}

double EEDFTwoTermApproximation::electronDiffusivity(const Eigen::VectorXd& f0)
{
    vector_fp y(options.m_points, 0.0);
    double nu = netProductionFreq(f0);
    for (size_t i = 0; i < options.m_points; i++) {
        if (m_gridCenter[i] != 0.0) {
            y[i] = m_gridCenter[i] * f0(i) /
                   (m_totalCrossSectionCenter[i] + nu / m_gamma / pow(m_gridCenter[i], 0.5));
        }
    }
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridCenter.data(), m_gridCenter.size());
    return 1./3. * m_gamma * simpson(f, x) / m_phase->N();
}

double EEDFTwoTermApproximation::electronMobility(const Eigen::VectorXd& f0)
{
    double nu = netProductionFreq(f0);
    vector_fp y(options.m_points + 1, 0.0);
    for (size_t i = 1; i < options.m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (f0(i) - f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        if (m_gridEdge[i] != 0.0) {
            y[i] = m_gridEdge[i] * df0 /
                   (m_totalCrossSectionEdge[i] + nu / m_gamma / pow(m_gridEdge[i], 0.5));
        }
    }
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridEdge.data(), m_gridEdge.size());
    return -1./3. * m_gamma * simpson(f, x) / m_phase->N();
}

void EEDFTwoTermApproximation::initSpeciesIndexCS()
{
    // set up target index
    m_kTargets.resize(m_phase->nElectronCrossSections());
    m_klocTargets.resize(m_phase->nElectronCrossSections());
    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++)
    {

        m_kTargets[k] = m_phase->targetIndex(k);
        // Check if it is a new target or not :
        auto it = find(m_k_lg_Targets.begin(), m_k_lg_Targets.end(), m_kTargets[k]);

        if (it == m_k_lg_Targets.end()){
            m_k_lg_Targets.push_back(m_kTargets[k]);
            m_klocTargets[k] = m_k_lg_Targets.size() - 1;
        } else {
            m_klocTargets[k] = distance(m_k_lg_Targets.begin(), it);
        }
    }

    m_X_targets.resize(m_k_lg_Targets.size());
    m_X_targets_prev.resize(m_k_lg_Targets.size());
    for (size_t k = 0; k < m_X_targets.size(); k++)
    {
        size_t k_glob = m_k_lg_Targets[k];
        m_X_targets[k] = m_phase->moleFraction(k_glob);
        m_X_targets_prev[k] = m_phase->moleFraction(k_glob);
    }

    // set up indices of species which has no cross-section data
    for (size_t k = 0; k < m_phase->nSpecies(); k++)
    {
        auto it = std::find(m_kTargets.begin(), m_kTargets.end(), k);
        if (it == m_kTargets.end()) {
            m_kOthers.push_back(k);
        }
    }
}

void EEDFTwoTermApproximation::checkSpeciesNoCrossSection()
{
    // warn that a specific species needs cross-section data.
    for (size_t k : m_kOthers) {
        if (m_phase->moleFraction(k) > options.m_moleFractionThreshold) {
            writelog("EEDFTwoTermApproximation:checkSpeciesNoCrossSection\n");
            writelog("Warning:The mole fraction of species {} is more than 0.01 (X = {:.3g}) but it has no cross-section data\n", m_phase->speciesName(k), m_phase->moleFraction(k));
        }
    }
}

void EEDFTwoTermApproximation::updateCS()
{
    // Compute sigma_m and sigma_\epsilon
    calculateTotalCrossSection();
    calculateTotalElasticCrossSection();
}

// Update the species mole fractions used for EEDF computation
void EEDFTwoTermApproximation::update_mole_fractions()
{

    double tmp_sum = 0.0;
    for (size_t k = 0; k < m_X_targets.size(); k++)
    {
        m_X_targets[k] = m_phase->moleFraction(m_k_lg_Targets[k]);
        tmp_sum = tmp_sum + m_phase->moleFraction(m_k_lg_Targets[k]);
    }

    // Normalize the mole fractions to unity:
    for (size_t k = 0; k < m_X_targets.size(); k++)
    {
        m_X_targets[k] = m_X_targets[k] / tmp_sum;
    }

}

void EEDFTwoTermApproximation::calculateTotalCrossSection()
{

    m_totalCrossSectionCenter.assign(options.m_points, 0.0);
    m_totalCrossSectionEdge.assign(options.m_points + 1, 0.0);
    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        vector_fp x = m_phase->energyLevels()[k];
        vector_fp y = m_phase->crossSections()[k];

        for (size_t i = 0; i < options.m_points; i++) {
            m_totalCrossSectionCenter[i] += m_X_targets[m_klocTargets[k]] *
                                            linearInterp(m_gridCenter[i], x, y);

        }
        for (size_t i = 0; i < options.m_points + 1; i++) {
            m_totalCrossSectionEdge[i] += m_X_targets[m_klocTargets[k]] *
                                          linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void EEDFTwoTermApproximation::calculateTotalElasticCrossSection()
{
    m_sigmaElastic.clear();
    m_sigmaElastic.resize(options.m_points, 0.0);
    for (size_t k : m_phase->kElastic()) {
        vector_fp x = m_phase->energyLevels()[k];
        vector_fp y = m_phase->crossSections()[k];
        // Note:
        // moleFraction(m_kTargets[k]) <=> m_X_targets[m_klocTargets[k]]
        double mass_ratio = ElectronMass / (m_phase->molecularWeight(m_kTargets[k]) / Avogadro);
        for (size_t i = 0; i < options.m_points; i++) {
            m_sigmaElastic[i] += 2.0 * mass_ratio * m_X_targets[m_klocTargets[k]] *
                                 linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void EEDFTwoTermApproximation::setGridCache()
{
    m_sigma.clear();
    m_sigma.resize(m_phase->nElectronCrossSections());
    m_sigma_offset.clear();
    m_sigma_offset.resize(m_phase->nElectronCrossSections());
    m_eps.clear();
    m_eps.resize(m_phase->nElectronCrossSections());
    m_j.clear();
    m_j.resize(m_phase->nElectronCrossSections());
    m_i.clear();
    m_i.resize(m_phase->nElectronCrossSections());
    for (size_t k = 0; k < m_phase->nElectronCrossSections(); k++) {
        auto x = m_phase->energyLevels()[k];
        auto y = m_phase->crossSections()[k];
        vector_fp eps1(options.m_points + 1);
        for (size_t i = 0; i < options.m_points + 1; i++) {
            eps1[i] = clip(m_phase->shiftFactor()[k] * m_gridEdge[i] + m_phase->threshold(k),
                           m_gridEdge[0] + 1e-9, m_gridEdge[options.m_points] - 1e-9);
        }
        vector_fp nodes = eps1;
        for (size_t i = 0; i < options.m_points + 1; i++) {
            if (m_gridEdge[i] >= eps1[0] && m_gridEdge[i] <= eps1[options.m_points]) {
                nodes.push_back(m_gridEdge[i]);
            }
        }
        for (size_t i = 0; i < x.size(); i++) {
            if (x[i] >= eps1[0] && x[i] <= eps1[options.m_points]) {
                nodes.push_back(x[i]);
            }
        }

        std::sort(nodes.begin(), nodes.end());

        auto last = std::unique(nodes.begin(), nodes.end());
        nodes.resize(std::distance(nodes.begin(), last));
        vector_fp sigma0(nodes.size());
        for (size_t i = 0; i < nodes.size(); i++) {
            sigma0[i] = linearInterp(nodes[i], x, y);
        }

        // search position of cell j
        for (size_t i = 1; i < nodes.size(); i++) {
            auto low = std::lower_bound(m_gridEdge.begin(), m_gridEdge.end(), nodes[i]);
            m_j[k].push_back(low - m_gridEdge.begin() - 1);
        }

        // search position of cell i
        for (size_t i = 1; i < nodes.size(); i++) {
            auto low = std::lower_bound(eps1.begin(), eps1.end(), nodes[i]);
            m_i[k].push_back(low - eps1.begin() - 1);
        }

        // construct sigma
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            vector_fp sigma{sigma0[i], sigma0[i+1]};
            m_sigma[k].push_back(sigma);
        }

        // construct eps
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            vector_fp eps{nodes[i], nodes[i+1]};
            m_eps[k].push_back(eps);
        }

        // construct sigma_offset
        auto x_offset = m_phase->energyLevels()[k];
        for (auto& element : x_offset) {
            element -= m_phase->threshold(k);
        }
        for (size_t i = 0; i < options.m_points; i++) {
            m_sigma_offset[k].push_back(linearInterp(m_gridCenter[i], x_offset, y));
        }
    }
}

double EEDFTwoTermApproximation::norm(const Eigen::VectorXd& f, const Eigen::VectorXd& grid)
{
    string m_quadratureMethod = "simpson";
    Eigen::VectorXd p(f.size());
    for (int i = 0; i < f.size(); i++) {
        p[i] = f(i) * pow(grid[i], 0.5);
    }
    return numericalQuadrature(m_quadratureMethod, p, grid);
}


void EEDFTwoTermApproximation::eeColIntegrals(vector_fp& A1, vector_fp& A2, vector_fp& A3,
                                              double& a, size_t nPoints)
{
    // Ensure vectors are initialized
    A1.assign(nPoints, 0.0);
    A2.assign(nPoints, 0.0);
    A3.assign(nPoints, 0.0);

    // Compute net production frequency
    double nu = netProductionFreq(m_f0);
    // simulations with repeated calls to update EEDF will produce numerical instability here
    double nu_floor = 1e-40; // adjust as needed for stability
    if (nu < nu_floor) {
        writelog("eeColIntegrals: nu = {:.3e} too small, applying floor\n", nu);
        nu = nu_floor;
    }

    // Compute effective cross-section term
    double sigma_tilde;
    for (size_t j = 1; j < nPoints; j++) {
        sigma_tilde = m_totalCrossSectionCenter[j] + nu / pow(m_gridEdge[j], 0.5) / m_gamma;
    }

    // Compute Coulomb logarithm
    double lnLambda;
    if (nu > 0.0) {
        lnLambda = log(sigma_tilde / nu);
    } else {
        lnLambda = log(4.0 * Pi * pow(m_gridEdge.back(), 3) / 3.0);
    }

    // Compute e-e collision prefactor
    a = 4.0 * Pi * ElectronCharge * ElectronCharge * lnLambda / (m_gamma * pow(nu, 2));

    // Compute integral terms A1, A2, A3
    for (size_t j = 1; j < nPoints; j++) {
        double eps_j = m_gridCenter[j]; // Electron energy level
        double integral_A1 = 0.0;
        double integral_A2 = 0.0;
        double integral_A3 = 0.0;

        for (size_t i = 1; i < nPoints; i++) {
            double eps_i = m_gridCenter[i];
            double f0_i = m_f0[i];

            double weight = f0_i * pow(eps_i, 0.5) * exp(-abs(eps_i - eps_j) / eps_i);

            integral_A1 += weight * pow(eps_i, 1.5);
            integral_A2 += weight * pow(eps_i, 0.5);
            integral_A3 += weight;
        }

        // Store computed values
        A1[j] = integral_A1;
        A2[j] = integral_A2;
        A3[j] = integral_A3;
    }

}

}
