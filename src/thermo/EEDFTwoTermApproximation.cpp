/**
 *  @file EEDFTwoTermApproximation.cpp
 *  EEDF Two-Term approximation solver.  Implementation file for class
 *  EEDFTwoTermApproximation.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/EEDFTwoTermApproximation.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/numerics/funcs.h"
#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"

namespace Cantera
{

typedef Eigen::SparseMatrix<double> SparseMat;

EEDFTwoTermApproximation::EEDFTwoTermApproximation(PlasmaPhase* s)
{
    // store a pointer to s.
    m_phase = s;
    m_first_call = true;
    m_has_EEDF = false;
    m_gamma = pow(2.0 * ElectronCharge / ElectronMass, 0.5);
}

void EEDFTwoTermApproximation::setLinearGrid(double& kTe_max, size_t& ncell)
{
    m_points = ncell;
    m_gridCenter.resize(m_points);
    m_gridEdge.resize(m_points + 1);
    m_f0.resize(m_points);
    m_f0_edge.resize(m_points + 1);
    for (size_t j = 0; j < m_points; j++) {
        m_gridCenter[j] = kTe_max * ( j + 0.5 ) / m_points;
        m_gridEdge[j] = kTe_max * j / m_points;
    }
    m_gridEdge[m_points] = kTe_max;
    setGridCache();
}

int EEDFTwoTermApproximation::calculateDistributionFunction()
{
    if (m_first_call) {
        initSpeciesIndexCrossSections();
        m_first_call = false;
    }

    updateMoleFractions();
    updateCrossSections();

    if (!m_has_EEDF) {
        if (m_firstguess == "maxwell") {
            for (size_t j = 0; j < m_points; j++) {
                m_f0(j) = 2.0 * pow(1.0 / Pi, 0.5) * pow(m_init_kTe, -3. / 2.) *
                          exp(-m_gridCenter[j] / m_init_kTe);
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
    for (size_t i = 0; i < m_points + 1; i++) {
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
    double delta = m_delta0;

    if (m_maxn == 0) {
        throw CanteraError("EEDFTwoTermApproximation::converge",
                           "m_maxn is zero; no iterations will occur.");
    }
    if (m_points == 0) {
        throw CanteraError("EEDFTwoTermApproximation::converge",
                           "m_points is zero; the EEDF grid is empty.");
    }
    if (isnan(delta) || delta == 0.0) {
        throw CanteraError("EEDFTwoTermApproximation::converge",
                           "m_delta0 is NaN or zero; solver cannot update.");
    }

    for (size_t n = 0; n < m_maxn; n++) {
        if (0.0 < err1 && err1 < err0) {
            delta *= log(m_factorM) / (log(err0) - log(err1));
        }

        Eigen::VectorXd f0_old = f0;
        f0 = iterate(f0_old, delta);
        checkFinite("EEDFTwoTermApproximation::converge: f0", f0.data(), f0.size());

        err0 = err1;
        Eigen::VectorXd Df0 = (f0_old - f0).cwiseAbs();
        err1 = norm(Df0, m_gridCenter);
        if (err1 < m_rtol) {
            break;
        } else if (n == m_maxn - 1) {
            throw CanteraError("WeaklyIonizedGas::converge", "Convergence failed; returning last iterate.\n");
            return;
        }
    }
}

Eigen::VectorXd EEDFTwoTermApproximation::iterate(const Eigen::VectorXd& f0, double delta)
{
    // CQM multiple call to vector_* and matrix_*
    // probably extremely ineficient
    // must be refactored!!

    SparseMat PQ(m_points, m_points);
    vector<double> g = vector_g(f0);

    for (size_t k : m_phase->kInelastic()) {
        SparseMat Q_k = matrix_Q(g, k);
        SparseMat P_k = matrix_P(g, k);
        PQ += (matrix_Q(g, k) - matrix_P(g, k)) * m_X_targets[m_klocTargets[k]];
    }

    SparseMat A = matrix_A(f0);
    SparseMat I(m_points, m_points);
    for (size_t i = 0; i < m_points; i++) {
        I.insert(i,i) = 1.0;
    }
    A -= PQ;
    A *= delta;
    A += I;

    // SparseLU :
    Eigen::SparseLU<SparseMat> solver(A);
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

    checkFinite("EEDFTwoTermApproximation::converge: f0", f1.data(), f1.size());
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

vector<double> EEDFTwoTermApproximation::vector_g(const Eigen::VectorXd& f0)
{
    vector<double> g(m_points, 0.0);
    const double f_min = 1e-300;  // Smallest safe floating-point value

    // Handle first point (i = 0)
    double f1 = std::max(f0(1), f_min);
    double f0_ = std::max(f0(0), f_min);
    g[0] = log(f1 / f0_) / (m_gridCenter[1] - m_gridCenter[0]);

    // Handle last point (i = N)
    size_t N = m_points - 1;
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

SparseMat EEDFTwoTermApproximation::matrix_P(const vector<double>& g, size_t k)
{
    SparseTriplets tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        auto j = static_cast<SparseMat::StorageIndex>(m_j[k][n]);
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double p = m_gamma * r;

        tripletList.emplace_back(j, j, p);
    }
    SparseMat P(m_points, m_points);
    P.setFromTriplets(tripletList.begin(), tripletList.end());
    return P;
}

SparseMat EEDFTwoTermApproximation::matrix_Q(const vector<double>& g, size_t k)
{
    SparseTriplets tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        auto i = static_cast<SparseMat::StorageIndex>(m_i[k][n]);
        auto j = static_cast<SparseMat::StorageIndex>(m_j[k][n]);
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double q = m_inFactor[k] * m_gamma * r;

        tripletList.emplace_back(i, j, q);
    }
    SparseMat Q(m_points, m_points);
    Q.setFromTriplets(tripletList.begin(), tripletList.end());
    return Q;
}

SparseMat EEDFTwoTermApproximation::matrix_A(const Eigen::VectorXd& f0)
{
    vector<double> a0(m_points + 1);
    vector<double> a1(m_points + 1);
    size_t N = m_points - 1;
    // Scharfetter-Gummel scheme
    double nu = netProductionFrequency(f0);
    a0[0] = NAN;
    a1[0] = NAN;
    a0[N+1] = NAN;
    a1[N+1] = NAN;

    double nDensity = m_phase->molarDensity() * Avogadro;
    double alpha;
    double E = m_phase->electricField();
    if (m_growth == "spatial") {
        double mu = electronMobility(f0);
        double D = electronDiffusivity(f0);
        alpha = (mu * E - sqrt(pow(mu * E, 2) - 4 * D * nu * nDensity)) / 2.0 / D / nDensity;
    } else {
        alpha = 0.0;
    }

    double sigma_tilde;
    double omega = 2 * Pi * m_phase->electricFieldFrequency();
    for (size_t j = 1; j < m_points; j++) {
        if (m_growth == "temporal") {
            sigma_tilde = m_totalCrossSectionEdge[j] + nu / pow(m_gridEdge[j], 0.5) / m_gamma;
        } else {
            sigma_tilde = m_totalCrossSectionEdge[j];
        }
        double q = omega / (nDensity * m_gamma * pow(m_gridEdge[j], 0.5));
        double W = -m_gamma * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double F = sigma_tilde * sigma_tilde / (sigma_tilde * sigma_tilde + q * q);
        double DA = m_gamma / 3.0 * pow(E / nDensity, 2.0) * m_gridEdge[j];
        double DB = m_gamma * m_phase->temperature() * Boltzmann / ElectronCharge * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double D = DA / sigma_tilde * F + DB;
        if (m_growth == "spatial") {
            W -= m_gamma / 3.0 * 2 * alpha * E / nDensity * m_gridEdge[j] / sigma_tilde;
        }

        double z = W * (m_gridCenter[j] - m_gridCenter[j-1]) / D;
        if (!std::isfinite(z)) {
            throw CanteraError("matrix_A", "Non-finite Peclet number encountered");
        }
        if (std::abs(z) > 500) {
            warn_user("EEDFTwoTermApproximation::matrix_A",
                "Large Peclet number z = {:.3e} at j = {}. "
                "W = {:.3e}, D = {:.3e}, E/N = {:.3e}\n",
                z, j, W, D, E / nDensity);
        }
        a0[j] = W / (1 - std::exp(-z));
        a1[j] = W / (1 - std::exp(z));
    }

    SparseTriplets tripletList;
    // center diagonal
    // zero flux b.c. at energy = 0
    tripletList.emplace_back(0, 0, a0[1]);

    for (size_t j = 1; j < m_points - 1; j++) {
        tripletList.emplace_back(j, j, a0[j+1] - a1[j]);
    }

    // upper diagonal
    for (size_t j = 0; j < m_points - 1; j++) {
        tripletList.emplace_back(j, j+1, a1[j+1]);
    }

    // lower diagonal
    for (size_t j = 1; j < m_points; j++) {
        tripletList.emplace_back(j, j-1, -a0[j]);
    }

    // zero flux b.c.
    tripletList.emplace_back(N, N, -a1[N]);

    SparseMat A(m_points, m_points);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    //plus G
    SparseMat G(m_points, m_points);
    if (m_growth == "temporal") {
        for (size_t i = 0; i < m_points; i++) {
            G.insert(i, i) = 2.0 / 3.0 * (pow(m_gridEdge[i+1], 1.5) - pow(m_gridEdge[i], 1.5)) * nu;
        }
    } else if (m_growth == "spatial") {
        double nDensity = m_phase->molarDensity() * Avogadro;
        for (size_t i = 0; i < m_points; i++) {
            double sigma_c = 0.5 * (m_totalCrossSectionEdge[i] + m_totalCrossSectionEdge[i + 1]);
            G.insert(i, i) = - alpha * m_gamma / 3 * (alpha * (pow(m_gridEdge[i + 1], 2) - pow(m_gridEdge[i], 2)) / sigma_c / 2
                 - E / nDensity * (m_gridEdge[i + 1] / m_totalCrossSectionEdge[i + 1] - m_gridEdge[i] / m_totalCrossSectionEdge[i]));
        }
    }
    return A + G;
}

double EEDFTwoTermApproximation::netProductionFrequency(const Eigen::VectorXd& f0)
{
    double nu = 0.0;
    vector<double> g = vector_g(f0);

    for (size_t k = 0; k < m_phase->nCollisions(); k++) {
        if (m_phase->collisionRate(k)->kind() == "ionization" ||
            m_phase->collisionRate(k)->kind() == "attachment") {
            SparseMat PQ = (matrix_Q(g, k) - matrix_P(g, k)) *
                              m_X_targets[m_klocTargets[k]];
            Eigen::VectorXd s = PQ * f0;
            checkFinite("EEDFTwoTermApproximation::netProductionFrequency: s",
                        s.data(), s.size());
            nu += s.sum();
        }
    }
    return nu;
}

double EEDFTwoTermApproximation::electronDiffusivity(const Eigen::VectorXd& f0)
{
    vector<double> y(m_points, 0.0);
    double nu = netProductionFrequency(f0);
    for (size_t i = 0; i < m_points; i++) {
        if (m_gridCenter[i] != 0.0) {
            y[i] = m_gridCenter[i] * f0(i) /
                   (m_totalCrossSectionCenter[i] + nu / m_gamma / pow(m_gridCenter[i], 0.5));
        }
    }
    double nDensity = m_phase->molarDensity() * Avogadro;
    auto f = Eigen::Map<const Eigen::ArrayXd>(y.data(), y.size());
    auto x = Eigen::Map<const Eigen::ArrayXd>(m_gridCenter.data(), m_gridCenter.size());
    return 1./3. * m_gamma * simpson(f, x) / nDensity;
}

double EEDFTwoTermApproximation::electronMobility(const Eigen::VectorXd& f0)
{
    double nu = netProductionFrequency(f0);
    vector<double> y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (f0(i) - f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        if (m_gridEdge[i] != 0.0) {
            y[i] = m_gridEdge[i] * df0 /
                   (m_totalCrossSectionEdge[i] + nu / m_gamma / pow(m_gridEdge[i], 0.5));
        }
    }
    double nDensity = m_phase->molarDensity() * Avogadro;
    auto f = ConstMappedVector(y.data(), y.size());
    auto x = ConstMappedVector(m_gridEdge.data(), m_gridEdge.size());
    return -1./3. * m_gamma * simpson(f, x) / nDensity;
}

void EEDFTwoTermApproximation::initSpeciesIndexCrossSections()
{
    // set up target index
    m_kTargets.resize(m_phase->nCollisions());
    m_klocTargets.resize(m_phase->nCollisions());
    m_inFactor.resize(m_phase->nCollisions());
    for (size_t k = 0; k < m_phase->nCollisions(); k++) {
        m_kTargets[k] = m_phase->targetIndex(k);
        // Check if it is a new target or not :
        auto it = find(m_k_lg_Targets.begin(), m_k_lg_Targets.end(), m_kTargets[k]);

        if (it == m_k_lg_Targets.end()){
            m_k_lg_Targets.push_back(m_kTargets[k]);
            m_klocTargets[k] = m_k_lg_Targets.size() - 1;
        } else {
            m_klocTargets[k] = distance(m_k_lg_Targets.begin(), it);
        }

        const auto& kind = m_phase->collisionRate(k)->kind();

        if (kind == "ionization") {
            m_inFactor[k] = 2;
        } else if (kind == "attachment") {
            m_inFactor[k] = 0;
        } else {
            m_inFactor[k] = 1;
        }
    }

    m_X_targets.resize(m_k_lg_Targets.size());
    m_X_targets_prev.resize(m_k_lg_Targets.size());
    for (size_t k = 0; k < m_X_targets.size(); k++) {
        size_t k_glob = m_k_lg_Targets[k];
        m_X_targets[k] = m_phase->moleFraction(k_glob);
        m_X_targets_prev[k] = m_phase->moleFraction(k_glob);
    }

    // set up indices of species which has no cross-section data
    for (size_t k = 0; k < m_phase->nSpecies(); k++) {
        auto it = std::find(m_kTargets.begin(), m_kTargets.end(), k);
        if (it == m_kTargets.end()) {
            m_kOthers.push_back(k);
        }
    }
}

void EEDFTwoTermApproximation::updateCrossSections()
{
    // Compute sigma_m and sigma_\epsilon
    calculateTotalCrossSection();
    calculateTotalElasticCrossSection();
}

// Update the species mole fractions used for EEDF computation
void EEDFTwoTermApproximation::updateMoleFractions()
{
    double tmp_sum = 0.0;
    for (size_t k = 0; k < m_X_targets.size(); k++) {
        m_X_targets[k] = m_phase->moleFraction(m_k_lg_Targets[k]);
        tmp_sum = tmp_sum + m_phase->moleFraction(m_k_lg_Targets[k]);
    }

    // Normalize the mole fractions to unity:
    for (size_t k = 0; k < m_X_targets.size(); k++) {
        m_X_targets[k] = m_X_targets[k] / tmp_sum;
    }
}

void EEDFTwoTermApproximation::calculateTotalCrossSection()
{
    m_totalCrossSectionCenter.assign(m_points, 0.0);
    m_totalCrossSectionEdge.assign(m_points + 1, 0.0);
    for (size_t k = 0; k < m_phase->nCollisions(); k++) {
        auto& x = m_phase->collisionRate(k)->energyLevels();
        auto& y = m_phase->collisionRate(k)->crossSections();

        for (size_t i = 0; i < m_points; i++) {
            m_totalCrossSectionCenter[i] += m_X_targets[m_klocTargets[k]] *
                                            linearInterp(m_gridCenter[i], x, y);
        }
        for (size_t i = 0; i < m_points + 1; i++) {
            m_totalCrossSectionEdge[i] += m_X_targets[m_klocTargets[k]] *
                                          linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void EEDFTwoTermApproximation::calculateTotalElasticCrossSection()
{
    m_sigmaElastic.clear();
    m_sigmaElastic.resize(m_points, 0.0);
    for (size_t k : m_phase->kElastic()) {
        auto& x = m_phase->collisionRate(k)->energyLevels();
        auto& y = m_phase->collisionRate(k)->crossSections();
        // Note:
        // moleFraction(m_kTargets[k]) <=> m_X_targets[m_klocTargets[k]]
        double mass_ratio = ElectronMass / (m_phase->molecularWeight(m_kTargets[k]) / Avogadro);
        for (size_t i = 0; i < m_points; i++) {
            m_sigmaElastic[i] += 2.0 * mass_ratio * m_X_targets[m_klocTargets[k]] *
                                 linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void EEDFTwoTermApproximation::setGridCache()
{
    m_sigma.clear();
    m_sigma.resize(m_phase->nCollisions());
    m_eps.clear();
    m_eps.resize(m_phase->nCollisions());
    m_j.clear();
    m_j.resize(m_phase->nCollisions());
    m_i.clear();
    m_i.resize(m_phase->nCollisions());
    for (size_t k = 0; k < m_phase->nCollisions(); k++) {
        auto& collision = m_phase->collisionRate(k);
        auto& x = collision->energyLevels();
        auto& y = collision->crossSections();
        vector<double> eps1(m_points + 1);
        int shiftFactor = (collision->kind() == "ionization") ? 2 : 1;

        for (size_t i = 0; i < m_points + 1; i++) {
            eps1[i] = clip(shiftFactor * m_gridEdge[i] + collision->threshold(),
                           m_gridEdge[0] + 1e-9, m_gridEdge[m_points] - 1e-9);
        }
        vector<double> nodes = eps1;
        for (size_t i = 0; i < m_points + 1; i++) {
            if (m_gridEdge[i] >= eps1[0] && m_gridEdge[i] <= eps1[m_points]) {
                nodes.push_back(m_gridEdge[i]);
            }
        }
        for (size_t i = 0; i < x.size(); i++) {
            if (x[i] >= eps1[0] && x[i] <= eps1[m_points]) {
                nodes.push_back(x[i]);
            }
        }

        std::sort(nodes.begin(), nodes.end());
        auto last = std::unique(nodes.begin(), nodes.end());
        nodes.resize(std::distance(nodes.begin(), last));
        vector<double> sigma0(nodes.size());
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
            m_sigma[k].push_back({sigma0[i], sigma0[i+1]});
        }

        // construct eps
        for (size_t i = 0; i < nodes.size() - 1; i++) {
            m_eps[k].push_back({nodes[i], nodes[i+1]});
        }

        // construct sigma_offset
        auto x_offset = collision->energyLevels();
        for (auto& element : x_offset) {
            element -= collision->threshold();
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

}
