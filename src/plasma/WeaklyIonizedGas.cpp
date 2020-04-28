// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/plasma/WeaklyIonizedGas.h"
#include "cantera/base/utilities.h"
#include "cantera/numerics/funcs.h"
#include <Eigen/SparseLU>

namespace Cantera {

typedef Eigen::Triplet<double> Triplet;

//! Calculate the norm of EEDF
double norm(const Eigen::VectorXd& f, const vector_fp& grid)
{
    vector_fp p(f.size());
    for (int i = 0; i < f.size(); i++) {
        p[i] = f(i) * pow(grid[i], 0.5);
    }
    return simpsonQuadrature(grid, p);
}

WeaklyIonizedGas::WeaklyIonizedGas()
    : m_chemionScatRate(0.0)
{
}

void WeaklyIonizedGas::calculateTotalCrossSection()
{
    m_totalCrossSectionCenter.assign(m_points, 0.0);
    m_totalCrossSectionEdge.assign(m_points + 1, 0.0);
    for (size_t k = 0; k < m_ncs; k++) {
        vector_fp& x = m_crossSections[k][0];
        vector_fp& y = m_crossSections[k][1];
        for (size_t i = 0; i < m_points; i++) {
            m_totalCrossSectionCenter[i] += moleFraction(m_kTargets[k]) *
                                            linearInterp(m_gridCenter[i], x, y);
        }
        for (size_t i = 0; i < m_points + 1; i++) {
            m_totalCrossSectionEdge[i] += moleFraction(m_kTargets[k]) *
                                          linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void WeaklyIonizedGas::calculateTotalElasticCrossSection()
{
    m_sigmaElastic.clear();
    m_sigmaElastic.resize(m_points, 0.0);
    for (size_t k : m_kElastic) {
        vector_fp& x = m_crossSections[k][0];
        vector_fp& y = m_crossSections[k][1];
        for (size_t i = 0; i < m_points; i++) {
            double mass_ratio = ElectronMass / (Dalton * molecularWeight(m_kTargets[k]));
            m_sigmaElastic[i] += 2.0 * mass_ratio * moleFraction(m_kTargets[k]) *
                                 linearInterp(m_gridEdge[i], x, y);
        }
    }
}

void WeaklyIonizedGas::calculateDistributionFunction()
{
    // Check density and update density-dependent properties
    double N = molarDensity() * Avogadro;
    if (m_N != N) {
        m_N = N;
        m_f0_ok = false;
    }

    // Skip the calculation if the independent variables did not change.
    // Note gas composition is automatically checked with compositionChanged().
    if (m_f0_ok == true) {
        return;
    }

    checkSpeciesNoCrossSection();
    calculateTotalCrossSection();
    calculateTotalElasticCrossSection();

    double kTe = m_kT;
    //! Use kTe for initial f0
    if (m_init_kTe != 0.0) {
        kTe = m_init_kTe;
    }

    for (size_t j = 0; j < m_points; j++) {
        m_f0(j) = 2.0 * pow(1.0/Pi, 0.5) * pow(kTe, -3./2.) *
                  std::exp(-m_gridCenter[j]/kTe);
    }

    if (m_E != 0.0) {
        m_f0 = converge(m_f0);
        double Te = electronTemperature(m_f0);
        // Evaluate the EEDF by comparing electron temperature to gas temperature,
        // and replace the EEDF with a Maxwellian distribution at gas temperature.
        if (Te < temperature()) {
            for (size_t j = 0; j < m_points; j++) {
                m_f0(j) = 2.0 * pow(1.0/Pi, 0.5) * pow(m_kT, -3./2.) *
                          std::exp(-m_gridCenter[j]/m_kT);
            }
        }
    }
    m_f0_ok = true;
}

double WeaklyIonizedGas::integralPQ(double a, double b, double u0, double u1,
                                      double g, double x0)
{
    double A1;
    double A2;
    if (g != 0.0) {
        double expm1a = std::expm1(g * (-a + x0));
        double expm1b = std::expm1(g * (-b + x0));
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

vector_fp WeaklyIonizedGas::vector_g(Eigen::VectorXd& f0)
{
    vector_fp g(m_points, 0.0);
    g[0] = std::log(f0(1)/f0(0)) / (m_gridCenter[1] - m_gridCenter[0]);
    double N = m_points - 1;
    g[N] = std::log(f0(N)/f0(N-1)) / (m_gridCenter[N] - m_gridCenter[N-1]);
    for (size_t i = 1; i < m_points - 1; i++) {
        g[i] = std::log(f0(i+1)/f0(i-1)) / (m_gridCenter[i+1] - m_gridCenter[i-1]);
    }
    return g;
}

SparseMat WeaklyIonizedGas::matrix_P(vector_fp& g, size_t k)
{
    std::vector<Triplet> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double p = m_gamma * r;
        tripletList.push_back(Triplet(j, j, p));
    }
    SparseMat P(m_points, m_points);
    P.setFromTriplets(tripletList.begin(), tripletList.end());
    return P;
}

SparseMat WeaklyIonizedGas::matrix_Q(vector_fp& g, size_t k)
{
    std::vector<Triplet> tripletList;
    for (size_t n = 0; n < m_eps[k].size(); n++) {
        double eps_a = m_eps[k][n][0];
        double eps_b = m_eps[k][n][1];
        double sigma_a = m_sigma[k][n][0];
        double sigma_b = m_sigma[k][n][1];
        double i = m_i[k][n];
        double j = m_j[k][n];
        double r = integralPQ(eps_a, eps_b, sigma_a, sigma_b, g[j], m_gridCenter[j]);
        double q = m_inFactor[k] * m_gamma * r;
        tripletList.push_back(Triplet(i, j, q));
    }
    SparseMat Q(m_points, m_points);
    Q.setFromTriplets(tripletList.begin(), tripletList.end());
    return Q;
}

SparseMat WeaklyIonizedGas::matrix_A(Eigen::VectorXd& f0)
{
    vector_fp a0(m_points + 1);
    vector_fp a1(m_points + 1);
    size_t N = m_points - 1;
    // Scharfetter-Gummel scheme
    double nu = netProductionFreq(f0);
    a0[0] = NAN;
    a1[0] = NAN;
    a0[N+1] = NAN;
    a1[N+1] = NAN;
    for (size_t j = 1; j < m_points; j++) {
        double sigma_tilde = m_totalCrossSectionEdge[j] + nu / pow(m_gridEdge[j], 0.5) / m_gamma;
        double omega = 2 * Pi * m_F;
        double q = omega / (m_N * m_gamma * pow(m_gridEdge[j], 0.5));
        double W = -m_gamma * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double F = sigma_tilde * sigma_tilde / (sigma_tilde * sigma_tilde + q * q);
        double DA = m_gamma / 3.0 * pow(m_E / m_N, 2.0) * m_gridEdge[j];
        double DB = m_gamma * m_kT * m_gridEdge[j] * m_gridEdge[j] * m_sigmaElastic[j];
        double D = DA / sigma_tilde * F + DB;
        double z = W * (m_gridCenter[j] - m_gridCenter[j-1]) / D;
        a0[j] = W / (1 - std::exp(-z));
        a1[j] = W / (1 - std::exp(z));
    }

    std::vector<Triplet> tripletList;
    // center diagonal
    // zero flux b.c. at energy = 0
    tripletList.push_back(Triplet(0, 0, a0[1]));
    for (size_t j = 1; j < m_points - 1; j++) {
        tripletList.push_back(Triplet(j, j, a0[j+1] - a1[j]));
    }

    // upper diagonal
    for (size_t j = 0; j < m_points - 1; j++) {
        tripletList.push_back(Triplet(j, j+1, a1[j+1]));
    }

    // lower diagonal
    for (size_t j = 1; j < m_points; j++) {
        tripletList.push_back(Triplet(j, j-1, -a0[j]));
    }

    // zero flux b.c.
    tripletList.push_back(Triplet(N, N, -a1[N]));

    SparseMat A(m_points, m_points);
    A.setFromTriplets(tripletList.begin(), tripletList.end());

    //plus G
    SparseMat G(m_points, m_points);
    for (size_t i = 0; i < m_points; i++) {
        G.insert(i,i) = 2.0 / 3.0 * (pow(m_gridEdge[i+1], 1.5) - pow(m_gridEdge[i], 1.5)) * nu;
    }
    return A + G;
}

Eigen::VectorXd WeaklyIonizedGas::iterate(Eigen::VectorXd& f0, double delta)
{
    SparseMat PQ(m_points, m_points);
    vector_fp g = vector_g(f0);
    for (size_t k : m_kInelastic) {
        PQ += (matrix_Q(g, k) - matrix_P(g, k)) * moleFraction(m_kTargets[k]);
    }

    SparseMat A = matrix_A(f0);
    SparseMat I(m_points, m_points);
    for (size_t i = 0; i < m_points; i++) {
        I.insert(i,i) = 1.0;
    }
    A -= PQ;
    A *= delta;
    A += I;

    // add chemionization scattering-in rate at the first grid
    f0(0) += m_chemionScatRate;

    // solve f0
    Eigen::SparseLU<SparseMat> solver(A);
    Eigen::VectorXd f1 = solver.solve(f0);

    f1 *= 1.0 / norm(f1, m_gridCenter);
    return f1;
}

Eigen::VectorXd WeaklyIonizedGas::converge(Eigen::VectorXd& f0)
{
    double err0 = 0.0;
    double err1 = 0.0;
    double delta = m_delta0;
    for (size_t n = 0; n < m_maxn; n++) {
        if (0.0 < err1 && err1 < err0) {
            // log extrapolation attempting to reduce the error a factor m
            delta *= std::log(m_factorM) / (std::log(err0) - std::log(err1));
        }
        Eigen::VectorXd f1 = iterate(f0, delta);
        err0 = err1;
        Eigen::VectorXd Df0(m_points);
        for (size_t i = 0; i < m_points; i++) {
            Df0(i) = std::abs(f0(i) - f1(i));
        }
        err1 = norm(Df0, m_gridCenter);
        if (err1 < m_rtol) {
            return f1;
        }
        f0 = f1;
    }
    throw CanteraError("WeaklyIonizedGas::converge", "Convergence failed");
}

double WeaklyIonizedGas::netProductionFreq(Eigen::VectorXd& f0)
{
    double nu = m_chemionScatRate;
    vector_fp g = vector_g(f0);

    for (size_t k = 0; k < m_ncs; k++) {
        if (kind(k) == "ionization" ||
            kind(k) == "attachment") {
            SparseMat PQ = (matrix_Q(g, k) - matrix_P(g, k)) *
                           moleFraction(m_kTargets[k]);
            Eigen::VectorXd s = PQ * f0;
            for (size_t i = 0; i < m_points; i++) {
                nu += s[i];
            }
        }
    }
    return nu;
}

double WeaklyIonizedGas::electronDiffusivity()
{
    calculateDistributionFunction();
    vector_fp y(m_points, 0.0);
    double nu = netProductionFreq(m_f0);
    for (size_t i = 0; i < m_points; i++) {
        if (m_gridCenter[i] != 0.0) {
            y[i] = m_gridCenter[i] * m_f0(i) /
                   (m_totalCrossSectionCenter[i] + nu / m_gamma / pow(m_gridCenter[i], 0.5));
        }
    }
    return 1./3. * m_gamma * simpsonQuadrature(m_gridCenter, y) / m_N;
}

double WeaklyIonizedGas::electronMobility()
{
    calculateDistributionFunction();
    double nu = netProductionFreq(m_f0);
    vector_fp y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        if (m_gridEdge[i] != 0.0) {
            y[i] = m_gridEdge[i] * df0 /
                   (m_totalCrossSectionEdge[i] + nu / m_gamma / pow(m_gridEdge[i], 0.5));
        }
    }
    return -1./3. * m_gamma * simpsonQuadrature(m_gridEdge, y) / m_N;
}

double WeaklyIonizedGas::realMobility()
{
    calculateDistributionFunction();
    double nu = netProductionFreq(m_f0);
    vector_fp y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        if (m_gridEdge[i] != 0.0) {
            double Q = m_totalCrossSectionEdge[i] + nu / m_gamma / pow(m_gridEdge[i], 0.5);
            double q = 2.0 * Pi * m_F / (m_N * m_gamma * pow(m_gridEdge[i], 0.5));
            y[i] = m_gridEdge[i] * Q / (Q * Q + q * q) * df0;
        }
    }
    return -1./3. * m_gamma * simpsonQuadrature(m_gridEdge, y) / m_N;
}

double WeaklyIonizedGas::powerGain()
{
    if (m_F != 0.0) {
        return m_E * m_E * electronMobility();
    } else {
        return m_E * m_E * realMobility();
    }
}

double WeaklyIonizedGas::elasticPowerLoss()
{
    calculateDistributionFunction();
    double sum = 0.0;
    vector_fp y(m_points + 1, 0.0);
    for (size_t i = 1; i < m_points; i++) {
        // calculate df0 at i-1/2
        double df0 = (m_f0(i) - m_f0(i-1)) / (m_gridCenter[i] - m_gridCenter[i-1]);
        double f0 = 0.5 * (m_f0(i-1) + m_f0(i));
        y[i] = m_sigmaElastic[i] * (m_gridEdge[i] * m_gridEdge[i] * f0 + m_kT *  df0);
    }
    sum += m_gamma * simpsonQuadrature(m_gridEdge, y);
    return sum * m_N;
}

double WeaklyIonizedGas::totalCollisionFreq()
{
    double sum = 0.0;
    for (size_t k = 0; k < m_ncs; k++) {
        sum += moleFraction(m_kTargets[k]) * rateCoefficient(k);
    }
    return sum * m_N;
}

double WeaklyIonizedGas::biMaxwellFraction(size_t k)
{
    return 1.0 / (1 + exp(-threshold(k) / m_kT));
}

double WeaklyIonizedGas::inelasticPowerLoss()
{
    calculateDistributionFunction();
    double sum = 0.0;
    for (size_t k : m_kInelastic) {
        if (kind(k) == "excitation") {
            double y_low = biMaxwellFraction(k);
            double y_up = 1.0 - y_low;
            sum += threshold(k) * moleFraction(m_kTargets[k]) *
                   (y_low * rateCoefficient(k) -
                    y_up * reverseRateCoefficient(k));
        }
    }
    return sum * m_N;
}

double WeaklyIonizedGas::rateCoefficient(size_t k)
{
    calculateDistributionFunction();
    vector_fp g = vector_g(m_f0);
    SparseMat P = matrix_P(g, k);
    Eigen::VectorXd s = P * m_f0;
    double sum = 0.0;
    for (size_t i = 0; i < m_points; i++) {
        sum += s(i);
    }
    return sum;
}

double WeaklyIonizedGas::reverseRateCoefficient(size_t k)
{
    calculateDistributionFunction();
    vector_fp g = vector_g(m_f0);
    SparseMat Q = matrix_Q(g, k);
    Eigen::VectorXd s = Q * m_f0;
    double sum = 0.0;
    for (size_t i = 0; i < m_points; i++) {
        sum += s(i);
    }
    return sum;
}

double WeaklyIonizedGas::meanElectronEnergy()
{
    calculateDistributionFunction();
    double sum = 0;
    for (size_t i = 0; i < m_points - 1; i++) {
        sum += (pow(m_gridEdge[i+1], 2.5) -  pow(m_gridEdge[i], 2.5)) * m_f0(i);
    }
    return 0.4 * sum;
}

double WeaklyIonizedGas::electronTemperature()
{
    double Te = 2./3. * meanElectronEnergy() / Boltzmann * ElectronCharge;
    if (Te < temperature()) {
        return temperature();
    } else {
        return Te;
    }
}

double WeaklyIonizedGas::electronTemperature(Eigen::VectorXd f0)
{
    double sum = 0;
    for (size_t i = 0; i < m_points - 1; i++) {
        sum += (pow(m_gridEdge[i+1], 2.5) -  pow(m_gridEdge[i], 2.5)) * f0(i);
    }
    return 2./3. * 0.4 * sum / Boltzmann * ElectronCharge;
}

}
