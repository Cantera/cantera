//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include <boost/math/special_functions/gamma.hpp>
#include "cantera/base/global.h"
#include "cantera/numerics/funcs.h"

namespace Cantera {

PlasmaPhase::PlasmaPhase(const std::string& inputFile, const std::string& id_)
    : m_x(2.0)
    , m_nPoints(1000)
    , m_electronName("E")
{
    initThermoFile(inputFile, id_);

    // initial grid
    m_electronEnergyLevels = Eigen::VectorXd::LinSpaced(m_nPoints, 0.001, 1.0);

    // initial electron temperature
    setElectronTemperature(temperature());
}

void PlasmaPhase::setIsotropicElectronEnergyDistribution()
{
    m_electronEnergyDist.resize(m_nPoints);
    double gamma1 = boost::math::tgamma(3.0 / 2.0 * m_x);
    double gamma2 = boost::math::tgamma(5.0 / 2.0 * m_x);
    double c1 = m_x * std::pow(gamma2, 1.5) / std::pow(gamma1, 2.5);
    double c2 = m_x * std::pow(gamma2 / gamma1, m_x);
    m_electronEnergyDist =
        c1 * m_electronEnergyLevels.array().sqrt() /
        std::pow(m_meanElectronEnergy, 1.5) *
        (-c2 * (m_electronEnergyLevels.array() /
        m_meanElectronEnergy).pow(m_x)).exp();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    m_electronTemp = Te;
    m_meanElectronEnergy = 3.0 / 2.0 * electronTemperature() *
                           Boltzmann / ElectronCharge;
    setIsotropicElectronEnergyDistribution();
}

void PlasmaPhase::setElectronEnergyLevels(const vector_fp& levels)
{
    m_nPoints = levels.size();
    m_electronEnergyLevels =
        Eigen::Map<const Eigen::VectorXd>(levels.data(), levels.size());
    setIsotropicElectronEnergyDistribution();
}

void PlasmaPhase::getElectronEnergyLevels(vector_fp& levels) const
{
    levels.resize(m_nPoints);
    Eigen::Map<Eigen::VectorXd>(levels.data(), m_nPoints) = m_electronEnergyLevels;
}

void PlasmaPhase::setElectronEnergyDistribution(const vector_fp& levels,
                                                const vector_fp& distrb)
{
    if (levels.size() != distrb.size()) {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistribution",
                           "Vector lengths need to be the same.");
    }
    setElectronEnergyLevels(levels);
    m_nPoints = levels.size();
    m_electronEnergyDistrb =
        Eigen::Map<const Eigen::VectorXd>(distrb.data(), distrb.size());
    // calculate mean electron energy and electron temperature
    Eigen::VectorXd eps52 = m_electronEnergyLevels.array().pow(5./2.);
    m_meanElectronEnergy =
        2.0 / 5.0 * trapezoidal(m_electronEnergyDist, eps52);
    Phase::setElectronTemperature(
        2.0 / 3.0 * m_meanElectronEnergy * Avogadro *
        ElectronCharge / GasConstant);
}

void PlasmaPhase::getElectronEnergyDistribution(vector_fp& distrb) const
{
    distrb.resize(m_nPoints);
    Eigen::Map<Eigen::VectorXd>(distrb.data(), m_nPoints) = m_electronEnergyDistrb;
}

void PlasmaPhase::updateThermo() const
{
    IdealGasPhase::updateThermo();
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);
    double tempNow = temperature();
    double electronTempNow = electronTemperature();
    size_t k = speciesIndex(m_electronName);
    // If the electron temperature has changed since the last time these
    // properties were computed, recompute them.
    if (cached.state1 != tempNow && cached.state2 != electronTempNow) {
        m_spthermo.update_single(k, electronTemperature(),
                &m_cp0_R[k], &m_h0_RT[k], &m_s0_R[k]);
        cached.state1 = tempNow;
        cached.state2 = electronTempNow;
    }
    // update the species Gibbs functions
    m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
}

}
