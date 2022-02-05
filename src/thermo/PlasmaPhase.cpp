//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include <boost/math/special_functions/gamma.hpp>

namespace Cantera {

PlasmaPhase::PlasmaPhase(const std::string& inputFile, const std::string& id_)
    : m_electronEnergyDistrbMethod("isotropic-velocity")
    , m_x(2.0)
    , m_nPoints(1000)
{
    initThermoFile(inputFile, id_);

    // initial grid
    m_EE_grid.resize(m_nPoints);
    for (size_t i = 0; i <= m_nPoints - 1; i++) {
        m_EE_grid(i) = 0.001 + 0.001 * i;
    }
    // initial electron temperature
    setElectronTemperature(temperature());
}

void PlasmaPhase::updateIsotropicElectronEnergyDistrb()
{
    m_EE_distrb_vec.resize(m_nPoints);
    m_EE_mean = 3.0 / 2.0 * electronTemperature() * GasConstant / (Avogadro * ElectronCharge);
    double gamma1 = boost::math::tgamma(3.0 / 2.0 * m_x);
    double gamma2 = boost::math::tgamma(5.0 / 2.0 * m_x);
    double c1 = m_x * std::pow(gamma2, 1.5) / std::pow(gamma1, 2.5);
    double c2 = m_x * std::pow(gamma2 / gamma1, m_x);
    m_EE_distrb_vec = c1 * m_EE_grid.array().sqrt() / std::pow(m_EE_mean, 1.5) *
                      (-c2 * (m_EE_grid.array() / m_EE_mean).pow(m_x)).exp();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    Phase::setElectronTemperature(Te);
    if (m_electronEnergyDistrbMethod == "isotropic-velocity") {
        updateIsotropicElectronEnergyDistrb();
    }
}

void PlasmaPhase::setElectronEnergyGrid(const vector_fp& grid)
{
    m_nPoints = grid.size();
    m_EE_grid.resize(m_nPoints);
    for (size_t i = 0; i < m_nPoints; i++) {
        m_EE_grid(i) = grid[i];
    }
}

void PlasmaPhase::getElectronEnergyGrid(vector_fp& grid) const
{
    grid.resize(m_nPoints);
    for (size_t i = 0; i < m_nPoints; i++) {
        grid[i] = m_EE_grid(i);
    }
}

}
