//! @file PlasmaPhase.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PlasmaPhase.h"
#include <boost/math/special_functions/gamma.hpp>
#include "cantera/base/global.h"

namespace Cantera {

PlasmaPhase::PlasmaPhase(const std::string& inputFile, const std::string& id_)
    : m_electronEnergyDistrbMethod("isotropic-velocity")
    , m_x(2.0)
    , m_nPoints(1000)
{
    initThermoFile(inputFile, id_);

    // initial grid
    m_electronEnergyGrid.resize(m_nPoints);
    for (size_t i = 0; i <= m_nPoints - 1; i++) {
        m_electronEnergyGrid(i) = 0.001 + 0.001 * i;
    }
    // initial electron temperature
    setElectronTemperature(temperature());
}

void PlasmaPhase::updateIsotropicElectronEnergyDistrb()
{
    m_electronEnergyDistrb.resize(m_nPoints);
    m_meanElectronEnergy = 3.0 / 2.0 * electronTemperature() * GasConstant / (Avogadro * ElectronCharge);
    double gamma1 = boost::math::tgamma(3.0 / 2.0 * m_x);
    double gamma2 = boost::math::tgamma(5.0 / 2.0 * m_x);
    double c1 = m_x * std::pow(gamma2, 1.5) / std::pow(gamma1, 2.5);
    double c2 = m_x * std::pow(gamma2 / gamma1, m_x);
    m_electronEnergyDistrb = c1 * m_electronEnergyGrid.array().sqrt() / std::pow(m_meanElectronEnergy, 1.5) *
                      (-c2 * (m_electronEnergyGrid.array() / m_meanElectronEnergy).pow(m_x)).exp();
}

void PlasmaPhase::setElectronTemperature(const double Te) {
    Phase::setElectronTemperature(Te);
    if (m_electronEnergyDistrbMethod == "isotropic-velocity") {
        updateIsotropicElectronEnergyDistrb();
    } else if (m_electronEnergyDistrbMethod == "user-specified") {
        warn_user("PlasmaPhase::setElectronTemperature",
            "The electron temperature is calculated automatically when"
            "electron energy distribution is user-specified.");
    }
}

void PlasmaPhase::setElectronEnergyDistrbMethod(std::string method)
{
    if (method == "isotropic-velocity" || method == "user-specified") {
        m_electronEnergyDistrbMethod = method;
    } else {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistrbMethod",
                           "Invalid method for electron energy distribution."
                           "Please use below method, 'isotropic-velocity' or"
                           "'user-specified'.");
    }
}

void PlasmaPhase::setElectronEnergyGrid(const vector_fp& grid)
{
    m_nPoints = grid.size();
    m_electronEnergyGrid.resize(m_nPoints);
    for (size_t i = 0; i < m_nPoints; i++) {
        m_electronEnergyGrid(i) = grid[i];
    }
}

void PlasmaPhase::getElectronEnergyGrid(vector_fp& grid) const
{
    grid.resize(m_nPoints);
    for (size_t i = 0; i < m_nPoints; i++) {
        grid[i] = m_electronEnergyGrid(i);
    }
}

void PlasmaPhase::setElectronEnergyDistrb(const vector_fp& grid,
                                          const vector_fp& distrb)
{
    if (grid.size() != distrb.size()) {
        throw CanteraError("PlasmaPhase::setElectronEnergyDistrb",
                           "Vector lengths need to be the same.");
    }
    setElectronEnergyGrid(grid);
    m_nPoints = grid.size();
    m_electronEnergyDistrb.resize(m_nPoints);
    for (size_t i = 0; i < m_nPoints; i++) {
        m_electronEnergyDistrb(i) = distrb[i];
    }
}

void PlasmaPhase::getElectronEnergyDistrb(vector_fp& distrb) const
{
    distrb.resize(m_nPoints);
    for (size_t i = 0; i < m_nPoints; i++) {
        distrb[i] = m_electronEnergyDistrb(i);
    }
}

}
