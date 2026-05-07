//! @file Wall.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/numerics/Func1.h"
#include "cantera/zeroD/Wall.h"

namespace Cantera
{

WallBase::WallBase(shared_ptr<ReactorBase> r0, shared_ptr<ReactorBase> r1,
                   const string& name) : ConnectorNode(r0, r1, name)
{
    if (!m_nodes.first || !m_nodes.second) {
        throw CanteraError("WallBase::WallBase",
            "Reactors must be provided to WallBase constructor.");
    }
    m_left = r0.get();
    m_right = r1.get();
    m_left->addWall(*this, 0);
    m_right->addWall(*this, 1);
}

void WallBase::setArea(double a) {
    m_area = a;
}

double Wall::velocity() const {
    if (m_vf) {
        return m_vf->eval(m_time);
    }
    return 0.;
}

double Wall::expansionRate()
{
    if (!ready()) {
        throw CanteraError("Wall::expansionRate",
                           "Wall is not ready; some parameters have not been set.");
    }
    double rate = m_k * m_area * (m_left->pressure() - m_right->pressure());

    if (m_vf) {
        rate += m_area * m_vf->eval(m_time);
    }
    return rate;
}

double Wall::heatFlux() const {
    if (m_qf) {
        return m_qf->eval(m_time);
    }
    return 0.;
}

double Wall::heatRate()
{
    if (!ready()) {
        throw CanteraError("Wall::heatRate",
                           "Wall is not ready; some parameters have not been set.");
    }
    double q1 = (m_area * m_rrth) *
                (m_left->temperature() - m_right->temperature());
    if (m_emiss > 0.0) {
        double tl = m_left->temperature();
        double tr = m_right->temperature();
        q1 += m_emiss * m_area * StefanBoltz * (tl*tl*tl*tl - tr*tr*tr*tr);
    }

    if (m_qf) {
        q1 += m_area * m_qf->eval(m_time);
    }
    return q1;
}

void Wall::addExpansionRateJacobian(
    SparseTriplets& trips, size_t row, double coeff, bool includePressureSpecies)
{
    double alpha = m_k * m_area;
    if (alpha == 0.0) {
        return;
    }
    // Chain rule for expansionRate(P_left - P_right): the wall supplies
    // d(expansionRate)/dDeltaP, while reactors supply dP/dy.
    m_left->addPressureJacobian(trips, row, coeff * alpha, includePressureSpecies);
    m_right->addPressureJacobian(trips, row, -coeff * alpha, includePressureSpecies);
}

void Wall::addHeatRateJacobian(SparseTriplets& trips, size_t row, double coeff)
{
    double leftCoeff = m_area * m_rrth;
    double rightCoeff = -leftCoeff;
    if (m_emiss > 0.0) {
        leftCoeff += 4.0 * m_emiss * m_area * StefanBoltz
                     * std::pow(m_left->temperature(), 3);
        rightCoeff -= 4.0 * m_emiss * m_area * StefanBoltz
                      * std::pow(m_right->temperature(), 3);
    }
    // Chain rule for heatRate(T_left, T_right): the wall supplies dQ/dT, while
    // reactors supply dT/dy.
    if (leftCoeff != 0.0) {
        m_left->addTemperatureJacobian(trips, row, coeff * leftCoeff);
    }
    if (rightCoeff != 0.0) {
        m_right->addTemperatureJacobian(trips, row, coeff * rightCoeff);
    }
}

}
