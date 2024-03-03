//! @file Wall.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/stringUtils.h"
#include "cantera/numerics/Func1.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/zeroD/ReactorNet.h"

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


void Wall::buildReactorJacobian(ReactorBase* r,
    vector<Eigen::Triplet<double>>& jacVector)
{
    // get derivative of heat transfer for both reactors
    vector<Eigen::Triplet<double>> network;
    size_t nsp = r->phase()->thermo()->nSpecies();
    size_t sidx = r->speciesOffset();
    size_t eidx = r->energyIndex();
    // define a scalar for direction based on left and right
    double direction = (r == m_left) ? 1.0 : -1.0;
    // elements within the current reactor
    // find dQdni for the current reactor w.r.t current reactor
    for (size_t i = sidx; i < nsp + sidx; i++) {
        double dQdni = m_rrth * m_area * direction * r->temperature_ddni(i);
        dQdni += m_emiss * m_area * direction * r->temperature_ddni(i) * 4
                 * pow(r->temperature(), 3);
        jacVector.emplace_back(eidx, i, dQdni);
    }
}

void Wall::buildNetworkJacobian(vector<Eigen::Triplet<double>>& jacVector)
{
    // No interdependent terms for reservoirs
    if (m_right->type() == "Reservoir" || m_left->type() == "Reservoir") {
        return;
    }
    // get derivatives for inter-dependent reactor terms
    //variables for the right side
    vector<Eigen::Triplet<double>> network;
    size_t r_nsp = m_right->phase()->thermo()->nSpecies();
    size_t r_sidx = m_right->speciesOffset();
    size_t r_net = m_right->network().globalStartIndex(m_right);
    size_t r_eidx = m_right->energyIndex();

    // variables for the left side
    size_t l_nsp = m_left->phase()->thermo()->nSpecies();
    size_t l_sidx = m_left->speciesOffset();
    size_t l_net = m_left->network().globalStartIndex(m_left);
    size_t l_eidx = m_left->energyIndex();

    if (m_right->energyEnabled()) {
        // find dQdni for the right reactor w.r.t left reactor
        for (size_t i = l_sidx; i < l_sidx + l_nsp; i++) {
            double dQdni = m_rrth * m_area * m_left->temperature_ddni(i);
            dQdni += m_emiss * m_area * m_left->temperature_ddni(i) * 4
                 * pow(m_left->temperature(), 3);
            jacVector.emplace_back(r_eidx + r_net, i + l_net, dQdni);
        }
    }

    if (m_left->energyEnabled()) {
        // find dQdni for the left reactor w.r.t right reactor
        for (size_t i = r_sidx; i < r_sidx + r_nsp; i++) {
            double dQdni = - m_rrth * m_area * m_right->temperature_ddni(i);
            dQdni -= m_emiss * m_area * m_right->temperature_ddni(i) * 4
                 * pow(m_right->temperature(), 3);
            jacVector.emplace_back(l_eidx + l_net, i + r_net, dQdni);
        }
    }
}

}
