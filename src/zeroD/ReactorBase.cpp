/**
 *  @file ReactorBase.cpp
 */

// Copyright 2001  California Institute of Technology

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorNet.h"

using namespace std;
namespace Cantera
{

ReactorBase::ReactorBase(const string& name) : m_nsp(0),
    m_thermo(0),
    m_vol(1.0),
    m_enthalpy(0.0),
    m_intEnergy(0.0),
    m_pressure(0.0),
    m_net(0)
{
    m_name = name;
}

void ReactorBase::setThermoMgr(thermo_t& thermo)
{
    m_thermo = &thermo;
    m_nsp = m_thermo->nSpecies();
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_pressure = m_thermo->pressure();
}

void ReactorBase::syncState()
{
    m_thermo->saveState(m_state);
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_pressure = m_thermo->pressure();
    if (m_net) {
        m_net->setNeedsReinit();
    }
}

void ReactorBase::addInlet(FlowDevice& inlet)
{
    m_inlet.push_back(&inlet);
}

void ReactorBase::addOutlet(FlowDevice& outlet)
{
    m_outlet.push_back(&outlet);
}

void ReactorBase::addWall(Wall& w, int lr)
{
    m_wall.push_back(&w);
    if (lr == 0) {
        m_lr.push_back(0);
    } else {
        m_lr.push_back(1);
    }
}

Wall& ReactorBase::wall(size_t n)
{
    return *m_wall[n];
}

ReactorNet& ReactorBase::network()
{
    if (m_net) {
        return *m_net;
    } else {
        throw CanteraError("ReactorBase::network",
                           "Reactor is not part of a ReactorNet");
    }
}

void ReactorBase::setNetwork(ReactorNet* net)
{
    m_net = net;
}

doublereal ReactorBase::residenceTime()
{
    doublereal mout = 0.0;
    for (size_t i = 0; i < m_outlet.size(); i++) {
        mout += m_outlet[i]->massFlowRate();
    }
    return mass()/mout;
}

FlowDevice& ReactorBase::inlet(size_t n)
{
    return *m_inlet[n];
}
FlowDevice& ReactorBase::outlet(size_t n)
{
    return *m_outlet[n];
}

}
