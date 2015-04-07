/**
 *  @file ReactorBase.cpp
 */

// Copyright 2001  California Institute of Technology

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"

using namespace std;
namespace Cantera
{

ReactorBase::ReactorBase(const string& name) : m_nsp(0),
    m_thermo(0),
    m_vol(1.0),
    m_vol0(1.0),
    m_init(false),
    m_nInlets(0),
    m_nOutlets(0),
    m_open(false),
    m_enthalpy(0.0),
    m_intEnergy(0.0),
    m_pressure(0.0),
    m_nwalls(0),
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

void ReactorBase::addInlet(FlowDevice& inlet)
{
    m_inlet.push_back(&inlet);
    m_open = true;
    m_nInlets++;
}

void ReactorBase::addOutlet(FlowDevice& outlet)
{
    m_outlet.push_back(&outlet);
    m_open = true;
    m_nOutlets++;
}

void ReactorBase::addWall(Wall& w, int lr)
{
    m_wall.push_back(&w);
    if (lr == 0) {
        m_lr.push_back(0);
    } else {
        m_lr.push_back(1);
    }
    m_nwalls++;
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
    int nout = static_cast<int>(m_outlet.size());
    doublereal mout = 0.0;
    for (int i = 0; i < nout; i++) {
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
