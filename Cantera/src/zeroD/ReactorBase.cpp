/**
 *  @file ReactorBase.cpp
 *
 * $Author: dggoodwin $
 * $Revision: 1.11 $
 * $Date: 2006/11/27 21:43:34 $
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ReactorBase.h"
#include "FlowDevice.h"
#include "Wall.h"

using namespace std;
namespace CanteraZeroD {

    ReactorBase::ReactorBase(string name) : m_nsp(0), 
                                 m_thermo(0), 
                                 m_time(0.0), 
                                 m_vol(1.0), 
                                 m_vol0(1.0), 
                                 m_init(false), 
                                 m_nInlets(0), 
                                 m_nOutlets(0),
                                 m_open(false), 
                                 m_enthalpy(0.0),
                                 m_intEnergy(0.0), 
                                 m_pressure(0.0),
                                 m_nwalls(0)
    {
        m_name = name;
    }
 
//     void ReactorBase::resetState() {
//         m_thermo->saveState(m_state);
//         m_enthalpy = m_thermo->enthalpy_mass();
//         m_intEnergy = m_thermo->intEnergy_mass();
//         m_pressure = m_thermo->pressure();
//         m_init = false;
//     }

    void ReactorBase::setThermoMgr(thermo_t& thermo){
        m_thermo = &thermo;
        m_nsp = m_thermo->nSpecies();
        m_thermo->saveState(m_state);
        m_enthalpy = m_thermo->enthalpy_mass();
        m_intEnergy = m_thermo->intEnergy_mass();
        m_pressure = m_thermo->pressure();
    }

    void ReactorBase::addInlet(FlowDevice& inlet) {
        m_inlet.push_back(&inlet);
        m_open = true;
        m_nInlets++;
    }

    void ReactorBase::addOutlet(FlowDevice& outlet) {
        m_outlet.push_back(&outlet);
        m_open = true;
        m_nOutlets++;
    }

    void ReactorBase::addWall(Wall& w, int lr) {
        m_wall.push_back(&w);
        if (lr == 0) m_lr.push_back(0);
        else m_lr.push_back(1);
        m_nwalls++;
    }

    Wall& ReactorBase::wall(int n) {
        return *m_wall[n];
    }

    doublereal ReactorBase::residenceTime() {
        int nout = static_cast<int>(m_outlet.size());
        doublereal mout = 0.0;
        for (int i = 0; i < nout; i++) 
            mout += m_outlet[i]->massFlowRate();
        return mass()/mout;
    }

    FlowDevice& ReactorBase::inlet(int n)  { return *m_inlet[n]; }
    FlowDevice& ReactorBase::outlet(int n) { return *m_outlet[n]; }

}
