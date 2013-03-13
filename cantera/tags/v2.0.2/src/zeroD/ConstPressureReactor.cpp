/**
 *  @file Reactor.cpp
 *
 *  A zero-dimensional reactor
 */

// Copyright 2001  California Institute of Technology

#include "cantera/zeroD/ConstPressureReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"

using namespace std;

namespace Cantera
{

ConstPressureReactor::ConstPressureReactor() : Reactor() {}

void ConstPressureReactor::
getInitialConditions(double t0, size_t leny, double* y)
{
    m_init = true;
    if (m_thermo == 0) {
        throw CanteraError("getInitialConditions",
                           "Error: reactor is empty.");
    }
    m_time = t0;
    m_thermo->restoreState(m_state);

    // total mass
    doublereal mass = m_thermo->density() * m_vol;

    // set components y + 2 ... y + K + 1 to the
    // mass M_k of each species
    m_thermo->getMassFractions(y+2);
    scale(y + 2, y + m_nsp + 2, y + 2, mass);

    // set the first component to the total enthalpy
    y[0] = m_thermo->enthalpy_mass() * mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // set the remaining components to the surface species
    // coverages on the walls
    size_t loc = m_nsp + 2;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            m_wall[m]->getCoverages(m_lr[m], y + loc);
            loc += surf->nSpecies();
        }
    }
}

void ConstPressureReactor::initialize(doublereal t0)
{
    m_thermo->restoreState(m_state);
    m_sdot.resize(m_nsp, 0.0);
    m_nv = m_nsp + 2;
    for (size_t w = 0; w < m_nwalls; w++)
        if (m_wall[w]->surface(m_lr[w])) {
            m_nv += m_wall[w]->surface(m_lr[w])->nSpecies();
        }
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();

    size_t nt = 0, maxnt = 0;
    for (size_t m = 0; m < m_nwalls; m++) {
        if (m_wall[m]->kinetics(m_lr[m])) {
            nt = m_wall[m]->kinetics(m_lr[m])->nTotalSpecies();
            if (nt > maxnt) {
                maxnt = nt;
            }
            if (m_wall[m]->kinetics(m_lr[m])) {
                if (&m_kin->thermo(0) !=
                        &m_wall[m]->kinetics(m_lr[m])->thermo(0)) {
                    throw CanteraError("ConstPressureReactor::initialize",
                                       "First phase of all kinetics managers must be"
                                       " the gas.");
                }
            }
        }
    }
    m_work.resize(maxnt);
    m_init = true;
}

void ConstPressureReactor::updateState(doublereal* y)
{

    // The components of y are the total enthalpy,
    // the total volume, and the mass of each species.
    doublereal h   = y[0];
    doublereal* mss = y + 2;
    doublereal mass = accumulate(y+2, y+2+m_nsp, 0.0);
    m_thermo->setMassFractions(mss);

    if (m_energy) {
        m_thermo->setState_HP(h/mass, m_pressure, 1.0e-4);
    } else {
        m_thermo->setPressure(m_pressure);
    }
    m_vol = mass / m_thermo->density();

    size_t loc = m_nsp + 2;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            m_wall[m]->setCoverages(m_lr[m], y+loc);
            loc += surf->nSpecies();
        }
    }

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}


/*
 * Called by the integrator to evaluate ydot given y at time 'time'.
 */
void ConstPressureReactor::evalEqs(doublereal time, doublereal* y,
                                   doublereal* ydot, doublereal* params)
{
    size_t nk;
    m_time = time;
    m_thermo->restoreState(m_state);

    Kinetics* kin;
    size_t npar, ploc;
    double mult;

    // process sensitivity parameters
    if (params) {

        npar = m_pnum.size();
        for (size_t n = 0; n < npar; n++) {
            mult = m_kin->multiplier(m_pnum[n]);
            m_kin->setMultiplier(m_pnum[n], mult*params[n]);
        }
        ploc = npar;
        for (size_t m = 0; m < m_nwalls; m++) {
            if (m_nsens_wall[m] > 0) {
                m_wall[m]->setSensitivityParameters(m_lr[m], params + ploc);
                ploc += m_nsens_wall[m];
            }
        }
    }

    m_vdot = 0.0;
    m_Q    = 0.0;

    // compute wall terms
    doublereal rs0, sum, wallarea;

    SurfPhase* surf;
    size_t lr, ns, loc = m_nsp+2, surfloc;
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    for (size_t i = 0; i < m_nwalls; i++) {
        lr = 1 - 2*m_lr[i];
        m_Q += lr*m_wall[i]->Q(time);
        kin = m_wall[i]->kinetics(m_lr[i]);
        surf = m_wall[i]->surface(m_lr[i]);
        if (surf && kin) {
            rs0 = 1.0/surf->siteDensity();
            nk = surf->nSpecies();
            sum = 0.0;
            surf->setTemperature(m_state[0]);
            m_wall[i]->syncCoverages(m_lr[i]);
            kin->getNetProductionRates(DATA_PTR(m_work));
            ns = kin->surfacePhaseIndex();
            surfloc = kin->kineticsSpeciesIndex(0,ns);
            for (size_t k = 1; k < nk; k++) {
                ydot[loc + k] = m_work[surfloc+k]*rs0*surf->size(k);
                sum -= ydot[loc + k];
            }
            ydot[loc] = sum;
            loc += nk;

            wallarea = m_wall[i]->area();
            for (size_t k = 0; k < m_nsp; k++) {
                m_sdot[k] += m_work[k]*wallarea;
            }
        }
    }

    // dummy equation
    ydot[1] = 0.0;

    /* species equations
     *  Equation is:
     *  \dot M_k = \hat W_k \dot\omega_k + \dot m_{in} Y_{k,in}
     *             - \dot m_{out} Y_{k} + A \dot s_k.
     */
    const doublereal* mw = DATA_PTR(m_thermo->molecularWeights());
    if (m_chem) {
        m_kin->getNetProductionRates(ydot+2);   // "omega dot"
    } else {
        fill(ydot + 2, ydot + 2 + m_nsp, 0.0);
    }
    for (size_t n = 0; n < m_nsp; n++) {
        ydot[n+2] *= m_vol;     //           moles/s/m^3 -> moles/s
        ydot[n+2] += m_sdot[n];
        ydot[n+2] *= mw[n];
    }


    /*
     *  Energy equation.
     *  \f[
     *  \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in}
     * - \dot m_{out} h.
     * \f]
     */
    if (m_energy) {
        ydot[0] = - m_Q;
    } else {
        ydot[0] = 0.0;
    }

    // add terms for open system
    if (m_open) {

        const doublereal* mf = m_thermo->massFractions();
        doublereal enthalpy = m_thermo->enthalpy_mass();

        // outlets

        doublereal mdot_out;
        for (size_t i = 0; i < m_nOutlets; i++) {
            mdot_out = m_outlet[i]->massFlowRate(time);
            for (size_t n = 0; n < m_nsp; n++) {
                ydot[2+n] -= mdot_out * mf[n];
            }
            if (m_energy) {
                ydot[0] -= mdot_out * enthalpy;
            }
        }


        // inlets

        doublereal mdot_in;
        for (size_t i = 0; i < m_nInlets; i++) {
            mdot_in = m_inlet[i]->massFlowRate(time);
            for (size_t n = 0; n < m_nsp; n++) {
                ydot[2+n] += m_inlet[i]->outletSpeciesMassFlowRate(n);
            }
            if (m_energy) {
                ydot[0] += mdot_in * m_inlet[i]->enthalpy_mass();
            }
        }
    }

    // reset sensitivity parameters
    if (params) {
        npar = m_pnum.size();
        for (size_t n = 0; n < npar; n++) {
            mult = m_kin->multiplier(m_pnum[n]);
            m_kin->setMultiplier(m_pnum[n], mult/params[n]);
        }
        ploc = npar;
        for (size_t m = 0; m < m_nwalls; m++) {
            if (m_nsens_wall[m] > 0) {
                m_wall[m]->resetSensitivityParameters(m_lr[m]);
                ploc += m_nsens_wall[m];
            }
        }
    }
}

size_t ConstPressureReactor::componentIndex(string nm) const
{
    if (nm == "H") {
        return 0;
    }
    if (nm == "V") {
        return 1;
    }
    // check for a gas species name
    size_t k = m_thermo->speciesIndex(nm);
    if (k != npos) {
        return k + 2;
    }

    // check for a wall species
    size_t walloffset = 0, kp = 0;
    thermo_t* th;
    for (size_t m = 0; m < m_nwalls; m++) {
        if (m_wall[m]->kinetics(m_lr[m])) {
            kp = m_wall[m]->kinetics(m_lr[m])->reactionPhaseIndex();
            th = &m_wall[m]->kinetics(m_lr[m])->thermo(kp);
            k = th->speciesIndex(nm);
            if (k != npos) {
                return k + 2 + m_nsp + walloffset;
            } else {
                walloffset += th->nSpecies();
            }
        }
    }
    return npos;
}

}
