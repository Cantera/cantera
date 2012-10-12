/**
 *  @file Reactor.cpp
 *
 *  A zero-dimensional reactor
 */

// Copyright 2001  California Institute of Technology

#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"

using namespace std;

namespace Cantera
{
Reactor::Reactor() : ReactorBase(),
    m_kin(0),
    m_vdot(0.0),
    m_Q(0.0),
    m_chem(true),
    m_energy(true),
    m_nsens(npos)
{}

// overloaded method of FuncEval. Called by the integrator to
// get the initial conditions.
void Reactor::getInitialConditions(double t0, size_t leny, double* y)
{
    m_init = true;
    if (m_thermo == 0) {
        cout << "Error: reactor is empty." << endl;
        return;
    }
    m_thermo->restoreState(m_state);

    // total mass
    doublereal mass = m_thermo->density() * m_vol;

    // set components y + 2 ... y + K + 1 to the
    // mass M_k of each species
    m_thermo->getMassFractions(y+2);
    scale(y + 2, y + m_nsp + 2, y + 2, mass);

    // set the first component to the total internal
    // energy
    y[0] = m_thermo->intEnergy_mass() * mass;

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
            //surf->getCoverages(y+loc);
            loc += surf->nSpecies();
        }
    }
}

/*
 *  Must be called before calling method 'advance'
 */
void Reactor::initialize(doublereal t0)
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
                    throw CanteraError("Reactor::initialize",
                                       "First phase of all kinetics managers must be"
                                       " the gas.");
                }
            }
        }
    }
    m_work.resize(maxnt);
    m_init = true;
}

size_t Reactor::nSensParams()
{
    if (m_nsens == npos) {
        // determine the number of sensitivity parameters
        size_t m, ns;
        m_nsens = m_pnum.size();
        for (m = 0; m < m_nwalls; m++) {
            ns = m_wall[m]->nSensParams(m_lr[m]);
            m_nsens_wall.push_back(ns);
            m_nsens += ns;
        }
    }
    return m_nsens;
}

void Reactor::updateState(doublereal* y)
{
    ThermoPhase& mix = *m_thermo;  // define for readability

    // The components of y are the total internal energy,
    // the total volume, and the mass of each species.

    // Set the mass fractions and  density of the mixture.


    doublereal u   = y[0];
    m_vol          = y[1];
    doublereal* mss = y + 2;
    doublereal mass = accumulate(y+2, y+2+m_nsp, 0.0);
    m_thermo->setMassFractions(mss);

    m_thermo->setDensity(mass/m_vol);

    doublereal temp = temperature();
    mix.setTemperature(temp);

    if (m_energy) {
        // Decreased the tolerance on delta_T to 1.0E-7 so that T is
        // accurate to 9 sig digits, because this is
        // used in the numerical jacobian routines where relative values
        // of 1.0E-7 are used in the deltas.
        m_thermo->setState_UV(u/mass,m_vol/mass, 1.0e-7);
        temp = mix.temperature(); //mix.setTemperature(temp);
    }
    //m_state[0] = temp;

    size_t loc = m_nsp + 2;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            //                surf->setTemperature(temp);
            //surf->setCoverages(y+loc);
            m_wall[m]->setCoverages(m_lr[m], y+loc);
            loc += surf->nSpecies();
        }
    }

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}


/*
 * Called by the integrator to evaluate ydot given y at time 'time'.
 */
void Reactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    m_thermo->restoreState(m_state);

    // process sensitivity parameters
    if (params) {
        size_t npar = m_pnum.size();
        for (size_t n = 0; n < npar; n++) {
            double mult = m_kin->multiplier(m_pnum[n]);
            m_kin->setMultiplier(m_pnum[n], mult*params[n]);
        }
        size_t ploc = npar;
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
    size_t loc = m_nsp+2;
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    for (size_t i = 0; i < m_nwalls; i++) {
        int lr = 1 - 2*m_lr[i];
        double vdot = lr*m_wall[i]->vdot(time);
        m_vdot += vdot;
        m_Q += lr*m_wall[i]->Q(time);
        Kinetics* kin = m_wall[i]->kinetics(m_lr[i]);
        SurfPhase* surf = m_wall[i]->surface(m_lr[i]);
        if (surf && kin) {
            double rs0 = 1.0/surf->siteDensity();
            size_t nk = surf->nSpecies();
            double sum = 0.0;
            surf->setTemperature(m_state[0]);
            m_wall[i]->syncCoverages(m_lr[i]);
            kin->getNetProductionRates(DATA_PTR(m_work));
            size_t ns = kin->surfacePhaseIndex();
            size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
            for (size_t k = 1; k < nk; k++) {
                ydot[loc + k] = m_work[surfloc+k]*rs0*surf->size(k);
                sum -= ydot[loc + k];
            }
            ydot[loc] = sum;
            loc += nk;

            double wallarea = m_wall[i]->area();
            for (size_t k = 0; k < m_nsp; k++) {
                m_sdot[k] += m_work[k]*wallarea;
            }
        }
    }

    // volume equation
    ydot[1] = m_vdot;

    /* species equations
     *  Equation is:
     *  \dot M_k = \hat W_k \dot\omega_k + \dot m_{in} Y_{k,in}
     *             - \dot m_{out} Y_{k} + A \dot s_k.
     */
    const vector_fp& mw = m_thermo->molecularWeights();
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
        ydot[0] = - m_thermo->pressure() * m_vdot - m_Q;
    } else {
        ydot[0] = 0.0;
    }

    // add terms for open system
    if (m_open) {
        const doublereal* mf = m_thermo->massFractions();
        doublereal enthalpy = m_thermo->enthalpy_mass();

        // outlets
        for (size_t i = 0; i < m_nOutlets; i++) {
            double mdot_out = m_outlet[i]->massFlowRate(time);
            for (size_t n = 0; n < m_nsp; n++) {
                ydot[2+n] -= mdot_out * mf[n];
            }
            if (m_energy) {
                ydot[0] -= mdot_out * enthalpy;
            }
        }

        // inlets
        for (size_t i = 0; i < m_nInlets; i++) {
            double mdot_in = m_inlet[i]->massFlowRate(time);
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
        size_t npar = m_pnum.size();
        for (size_t n = 0; n < npar; n++) {
            double mult = m_kin->multiplier(m_pnum[n]);
            m_kin->setMultiplier(m_pnum[n], mult/params[n]);
        }
        size_t ploc = npar;
        for (size_t m = 0; m < m_nwalls; m++) {
            if (m_nsens_wall[m] > 0) {
                m_wall[m]->resetSensitivityParameters(m_lr[m]);
                ploc += m_nsens_wall[m];
            }
        }
    }
}

void Reactor::addSensitivityReaction(size_t rxn)
{
    if (rxn >= m_kin->nReactions())
        throw CanteraError("Reactor::addSensitivityReaction",
                           "Reaction number out of range ("+int2str(rxn)+")");
    m_pnum.push_back(rxn);
    m_pname.push_back(name()+": "+m_kin->reactionString(rxn));
    m_mult_save.push_back(1.0);
}


size_t Reactor::componentIndex(const string& nm) const
{
    if (nm == "U") {
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
