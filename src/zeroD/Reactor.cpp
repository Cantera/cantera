/**
 *  @file Reactor.cpp A zero-dimensional reactor
 */

// Copyright 2001  California Institute of Technology

#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/ReactorNet.h"

#include <cfloat>

using namespace std;

namespace Cantera
{
Reactor::Reactor() : ReactorBase(),
    m_kin(0),
    m_vdot(0.0),
    m_Q(0.0),
    m_chem(false),
    m_energy(true),
    m_nsens(npos)
{}

void Reactor::getInitialConditions(double t0, size_t leny, double* y)
{
    m_init = true;
    if (m_thermo == 0) {
        cout << "Error: reactor is empty." << endl;
        return;
    }
    m_thermo->restoreState(m_state);

    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // set the third component to the total internal energy
    y[2] = m_thermo->intEnergy_mass() * m_mass;

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+3);

    // set the remaining components to the surface species
    // coverages on the walls
    size_t loc = m_nsp + 3;
    SurfPhase* surf;
    for (size_t m = 0; m < m_nwalls; m++) {
        surf = m_wall[m]->surface(m_lr[m]);
        if (surf) {
            m_wall[m]->getCoverages(m_lr[m], y + loc);
            loc += surf->nSpecies();
        }
    }
}

void Reactor::initialize(doublereal t0)
{
    m_thermo->restoreState(m_state);
    m_sdot.resize(m_nsp, 0.0);
    m_wdot.resize(m_nsp, 0.0);
    m_nv = m_nsp + 3;
    for (size_t w = 0; w < m_nwalls; w++)
        if (m_wall[w]->surface(m_lr[w])) {
            m_nv += m_wall[w]->surface(m_lr[w])->nSpecies();
        }

    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();

    size_t nt = 0, maxnt = 0;
    for (size_t m = 0; m < m_nwalls; m++) {
        m_wall[m]->initialize();
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
    std::sort(m_pnum.begin(), m_pnum.end());
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
    for (size_t i = 0; i < m_nv; i++) {
        AssertFinite(y[i], "Reactor::updateState",
                     "y[" + int2str(i) + "] is not finite");
    }

    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the total internal energy, [3...K+3] are the mass fractions of each
    // species, and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];

    m_thermo->setMassFractions_NoNorm(y+3);

    if (m_energy) {
        // Use Newton's method to determine the mixture temperature. Tight
        // tolerances are required both for Jacobian evaluation and for
        // sensitivity analysis to work correctly.

        doublereal U = y[2];
        doublereal T = temperature();
        double dT = 100;

        int i = 0;
        while (abs(dT / T) > 10 * DBL_EPSILON) {
            m_thermo->setState_TR(T, m_mass / m_vol);
            double dUdT = m_thermo->cv_mass() * m_mass;
            dT = (m_thermo->intEnergy_mass() * m_mass - U) / dUdT;
            dT = std::min(dT, 0.5 * T);
            T -= dT;
            i++;
            if (i > 100) {
                std::string message = "no convergence";
                message += "\nU/m = " + fp2str(U / m_mass);
                message += "\nT = " + fp2str(T);
                message += "\nrho = " + fp2str(m_mass / m_vol);
                message += "\n";
                throw CanteraError("Reactor::updateState", message);
            }
        }
    } else {
        m_thermo->setDensity(m_mass/m_vol);
    }

    size_t loc = m_nsp + 3;
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
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

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
    double dmdt = 0.0; // dm/dt (gas phase)
    double* dYdt = ydot + 3;

    // compute wall terms
    size_t loc = m_nsp+3;
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

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    double mdot_surf = 0.0; // net mass flux from surfaces
    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        dYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k] / m_mass;
        mdot_surf += m_sdot[k] * mw[k];
    }
    dmdt += mdot_surf; // mass added to gas phase from surface reations

    for (size_t k = 0; k < m_nsp; k++) {
        // dilution by net surface mass flux
        dYdt[k] -= Y[k] * mdot_surf / m_mass;
    }

    /*
     *  Energy equation.
     *  \f[
     *  \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in}
     * - \dot m_{out} h.
     * \f]
     */
    if (m_energy) {
        ydot[2] = - m_thermo->pressure() * m_vdot - m_Q;
    } else {
        ydot[2] = 0.0;
    }

    // add terms for open system
    if (m_open) {
        doublereal enthalpy = m_thermo->enthalpy_mass();

        // outlets
        for (size_t i = 0; i < m_nOutlets; i++) {
            double mdot_out = m_outlet[i]->massFlowRate(time);
            dmdt -= mdot_out; // mass flow out of system
            if (m_energy) {
                ydot[2] -= mdot_out * enthalpy;
            }
        }

        // inlets
        for (size_t i = 0; i < m_nInlets; i++) {
            double mdot_in = m_inlet[i]->massFlowRate(time);
            dmdt += mdot_in; // mass flow into system
            for (size_t n = 0; n < m_nsp; n++) {
                double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
                // flow of species into system and dilution by other species
                dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;
            }
            if (m_energy) {
                ydot[2] += mdot_in * m_inlet[i]->enthalpy_mass();
            }
        }
    }

    ydot[0] = dmdt;

    for (size_t i = 0; i < m_nv; i++) {
        AssertFinite(ydot[i], "Reactor::evalEqs",
                     "ydot[" + int2str(i) + "] is not finite");
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

    network().registerSensitivityReaction(this, rxn,
                                          name()+": "+m_kin->reactionString(rxn));
    m_pnum.push_back(rxn);
    m_mult_save.push_back(1.0);
}

std::vector<std::pair<void*, int> > Reactor::getSensitivityOrder() const
{
    std::vector<std::pair<void*, int> > order;
    order.push_back(std::make_pair(const_cast<Reactor*>(this), 0));
    for (size_t n = 0; n < m_nwalls; n++) {
        if (m_nsens_wall[n]) {
            order.push_back(std::make_pair(m_wall[n], m_lr[n]));
        }
    }
    return order;
}

size_t Reactor::componentIndex(const string& nm) const
{
    if (nm == "m") {
        return 0;
    }
    if (nm == "V") {
        return 1;
    }
    if (nm == "U") {
        return 2;
    }

    // check for a gas species name
    size_t k = m_thermo->speciesIndex(nm);
    if (k != npos) {
        return k + 3;
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
                return k + 3 + m_nsp + walloffset;
            } else {
                walloffset += th->nSpecies();
            }
        }
    }
    return npos;
}

}
