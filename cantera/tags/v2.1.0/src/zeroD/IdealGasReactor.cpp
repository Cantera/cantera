/**
 *  @file IdealGasReactor.cpp A zero-dimensional reactor
 */

#include "cantera/zeroD/IdealGasReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/ReactorNet.h"

#include <cfloat>

using namespace std;

namespace Cantera
{

void IdealGasReactor::setThermoMgr(ThermoPhase& thermo)
{
    //! @TODO: Add a method to ThermoPhase that indicates whether a given
    //! subclass is compatible with this reactor model
    if (thermo.eosType() != cIdealGas) {
        throw CanteraError("IdealGasReactor::setThermoMgr",
                           "Incompatible phase type provided");
    }
    Reactor::setThermoMgr(thermo);
}

void IdealGasReactor::getInitialConditions(double t0, size_t leny, double* y)
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

    // Set the third component to the temperature
    y[2] = m_thermo->temperature();

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

void IdealGasReactor::initialize(doublereal t0)
{
    m_thermo->restoreState(m_state);
    m_sdot.resize(m_nsp, 0.0);
    m_wdot.resize(m_nsp, 0.0);
    m_uk.resize(m_nsp, 0.0);
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
                    throw CanteraError("IdealGasReactor::initialize",
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

void IdealGasReactor::updateState(doublereal* y)
{
    for (size_t i = 0; i < m_nv; i++) {
        AssertFinite(y[i], "IdealGasReactor::updateState",
                     "y[" + int2str(i) + "] is not finite");
    }

    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the temperature, [3...K+3] are the mass fractions of each species,
    // and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];

    m_thermo->setMassFractions_NoNorm(y+3);
    m_thermo->setState_TR(y[2], m_mass / m_vol);

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

void IdealGasReactor::evalEqs(doublereal time, doublereal* y,
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
    double mcvdTdt = 0.0; // m * c_v * dT/dt
    double dmdt = 0.0; // dm/dt (gas phase)
    double* dYdt = ydot + 3;

    m_thermo->getPartialMolarIntEnergies(&m_uk[0]);

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
    dmdt += mdot_surf;

    // compression work and external heat transfer
    mcvdTdt += - m_pressure * m_vdot - m_Q;

    for (size_t n = 0; n < m_nsp; n++) {
        // heat release from gas phase and surface reations
        mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol;
        mcvdTdt -= m_sdot[n] * m_uk[n];
        // dilution by net surface mass flux
        dYdt[n] -= Y[n] * mdot_surf / m_mass;
    }

    // add terms for open system
    if (m_open) {
        // outlets
        for (size_t i = 0; i < m_nOutlets; i++) {
            double mdot_out = m_outlet[i]->massFlowRate(time);
            dmdt -= mdot_out; // mass flow out of system
            mcvdTdt -= mdot_out * m_pressure * m_vol / m_mass; // flow work
        }

        // inlets
        for (size_t i = 0; i < m_nInlets; i++) {
            double mdot_in = m_inlet[i]->massFlowRate(time);
            dmdt += mdot_in; // mass flow into system
            mcvdTdt += m_inlet[i]->enthalpy_mass() * mdot_in;
            for (size_t n = 0; n < m_nsp; n++) {
                double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
                // flow of species into system and dilution by other species
                dYdt[n] += (mdot_spec - mdot_in * Y[n]) / m_mass;

                // In combintion with h_in*mdot_in, flow work plus thermal
                // energy carried with the species
                mcvdTdt -= m_uk[n] / mw[n] * mdot_spec;
            }
        }
    }

    ydot[0] = dmdt;
    ydot[1] = m_vdot;
    if (m_energy) {
        ydot[2] = mcvdTdt / (m_mass * m_thermo->cv_mass());
    } else {
        ydot[2] = 0;
    }

    for (size_t i = 0; i < m_nv; i++) {
        AssertFinite(ydot[i], "IdealGasReactor::evalEqs",
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

size_t IdealGasReactor::componentIndex(const string& nm) const
{
    if (nm == "T") {
        return 2;
    } else {
        return Reactor::componentIndex(nm);
    }
}

}
