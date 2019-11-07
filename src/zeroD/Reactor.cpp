//! @file Reactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"

#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{
Reactor::Reactor() :
    m_kin(0),
    m_vdot(0.0),
    m_Q(0.0),
    m_mass(0.0),
    m_chem(false),
    m_energy(true),
    m_nv(0)
{}

void Reactor::setKineticsMgr(Kinetics& kin)
{
    m_kin = &kin;
    if (m_kin->nReactions() == 0) {
        setChemistry(false);
    } else {
        setChemistry(true);
    }
}

void Reactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("Reactor::getState",
                           "Error: reactor is empty.");
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
    getSurfaceInitialConditions(y + m_nsp + 3);
}

void Reactor::getSurfaceInitialConditions(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        S->getCoverages(y + loc);
        loc += S->thermo()->nSpecies();
    }
}

void Reactor::initialize(doublereal t0)
{
    if (!m_thermo || (m_chem && !m_kin)) {
        throw CanteraError("Reactor::initialize", "Reactor contents not set"
                " for reactor '" + m_name + "'.");
    }
    m_thermo->restoreState(m_state);
    m_sdot.resize(m_nsp, 0.0);
    m_wdot.resize(m_nsp, 0.0);

    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();

    for (size_t n = 0; n < m_wall.size(); n++) {
        WallBase* W = m_wall[n];
        W->initialize();
    }

    m_nv = m_nsp + 3;
    size_t maxnt = 0;
    for (auto& S : m_surfaces) {
        m_nv += S->thermo()->nSpecies();
        size_t nt = S->kinetics()->nTotalSpecies();
        maxnt = std::max(maxnt, nt);
        if (m_chem && &m_kin->thermo(0) != &S->kinetics()->thermo(0)) {
            throw CanteraError("Reactor::initialize",
                "First phase of all kinetics managers must be the gas.");
        }
    }
    m_work.resize(maxnt);
}

size_t Reactor::nSensParams()
{
    size_t ns = m_sensParams.size();
    for (auto& S : m_surfaces) {
        ns += S->nSensParams();
    }
    return ns;
}

void Reactor::syncState()
{
    ReactorBase::syncState();
    m_mass = m_thermo->density() * m_vol;
}

void Reactor::updateState(doublereal* y)
{
    // The components of y are [0] the total mass, [1] the total volume,
    // [2] the total internal energy, [3...K+3] are the mass fractions of each
    // species, and [K+3...] are the coverages of surface species on each wall.
    m_mass = y[0];
    m_vol = y[1];
    m_thermo->setMassFractions_NoNorm(y+3);

    if (m_energy) {
        double U = y[2];
        // Residual function: error in internal energy as a function of T
        auto u_err = [this, U](double T) {
            m_thermo->setState_TR(T, m_mass / m_vol);
            return m_thermo->intEnergy_mass() * m_mass - U;
        };

        double T = m_thermo->temperature();
        boost::uintmax_t maxiter = 100;
        std::pair<double, double> TT;
        try {
            TT = bmt::bracket_and_solve_root(
                u_err, T, 1.2, true, bmt::eps_tolerance<double>(48), maxiter);
        } catch (std::exception&) {
            // Try full-range bisection if bracketing fails (e.g. near
            // temperature limits for the phase's equation of state)
            try {
                TT = bmt::bisect(u_err, m_thermo->minTemp(), m_thermo->maxTemp(),
                    bmt::eps_tolerance<double>(48), maxiter);
            } catch (std::exception& err2) {
                // Set m_thermo back to a reasonable state if root finding fails
                m_thermo->setState_TR(T, m_mass / m_vol);
                throw CanteraError("Reactor::updateState",
                    "{}\nat U = {}, rho = {}", err2.what(), U, m_mass / m_vol);
            }
        }
        if (fabs(TT.first - TT.second) > 1e-7*TT.first) {
            throw CanteraError("Reactor::updateState", "root finding failed");
        }
        m_thermo->setState_TR(TT.second, m_mass / m_vol);
    } else {
        m_thermo->setDensity(m_mass/m_vol);
    }

    updateSurfaceState(y + m_nsp + 3);

    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    m_pressure = m_thermo->pressure();
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);
}

void Reactor::updateSurfaceState(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        S->setCoverages(y+loc);
        loc += S->thermo()->nSpecies();
    }
}

void Reactor::evalEqs(doublereal time, doublereal* y,
                      doublereal* ydot, doublereal* params)
{
    double dmdt = 0.0; // dm/dt (gas phase)
    double* dYdt = ydot + 3;

    evalFlowDevices(time);
    evalWalls(time);
    applySensitivity(params);
    m_thermo->restoreState(m_state);
    double mdot_surf = evalSurfaces(time, ydot + m_nsp + 3);
    dmdt += mdot_surf; // mass added to gas phase from surface reactions

    // volume equation
    ydot[1] = m_vdot;

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        dYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k] / m_mass;
        // dilution by net surface mass flux
        dYdt[k] -= Y[k] * mdot_surf / m_mass;
    }

    // Energy equation.
    // \f[
    //     \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in} - \dot m_{out} h.
    // \f]
    if (m_energy) {
        ydot[2] = - m_thermo->pressure() * m_vdot - m_Q;
    } else {
        ydot[2] = 0.0;
    }

    // add terms for outlets
    for (size_t i = 0; i < m_outlet.size(); i++) {
        dmdt -= m_mdot_out[i]; // mass flow out of system
        if (m_energy) {
            ydot[2] -= m_mdot_out[i] * m_enthalpy;
        }
    }

    // add terms for inlets
    for (size_t i = 0; i < m_inlet.size(); i++) {
        dmdt += m_mdot_in[i]; // mass flow into system
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = m_inlet[i]->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            dYdt[n] += (mdot_spec - m_mdot_in[i] * Y[n]) / m_mass;
        }
        if (m_energy) {
            ydot[2] += m_mdot_in[i] * m_inlet[i]->enthalpy_mass();
        }
    }

    ydot[0] = dmdt;
    resetSensitivity(params);
}

void Reactor::evalWalls(double t)
{
    m_vdot = 0.0;
    m_Q = 0.0;
    for (size_t i = 0; i < m_wall.size(); i++) {
        int lr = 1 - 2*m_lr[i];
        m_vdot += lr*m_wall[i]->vdot(t);
        m_Q += lr*m_wall[i]->Q(t);
    }
}

void Reactor::evalFlowDevices(double t)
{
    for (size_t i = 0; i < m_outlet.size(); i++) {
        m_mdot_out[i] = m_outlet[i]->massFlowRate(t);
    }
    for (size_t i = 0; i < m_inlet.size(); i++) {
        m_mdot_in[i] = m_inlet[i]->massFlowRate(t);
    }
}

double Reactor::evalSurfaces(double t, double* ydot)
{
    const vector_fp& mw = m_thermo->molecularWeights();
    fill(m_sdot.begin(), m_sdot.end(), 0.0);
    size_t loc = 0; // offset into ydot
    double mdot_surf = 0.0; // net mass flux from surface

    for (auto S : m_surfaces) {
        Kinetics* kin = S->kinetics();
        SurfPhase* surf = S->thermo();

        double rs0 = 1.0/surf->siteDensity();
        size_t nk = surf->nSpecies();
        double sum = 0.0;
        surf->setTemperature(m_state[0]);
        S->syncCoverages();
        kin->getNetProductionRates(&m_work[0]);
        size_t ns = kin->surfacePhaseIndex();
        size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
        for (size_t k = 1; k < nk; k++) {
            ydot[loc + k] = m_work[surfloc+k]*rs0*surf->size(k);
            sum -= ydot[loc + k];
        }
        ydot[loc] = sum;
        loc += nk;

        double wallarea = S->area();
        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += m_work[k]*wallarea;
            mdot_surf += m_sdot[k] * mw[k];
        }
    }
    return mdot_surf;
}

void Reactor::addSensitivityReaction(size_t rxn)
{
    if (!m_chem || rxn >= m_kin->nReactions()) {
        throw CanteraError("Reactor::addSensitivityReaction",
                           "Reaction number out of range ({})", rxn);
    }

    size_t p = network().registerSensitivityParameter(
        name()+": "+m_kin->reactionString(rxn), 1.0, 1.0);
    m_sensParams.emplace_back(
        SensitivityParameter{rxn, p, 1.0, SensParameterType::reaction});
}

void Reactor::addSensitivitySpeciesEnthalpy(size_t k)
{
    if (k >= m_thermo->nSpecies()) {
        throw CanteraError("Reactor::addSensitivitySpeciesEnthalpy",
                           "Species index out of range ({})", k);
    }

    size_t p = network().registerSensitivityParameter(
        name() + ": " + m_thermo->speciesName(k) + " enthalpy",
        0.0, GasConstant * 298.15);
    m_sensParams.emplace_back(
        SensitivityParameter{k, p, m_thermo->Hf298SS(k),
                             SensParameterType::enthalpy});
}

size_t Reactor::speciesIndex(const string& nm) const
{
    // check for a gas species name
    size_t k = m_thermo->speciesIndex(nm);
    if (k != npos) {
        return k;
    }

    // check for a wall species
    size_t offset = m_nsp;
    for (auto& S : m_surfaces) {
        ThermoPhase* th = S->thermo();
        k = th->speciesIndex(nm);
        if (k != npos) {
            return k + offset;
        } else {
            offset += th->nSpecies();
        }
    }
    return npos;
}

size_t Reactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + 3;
    } else if (nm == "mass") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else if (nm == "int_energy") {
        return 2;
    } else {
        return npos;
    }
}

std::string Reactor::componentName(size_t k) {
    if (k == 0) {
        return "mass";
    } else if (k == 1) {
        return "volume";
    } else if (k == 2) {
        return "int_energy";
    } else if (k >= 3 && k < neq()) {
        k -= 3;
        if (k < m_thermo->nSpecies()) {
            return m_thermo->speciesName(k);
        } else {
            k -= m_thermo->nSpecies();
        }
        for (auto& S : m_surfaces) {
            ThermoPhase* th = S->thermo();
            if (k < th->nSpecies()) {
                return th->speciesName(k);
            } else {
                k -= th->nSpecies();
            }
        }
    }
    throw CanteraError("Reactor::componentName", "Index is out of bounds.");
}

void Reactor::applySensitivity(double* params)
{
    if (!params) {
        return;
    }
    for (auto& p : m_sensParams) {
        if (p.type == SensParameterType::reaction) {
            p.value = m_kin->multiplier(p.local);
            m_kin->setMultiplier(p.local, p.value*params[p.global]);
        } else if (p.type == SensParameterType::enthalpy) {
            m_thermo->modifyOneHf298SS(p.local, p.value + params[p.global]);
        }
    }
    for (auto& S : m_surfaces) {
        S->setSensitivityParameters(params);
    }
    m_thermo->invalidateCache();
    if (m_kin) {
        m_kin->invalidateCache();
    }
}

void Reactor::resetSensitivity(double* params)
{
    if (!params) {
        return;
    }
    for (auto& p : m_sensParams) {
        if (p.type == SensParameterType::reaction) {
            m_kin->setMultiplier(p.local, p.value);
        } else if (p.type == SensParameterType::enthalpy) {
            m_thermo->resetHf298(p.local);
        }
    }
    for (auto& S : m_surfaces) {
        S->resetSensitivityParameters();
    }
    m_thermo->invalidateCache();
    if (m_kin) {
        m_kin->invalidateCache();
    }
}

void Reactor::setAdvanceLimits(const double *limits)
{
    if (m_thermo == 0) {
        throw CanteraError("Reactor::setAdvanceLimits",
                           "Error: reactor is empty.");
    }
    m_advancelimits.assign(limits, limits + m_nv);

    // resize to zero length if no limits are set
    if (std::none_of(m_advancelimits.begin(), m_advancelimits.end(),
                     [](double val){return val>0;})) {
        m_advancelimits.resize(0);
    }
}

bool Reactor::getAdvanceLimits(double *limits)
{
    bool has_limit = hasAdvanceLimits();
    if (has_limit) {
        std::copy(m_advancelimits.begin(), m_advancelimits.end(), limits);
    } else {
        std::fill(limits, limits + m_nv, -1.0);
    }
    return has_limit;
}

void Reactor::setAdvanceLimit(const string& nm, const double limit)
{
    size_t k = componentIndex(nm);

    if (m_thermo == 0) {
        throw CanteraError("Reactor::setAdvanceLimit",
                           "Error: reactor is empty.");
    }
    if (m_nv == 0) {
        if (m_net == 0) {
            throw CanteraError("Reactor::setAdvanceLimit",
                               "Cannot set limit on a reactor that is not "
                               "assigned to a ReactorNet object.");
        } else {
            m_net->initialize();
        }
    } else if (k > m_nv) {
        throw CanteraError("Reactor::setAdvanceLimit",
                           "Index out of bounds.");
    }
    m_advancelimits.resize(m_nv, -1.0);
    m_advancelimits[k] = limit;

    // resize to zero length if no limits are set
    if (std::none_of(m_advancelimits.begin(), m_advancelimits.end(),
                     [](double val){return val>0;})) {
        m_advancelimits.resize(0);
    }
}

}
