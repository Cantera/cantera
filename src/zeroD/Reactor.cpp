//! @file Reactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/Reactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/ReactorNet.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/Solution.h"
#include "cantera/base/utilities.h"

#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{
Reactor::Reactor() :
    m_kin(0),
    m_vdot(0.0),
    m_Qdot(0.0),
    m_mass(0.0),
    m_chem(false),
    m_energy(true),
    m_nv(0)
{}

void Reactor::insert(shared_ptr<Solution> sol) {
    setThermoMgr(*sol->thermo());
    setKineticsMgr(*sol->kinetics());
}

void Reactor::setDerivativeSettings(AnyMap& settings)
{
    m_kin->setDerivativeSettings(settings);
}

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
    updateConnected(true);

    for (size_t n = 0; n < m_wall.size(); n++) {
        WallBase* W = m_wall[n];
        W->initialize();
    }

    m_nv = m_nsp + 3;
    m_nv_surf = 0;
    size_t maxnt = 0;
    for (auto& S : m_surfaces) {
        m_nv_surf += S->thermo()->nSpecies();
        size_t nt = S->kinetics()->nTotalSpecies();
        maxnt = std::max(maxnt, nt);
    }
    m_nv += m_nv_surf;
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
            // Try full-range bisection if bracketing fails (for example, near
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

    updateConnected(true);
    updateSurfaceState(y + m_nsp + 3);
}

void Reactor::updateSurfaceState(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        S->setCoverages(y+loc);
        loc += S->thermo()->nSpecies();
    }
}

void Reactor::updateConnected(bool updatePressure) {
    // save parameters needed by other connected reactors
    m_enthalpy = m_thermo->enthalpy_mass();
    if (updatePressure) {
        m_pressure = m_thermo->pressure();
    }
    m_intEnergy = m_thermo->intEnergy_mass();
    m_thermo->saveState(m_state);

    // Update the mass flow rate of connected flow devices
    double time = (m_net != nullptr) ? m_net->time() : 0.0;
    for (size_t i = 0; i < m_outlet.size(); i++) {
        m_outlet[i]->updateMassFlowRate(time);
    }
    for (size_t i = 0; i < m_inlet.size(); i++) {
        m_inlet[i]->updateMassFlowRate(time);
    }
}

void Reactor::eval(double time, double* LHS, double* RHS)
{
    double& dmdt = RHS[0];
    double* mdYdt = RHS + 3; // mass * dY/dt

    evalWalls(time);
    m_thermo->restoreState(m_state);
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* Y = m_thermo->massFractions();

    evalSurfaces(LHS + m_nsp + 3, RHS + m_nsp + 3, m_sdot.data());
     // mass added to gas phase from surface reactions
    double mdot_surf = dot(m_sdot.begin(), m_sdot.end(), mw.begin());
    dmdt = mdot_surf;

    // volume equation
    RHS[1] = m_vdot;

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        mdYdt[k] = (m_wdot[k] * m_vol + m_sdot[k]) * mw[k];
        // dilution by net surface mass flux
        mdYdt[k] -= Y[k] * mdot_surf;
        LHS[k+3] = m_mass;
    }

    // Energy equation.
    // \f[
    //     \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in} - \dot m_{out} h.
    // \f]
    if (m_energy) {
        RHS[2] = - m_thermo->pressure() * m_vdot + m_Qdot;
    } else {
        RHS[2] = 0.0;
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        double mdot = outlet->massFlowRate();
        dmdt -= mdot; // mass flow out of system
        if (m_energy) {
            RHS[2] -= mdot * m_enthalpy;
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        dmdt += mdot; // mass flow into system
        for (size_t n = 0; n < m_nsp; n++) {
            double mdot_spec = inlet->outletSpeciesMassFlowRate(n);
            // flow of species into system and dilution by other species
            mdYdt[n] += (mdot_spec - mdot * Y[n]);
        }
        if (m_energy) {
            RHS[2] += mdot * inlet->enthalpy_mass();
        }
    }
}

void Reactor::evalWalls(double t)
{
    m_vdot = 0.0;
    m_Qdot = 0.0;
    for (size_t i = 0; i < m_wall.size(); i++) {
        int f = 2 * m_lr[i] - 1;
        m_vdot -= f * m_wall[i]->vdot(t);
        m_Qdot += f * m_wall[i]->Q(t);
    }
}

void Reactor::evalSurfaces(double* LHS, double* RHS, double* sdot)
{
    fill(sdot, sdot + m_nsp, 0.0);
    size_t loc = 0; // offset into ydot

    for (auto S : m_surfaces) {
        Kinetics* kin = S->kinetics();
        SurfPhase* surf = S->thermo();

        double rs0 = 1.0/surf->siteDensity();
        size_t nk = surf->nSpecies();
        double sum = 0.0;
        S->syncState();
        kin->getNetProductionRates(&m_work[0]);
        size_t ns = kin->surfacePhaseIndex();
        size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
        for (size_t k = 1; k < nk; k++) {
            LHS[loc] = 1.0;
            RHS[loc + k] = m_work[surfloc + k] * rs0 * surf->size(k);
            sum -= RHS[loc + k];
        }
        LHS[loc] = 1.0;
        RHS[loc] = sum;
        loc += nk;

        size_t bulkloc = kin->kineticsSpeciesIndex(m_thermo->speciesName(0));
        double wallarea = S->area();
        for (size_t k = 0; k < m_nsp; k++) {
            sdot[k] += m_work[bulkloc + k] * wallarea;
        }
    }
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
    if (k == npos) {
        throw CanteraError("Reactor::setAdvanceLimit", "No component named '{}'", nm);
    }

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
