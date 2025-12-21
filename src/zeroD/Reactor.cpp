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
#include "cantera/kinetics/Reaction.h"
#include "cantera/base/Solution.h"
#include "cantera/base/utilities.h"

#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

Reactor::Reactor(shared_ptr<Solution> sol, const string& name)
    : Reactor(sol, true, name)
{
}

Reactor::Reactor(shared_ptr<Solution> sol, bool clone, const string& name)
    : ReactorBase(sol, clone, name)
{
    m_kin = m_solution->kinetics().get();
    setChemistryEnabled(m_kin->nReactions() > 0);
    m_vol = 1.0; // By default, the volume is set to 1.0 m^3.
    m_sdot.resize(m_nsp, 0.0);
}

void Reactor::setDerivativeSettings(AnyMap& settings)
{
    m_kin->setDerivativeSettings(settings);
}

void Reactor::getState(double* y)
{
    // set the first component to the total mass
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_mass;

    // set the second component to the total volume
    y[1] = m_vol;

    // set the third component to the total internal energy
    y[2] = m_thermo->intEnergy_mass() * m_mass;

    // set components y+3 ... y+K+2 to the mass fractions of each species
    m_thermo->getMassFractions(y+3);
}

void Reactor::initialize(double t0)
{
    if (!m_thermo || (m_chem && !m_kin)) {
        throw CanteraError("Reactor::initialize", "Reactor contents not set"
                " for reactor '" + m_name + "'.");
    }
    m_wdot.resize(m_nsp, 0.0);
    updateConnected(true);

    for (size_t n = 0; n < m_wall.size(); n++) {
        WallBase* W = m_wall[n];
        W->initialize();
    }

    m_nv = m_nsp + 3;
}

void Reactor::updateState(double* y)
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
            m_thermo->setState_TD(T, m_mass / m_vol);
            return m_thermo->intEnergy_mass() * m_mass - U;
        };

        double T = m_thermo->temperature();
        boost::uintmax_t maxiter = 100;
        pair<double, double> TT;
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
                m_thermo->setState_TD(T, m_mass / m_vol);
                throw CanteraError("Reactor::updateState",
                    "{}\nat U = {}, rho = {}", err2.what(), U, m_mass / m_vol);
            }
        }
        if (fabs(TT.first - TT.second) > 1e-7*TT.first) {
            throw CanteraError("Reactor::updateState", "root finding failed");
        }
        m_thermo->setState_TD(TT.second, m_mass / m_vol);
    } else {
        m_thermo->setDensity(m_mass/m_vol);
    }

    updateConnected(true);
}

void Reactor::eval(double time, double* LHS, double* RHS)
{
    double& dmdt = RHS[0];
    double* mdYdt = RHS + 3; // mass * dY/dt

    evalWalls(time);
    updateSurfaceProductionRates();
    const vector<double>& mw = m_thermo->molecularWeights();
    const double* Y = m_thermo->massFractions();

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
    // @f[
    //     \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in} - \dot m_{out} h.
    // @f]
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
    // time is currently unused
    m_vdot = 0.0;
    m_Qdot = 0.0;
    for (size_t i = 0; i < m_wall.size(); i++) {
        int f = 2 * m_lr[i] - 1;
        m_vdot -= f * m_wall[i]->expansionRate();
        m_Qdot += f * m_wall[i]->heatRate();
    }
}

vector<size_t> Reactor::steadyConstraints() const
{
    if (!energyEnabled()) {
        throw CanteraError("Reactor::steadyConstraints", "Steady state solver cannot"
            " be used with {0} when energy equation is disabled."
            "\nConsider using IdealGas{0} instead.\n"
            "See https://github.com/Cantera/enhancements/issues/234", type());
    }
    if (nSurfs() != 0) {
        throw CanteraError("Reactor::steadyConstraints", "Steady state solver cannot"
            " currently be used when reactor surfaces are present.\n"
            "See https://github.com/Cantera/enhancements/issues/234.");
    }
    return {1}; // volume
}

Eigen::SparseMatrix<double> Reactor::finiteDifferenceJacobian()
{
    if (m_nv == 0) {
        throw CanteraError("Reactor::finiteDifferenceJacobian",
                           "Reactor must be initialized first.");
    }
    vector<Eigen::Triplet<double>> trips;
    Eigen::ArrayXd yCurrent(m_nv);
    getState(yCurrent.data());
    double time = (m_net != nullptr) ? m_net->time() : 0.0;

    Eigen::ArrayXd yPerturbed = yCurrent;
    Eigen::ArrayXd lhsPerturbed(m_nv), lhsCurrent(m_nv);
    Eigen::ArrayXd rhsPerturbed(m_nv), rhsCurrent(m_nv);
    lhsCurrent = 1.0;
    rhsCurrent = 0.0;
    updateState(yCurrent.data());
    eval(time, lhsCurrent.data(), rhsCurrent.data());

    double rel_perturb = std::sqrt(std::numeric_limits<double>::epsilon());
    double atol = (m_net != nullptr) ? m_net->atol() : 1e-15;

    for (size_t j = 0; j < m_nv; j++) {
        yPerturbed = yCurrent;
        double delta_y = std::max(std::abs(yCurrent[j]), 1000 * atol) * rel_perturb;
        yPerturbed[j] += delta_y;

        updateState(yPerturbed.data());
        lhsPerturbed = 1.0;
        rhsPerturbed = 0.0;
        eval(time, lhsPerturbed.data(), rhsPerturbed.data());

        // d ydot_i/dy_j
        for (size_t i = 0; i < m_nv; i++) {
            double ydotPerturbed = rhsPerturbed[i] / lhsPerturbed[i];
            double ydotCurrent = rhsCurrent[i] / lhsCurrent[i];
            if (ydotCurrent != ydotPerturbed) {
                trips.emplace_back(static_cast<int>(i), static_cast<int>(j),
                                   (ydotPerturbed - ydotCurrent) / delta_y);
            }
        }
    }
    updateState(yCurrent.data());

    Eigen::SparseMatrix<double> jac(m_nv, m_nv);
    jac.setFromTriplets(trips.begin(), trips.end());
    return jac;
}

void Reactor::addSensitivityReaction(size_t rxn)
{
    if (!m_chem || rxn >= m_kin->nReactions()) {
        throw CanteraError("Reactor::addSensitivityReaction",
                           "Reaction number out of range ({})", rxn);
    }

    size_t p = network().registerSensitivityParameter(
        name()+": "+m_kin->reaction(rxn)->equation(), 1.0, 1.0);
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

void Reactor::updateSurfaceProductionRates()
{
    m_sdot.assign(m_nsp, 0.0);
    for (auto& S : m_surfaces) {
        const auto& sdot = S->surfaceProductionRates();
        size_t offset = S->kinetics()->kineticsSpeciesIndex(m_thermo->speciesName(0));
        for (size_t k = 0; k < m_nsp; k++) {
            m_sdot[k] += sdot[offset + k] * S->area();
        }
    }
}

size_t Reactor::componentIndex(const string& nm) const
{
    if (nm == "mass") {
        return 0;
    }
    if (nm == "volume") {
        return 1;
    }
    if (nm == "int_energy") {
        return 2;
    }
    try {
        return m_thermo->speciesIndex(nm) + 3;
    } catch (const CanteraError&) {
        throw CanteraError("Reactor::componentIndex",
            "Component '{}' not found", nm);
    }
}

string Reactor::componentName(size_t k) {
    if (k == 0) {
        return "mass";
    } else if (k == 1) {
        return "volume";
    } else if (k == 2) {
        return "int_energy";
    } else if (k >= 3 && k < neq()) {
        return m_thermo->speciesName(k - 3);
    }
    throw IndexError("Reactor::componentName", "component", k, m_nv);
}

double Reactor::upperBound(size_t k) const {
    if (k == 0) {
        return BigNumber; // mass
    } else if (k == 1) {
        return BigNumber; // volume
    } else if (k == 2) {
        return BigNumber; // internal energy
    } else if (k >= 3 && k < m_nv) {
        return 1.0; // species mass fraction or surface coverage
    } else {
        throw CanteraError("Reactor::upperBound", "Index {} is out of bounds.", k);
    }
}

double Reactor::lowerBound(size_t k) const {
    if (k == 0) {
        return 0; // mass
    } else if (k == 1) {
        return 0; // volume
    } else if (k == 2) {
        return -BigNumber; // internal energy
    } else if (k >= 3 && k < m_nv) {
        return -Tiny; // species mass fraction or surface coverage
    } else {
        throw CanteraError("Reactor::lowerBound", "Index {} is out of bounds.", k);
    }
}

void Reactor::resetBadValues(double* y) {
    for (size_t k = 3; k < m_nv; k++) {
        y[k] = std::max(y[k], 0.0);
    }
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
    m_thermo->invalidateCache();
    if (m_kin) {
        m_kin->invalidateCache();
    }
}

void Reactor::setAdvanceLimits(const double *limits)
{
    m_advancelimits.assign(limits, limits + m_nv);

    // resize to zero length if no limits are set
    if (std::none_of(m_advancelimits.begin(), m_advancelimits.end(),
                     [](double val){return val>0;})) {
        m_advancelimits.resize(0);
    }
}

bool Reactor::getAdvanceLimits(double *limits) const
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
