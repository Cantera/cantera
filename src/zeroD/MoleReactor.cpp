//! @file MoleReactor.cpp A zero-dimensional reactor with a moles as the
//! state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/MoleReactor.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/base/utilities.h"
#include <boost/math/tools/roots.hpp>

using namespace std;
namespace bmt = boost::math::tools;

namespace Cantera
{

void MoleReactor::getSurfaceInitialConditions(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        double area = S->area();
        auto currPhase = S->thermo();
        size_t tempLoc = currPhase->nSpecies();
        double surfDensity = currPhase->siteDensity();
        S->getCoverages(y + loc);
        // convert coverages to moles
        for (size_t i = 0; i < tempLoc; i++) {
            y[i + loc] = y[i + loc] * area * surfDensity / currPhase->size(i);
        }
        loc += tempLoc;
    }
}

void MoleReactor::initialize(double t0)
{
    Reactor::initialize(t0);
    m_nv -= 1; // moles gives the state one fewer variables
}

void MoleReactor::updateSurfaceState(double* y)
{
    size_t loc = 0;
    vector<double> coverages(m_nv_surf, 0.0);
    for (auto& S : m_surfaces) {
        auto surf = S->thermo();
        double invArea = 1/S->area();
        double invSurfDensity = 1/surf->siteDensity();
        size_t tempLoc = surf->nSpecies();
        for (size_t i = 0; i < tempLoc; i++) {
            coverages[i + loc] = y[i + loc] * invArea * surf->size(i) * invSurfDensity;
        }
        S->setCoverages(coverages.data()+loc);
        loc += tempLoc;
    }
}

void MoleReactor::evalSurfaces(double* LHS, double* RHS, double* sdot)
{
    fill(sdot, sdot + m_nsp, 0.0);
    size_t loc = 0; // offset into ydot
    for (auto S : m_surfaces) {
        Kinetics* kin = S->kinetics();
        SurfPhase* surf = S->thermo();
        double wallarea = S->area();
        size_t nk = surf->nSpecies();
        S->restoreState();
        kin->getNetProductionRates(&m_work[0]);
        for (size_t k = 0; k < nk; k++) {
            RHS[loc + k] = m_work[k] * wallarea / surf->size(k);
        }
        loc += nk;

        size_t bulkloc = kin->kineticsSpeciesIndex(m_thermo->speciesName(0));

        for (size_t k = 0; k < m_nsp; k++) {
            sdot[k] += m_work[bulkloc + k] * wallarea;
        }
    }
}

void MoleReactor::addSurfaceJacobian(vector<Eigen::Triplet<double>> &triplets)
{
    size_t offset = m_nsp;
    for (auto& S : m_surfaces) {
        S->restoreState();
        double A = S->area();
        auto kin = S->kinetics();
        size_t nk = S->thermo()->nSpecies();
        // index of gas and surface phases to check if the species is in gas or surface
        size_t spi = 0;
        size_t gpi = kin->speciesPhaseIndex(kin->kineticsSpeciesIndex(
            m_thermo->speciesName(0)));
        // get surface jacobian in concentration units
        Eigen::SparseMatrix<double> surfJac = kin->netProductionRates_ddCi();
        // loop through surface specific jacobian and add elements to triplets vector
        // accordingly
        for (int k=0; k<surfJac.outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(surfJac, k); it; ++it) {
                size_t row = it.row();
                size_t col = it.col();
                auto& rowPhase = kin->speciesPhase(row);
                auto& colPhase = kin->speciesPhase(col);
                size_t rpi = kin->phaseIndex(rowPhase.name());
                size_t cpi = kin->phaseIndex(colPhase.name());
                // check if the reactor kinetics object contains both phases to avoid
                // any solid phases which may be included then use phases to map surf
                // kinetics indicies to reactor kinetic indices
                if ((rpi == spi || rpi == gpi) && (cpi == spi || cpi == gpi) ) {
                    // subtract start of phase
                    row -= kin->kineticsSpeciesIndex(0, rpi);
                    col -= kin->kineticsSpeciesIndex(0, cpi);
                    // since the gas phase is the first phase in the reactor state
                    // vector add the offset only if it is a surf species index to
                    // both row and col
                    row = (rpi != spi) ? row : row + offset;
                    // determine appropriate scalar to account for dimensionality
                    // gas phase species indices will be less than m_nsp
                    // so use volume if that is the case or area otherwise
                    double scalar = A;
                    if (cpi == spi) {
                        col += offset;
                        scalar /= A;
                    } else {
                        scalar /= m_vol;
                    }
                    // push back scaled value triplet
                    triplets.emplace_back(static_cast<int>(row), static_cast<int>(col),
                                          scalar * it.value());
                }
            }
        }
        // add species in this surface to the offset
        offset += nk;
    }
}

void MoleReactor::getMoles(double* y)
{
    // Use inverse molecular weights to convert to moles
    const double* Y = m_thermo->massFractions();
    const vector<double>& imw = m_thermo->inverseMolecularWeights();
    for (size_t i = 0; i < m_nsp; i++) {
        y[i] = m_mass * imw[i] * Y[i];
    }
}

void MoleReactor::setMassFromMoles(double* y)
{
    const vector<double>& mw = m_thermo->molecularWeights();
    // calculate mass from moles
    m_mass = 0;
    for (size_t i = 0; i < m_nsp; i++) {
        m_mass += y[i] * mw[i];
    }
}

void MoleReactor::getState(double* y)
{
    if (m_thermo == 0) {
        throw CanteraError("MoleReactor::getState",
                           "Error: reactor is empty.");
    }
    m_thermo->restoreState(m_state);
    // set the first component to the internal energy
    m_mass = m_thermo->density() * m_vol;
    y[0] = m_thermo->intEnergy_mass() * m_mass;
    // set the second component to the total volume
    y[1] = m_vol;
    // set components y+2 ... y+K+2 to the moles of each species
    getMoles(y + m_sidx);
    // set the remaining components to the surface species
    // moles on walls
    getSurfaceInitialConditions(y+m_nsp+m_sidx);
}

void MoleReactor::updateState(double* y)
{
    // The components of y are [0] total internal energy, [1] the total volume, and
    // [2...K+3] are the moles of each species, and [K+3...] are the moles
    // of surface species on each wall.
    setMassFromMoles(y + m_sidx);
    m_vol = y[1];
    m_thermo->setMolesNoTruncate(y + m_sidx);
    if (m_energy) {
        double U = y[0];
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
                throw CanteraError("MoleReactor::updateState",
                    "{}\nat U = {}, rho = {}", err2.what(), U, m_mass / m_vol);
            }
        }
        if (fabs(TT.first - TT.second) > 1e-7*TT.first) {
            throw CanteraError("MoleReactor::updateState", "root finding failed");
        }
        m_thermo->setState_TD(TT.second, m_mass / m_vol);
    } else {
        m_thermo->setDensity(m_mass / m_vol);
    }
    updateConnected(true);
    updateSurfaceState(y + m_nsp + m_sidx);
}

void MoleReactor::eval(double time, double* LHS, double* RHS)
{
    double* dndt = RHS + m_sidx; // moles per time

    evalWalls(time);
    m_thermo->restoreState(m_state);

    evalSurfaces(LHS + m_nsp + m_sidx, RHS + m_nsp + m_sidx, m_sdot.data());
    // inverse molecular weights for conversion
    const vector<double>& imw = m_thermo->inverseMolecularWeights();
    // volume equation
    RHS[1] = m_vdot;

    if (m_chem) {
        m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
    }

    // Energy equation.
    // @f[
    //     \dot U = -P\dot V + A \dot q + \dot m_{in} h_{in} - \dot m_{out} h.
    // @f]
    if (m_energy) {
        RHS[0] = - m_thermo->pressure() * m_vdot + m_Qdot;
    } else {
        RHS[0] = 0.0;
    }

    for (size_t k = 0; k < m_nsp; k++) {
        // production in gas phase and from surfaces
        dndt[k] = m_wdot[k] * m_vol + m_sdot[k];
    }

    // add terms for outlets
    for (auto outlet : m_outlet) {
        // flow of species into system and dilution by other species
        for (size_t n = 0; n < m_nsp; n++) {
            dndt[n] -= outlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
        // energy update based on mass flow
        double mdot = outlet->massFlowRate();
        if (m_energy) {
            RHS[0] -= mdot * m_enthalpy;
        }
    }

    // add terms for inlets
    for (auto inlet : m_inlet) {
        double mdot = inlet->massFlowRate();
        for (size_t n = 0; n < m_nsp; n++) {
            // flow of species into system and dilution by other species
            dndt[n] += inlet->outletSpeciesMassFlowRate(n) * imw[n];
        }
        if (m_energy) {
            RHS[0] += mdot * inlet->enthalpy_mass();
        }
    }
}


size_t MoleReactor::componentIndex(const string& nm) const
{
    size_t k = speciesIndex(nm);
    if (k != npos) {
        return k + m_sidx;
    } else if (nm == "int_energy") {
        return 0;
    } else if (nm == "volume") {
        return 1;
    } else {
        return npos;
    }
}

string MoleReactor::componentName(size_t k) {
    if (k == 0) {
        return "int_energy";
    } else if (k == 1) {
        return "volume";
    } else if (k >= m_sidx && k < neq()) {
        k -= m_sidx;
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
    throw CanteraError("MoleReactor::componentName", "Index is out of bounds.");
}

double MoleReactor::upperBound(size_t k) const {
    // Component is either int_energy, volume, or moles of a bulk or surface species
    return BigNumber;
}

double MoleReactor::lowerBound(size_t k) const {
    if (k == 0) {
        return -BigNumber; // int_energy
    } else if (k == 1) {
        return 0; // volume
    } else if (k >= 2 && k < m_nv) {
        return -Tiny; // moles of bulk or surface species
    } else {
        throw CanteraError("MoleReactor::lowerBound", "Index {} is out of bounds.", k);
    }
}

void MoleReactor::resetBadValues(double* y) {
    for (size_t k = m_sidx; k < m_nv; k++) {
        y[k] = std::max(y[k], 0.0);
    }
}

}
