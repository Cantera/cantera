//! @file MoleReactor.cpp A zero-dimensional reactor with a moles as the
//! state

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/MoleReactor.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/zeroD/ReactorSurface.h"
#include "cantera/kinetics/Kinetics.h"

using namespace std;

namespace Cantera
{

void MoleReactor::getSurfaceInitialConditions(double* y)
{
    size_t loc = 0;
    for (auto& S : m_surfaces) {
        double area = S->area();
        auto currPhase = S->thermo();
        double tempLoc = currPhase->nSpecies();
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
    vector_fp coverages(m_nv_surf, 0.0);
    for (auto& S : m_surfaces) {
        auto surf = S->thermo();
        double invArea = 1/S->area();
        double invSurfDensity = 1/surf->siteDensity();
        double tempLoc = surf->nSpecies();
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
        S->syncState();
        kin->getNetProductionRates(&m_work[0]);
        size_t ns = kin->surfacePhaseIndex();
        size_t surfloc = kin->kineticsSpeciesIndex(0,ns);
        for (size_t k = 0; k < nk; k++) {
            RHS[loc + k] = m_work[surfloc + k] * wallarea / surf->size(k);
        }
        loc += nk;

        size_t bulkloc = kin->kineticsSpeciesIndex(m_thermo->speciesName(0));

        for (size_t k = 0; k < m_nsp; k++) {
            sdot[k] += m_work[bulkloc + k] * wallarea;
        }
    }
}

}
