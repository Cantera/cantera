/**
 *  @file ImplicitSurfChem.cpp
 * Definitions for the implicit integration of surface site density equations
 *  (see \ref  kineticsmgr and class
 *  \link Cantera::ImplicitSurfChem ImplicitSurfChem\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/ImplicitSurfChem.h"
#include "cantera/kinetics/solveSP.h"
#include "cantera/thermo/SurfPhase.h"

using namespace std;

namespace Cantera
{

ImplicitSurfChem::ImplicitSurfChem(
        vector<InterfaceKinetics*> k, double rtol, double atol,
        double maxStepSize, size_t maxSteps,
        size_t maxErrTestFails) :
    m_nv(0),
    m_numTotalBulkSpecies(0),
    m_numTotalSpecies(0),
    m_atol(atol),
    m_rtol(rtol),
    m_maxstep(maxStepSize),
    m_nmax(maxSteps),
    m_maxErrTestFails(maxErrTestFails),
    m_mediumSpeciesStart(-1),
    m_bulkSpeciesStart(-1),
    m_surfSpeciesStart(-1),
    m_commonTempPressForPhases(true),
    m_ioFlag(0)
{
    size_t ntmax = 0;
    size_t kinSpIndex = 0;
    // Loop over the number of surface kinetics objects
    for (size_t n = 0; n < k.size(); n++) {
        InterfaceKinetics* kinPtr = k[n];
        m_vecKinPtrs.push_back(kinPtr);
        size_t ns = k[n]->surfacePhaseIndex();
        if (ns == npos) {
            throw CanteraError("ImplicitSurfChem::ImplicitSurfChem",
                               "kinetics manager contains no surface phase");
        }
        m_surfindex.push_back(ns);
        m_surf.push_back((SurfPhase*)&k[n]->thermo(ns));
        size_t nsp = m_surf.back()->nSpecies();
        m_nsp.push_back(nsp);
        m_nv += m_nsp.back();
        size_t nt = k[n]->nTotalSpecies();
        ntmax = std::max(nt, ntmax);
        m_specStartIndex.push_back(kinSpIndex);
        kinSpIndex += nsp;
        size_t nPhases = kinPtr->nPhases();
        vector_int pLocTmp(nPhases);
        size_t imatch = npos;
        for (size_t ip = 0; ip < nPhases; ip++) {
            if (ip != ns) {
                ThermoPhase* thPtr = & kinPtr->thermo(ip);
                if ((imatch = checkMatch(m_bulkPhases, thPtr)) == npos) {
                    m_bulkPhases.push_back(thPtr);
                    nsp = thPtr->nSpecies();
                    m_numTotalBulkSpecies += nsp;
                    imatch = m_bulkPhases.size() - 1;
                }
                pLocTmp[ip] = int(imatch);
            } else {
                pLocTmp[ip] = -int(n);
            }
        }
        pLocVec.push_back(pLocTmp);
    }
    m_numTotalSpecies = m_nv + m_numTotalBulkSpecies;
    m_concSpecies.resize(m_numTotalSpecies, 0.0);
    m_concSpeciesSave.resize(m_numTotalSpecies, 0.0);

    m_integ.reset(newIntegrator("CVODE"));

    // use backward differencing, with a full Jacobian computed
    // numerically, and use a Newton linear iterator
    m_integ->setMethod(BDF_Method);
    m_integ->setProblemType(DENSE + NOJAC);
    m_work.resize(ntmax);
}

int ImplicitSurfChem::checkMatch(std::vector<ThermoPhase*> m_vec, ThermoPhase* thPtr)
{
    int retn = -1;
    for (int i = 0; i < (int) m_vec.size(); i++) {
        ThermoPhase* th = m_vec[i];
        if (th == thPtr) {
            return i;
        }
    }
    return retn;
}

void ImplicitSurfChem::getState(doublereal* c)
{
    size_t loc = 0;
    for (size_t n = 0; n < m_surf.size(); n++) {
        m_surf[n]->getCoverages(c + loc);
        loc += m_nsp[n];
    }
}

void ImplicitSurfChem::setMaxStepSize(double maxstep)
{
    m_maxstep = maxstep;
    if (m_maxstep > 0) {
        m_integ->setMaxStepSize(m_maxstep);
    }
}

void ImplicitSurfChem::setTolerances(double rtol, double atol)
{
    m_rtol = rtol;
    m_atol = atol;
    m_integ->setTolerances(m_rtol, m_atol);
}

void ImplicitSurfChem::setMaxSteps(size_t maxsteps)
{
    m_nmax = maxsteps;
    m_integ->setMaxSteps(static_cast<int>(m_nmax));
}

void ImplicitSurfChem::setMaxErrTestFails(size_t maxErrTestFails)
{
    m_maxErrTestFails = maxErrTestFails;
    m_integ->setMaxErrTestFails(static_cast<int>(m_maxErrTestFails));
}

void ImplicitSurfChem::initialize(doublereal t0)
{
    this->setTolerances(m_rtol, m_atol);
    this->setMaxStepSize(m_maxstep);
    this->setMaxSteps(m_nmax);
    this->setMaxErrTestFails(m_maxErrTestFails);
    m_integ->initialize(t0, *this);
}

void ImplicitSurfChem::integrate(doublereal t0, doublereal t1)
{
    this->initialize(t0);
    if (fabs(t1 - t0) < m_maxstep || m_maxstep == 0) {
        // limit max step size on this run to t1 - t0
        m_integ->setMaxStepSize(t1 - t0);
    }
    m_integ->integrate(t1);
    updateState(m_integ->solution());
}

void ImplicitSurfChem::integrate0(doublereal t0, doublereal t1)
{
    m_integ->integrate(t1);
    updateState(m_integ->solution());
}

void ImplicitSurfChem::updateState(doublereal* c)
{
    size_t loc = 0;
    for (size_t n = 0; n < m_surf.size(); n++) {
        m_surf[n]->setCoverages(c + loc);
        loc += m_nsp[n];
    }
}

void ImplicitSurfChem::eval(doublereal time, doublereal* y,
                            doublereal* ydot, doublereal* p)
{
    updateState(y); // synchronize the surface state(s) with y
    size_t loc = 0;
    for (size_t n = 0; n < m_surf.size(); n++) {
        double rs0 = 1.0/m_surf[n]->siteDensity();
        m_vecKinPtrs[n]->getNetProductionRates(m_work.data());
        size_t kstart = m_vecKinPtrs[n]->kineticsSpeciesIndex(0,m_surfindex[n]);
        double sum = 0.0;
        for (size_t k = 1; k < m_nsp[n]; k++) {
            ydot[k + loc] = m_work[kstart + k] * rs0 * m_surf[n]->size(k);
            sum -= ydot[k];
        }
        ydot[loc] = sum;
        loc += m_nsp[n];
    }
}

void ImplicitSurfChem::solvePseudoSteadyStateProblem(int ifuncOverride,
        doublereal timeScaleOverride)
{
    int ifunc;
    // set bulkFunc. We assume that the bulk concentrations are constant.
    int bulkFunc = BULK_ETCH;
    // time scale - time over which to integrate equations
    doublereal time_scale = timeScaleOverride;
    if (!m_surfSolver) {
        m_surfSolver.reset(new solveSP(this, bulkFunc));
        // set ifunc, which sets the algorithm.
        ifunc = SFLUX_INITIALIZE;
    } else {
        ifunc = SFLUX_RESIDUAL;
    }

    // Possibly override the ifunc value
    if (ifuncOverride >= 0) {
        ifunc = ifuncOverride;
    }

    // Get the specifications for the problem from the values
    // in the ThermoPhase objects for all phases.
    //
    //  1) concentrations of all species in all phases, m_concSpecies[]
    //  2) Temperature and pressure
    getConcSpecies(m_concSpecies.data());
    InterfaceKinetics* ik = m_vecKinPtrs[0];
    ThermoPhase& tp = ik->thermo(ik->reactionPhaseIndex());
    doublereal TKelvin = tp.temperature();
    doublereal PGas = tp.pressure();

    // Make sure that there is a common temperature and pressure for all
    // ThermoPhase objects belonging to the interfacial kinetics object, if it
    // is required by the problem statement.
    if (m_commonTempPressForPhases) {
        setCommonState_TP(TKelvin, PGas);
    }

    doublereal reltol = 1.0E-6;
    doublereal atol = 1.0E-20;

    // Install a filter for negative concentrations. One of the few ways solveSS
    // can fail is if concentrations on input are below zero.
    bool rset = false;
    for (size_t k = 0; k < m_nv; k++) {
        if (m_concSpecies[k] < 0.0) {
            rset = true;
            m_concSpecies[k] = 0.0;
        }
    }
    if (rset) {
        setConcSpecies(m_concSpecies.data());
    }

    m_surfSolver->m_ioflag = m_ioFlag;

    // Save the current solution
    m_concSpeciesSave = m_concSpecies;

    int retn = m_surfSolver->solveSurfProb(ifunc, time_scale, TKelvin, PGas,
                                           reltol, atol);
    if (retn != 1) {
        // reset the concentrations
        m_concSpecies = m_concSpeciesSave;
        setConcSpecies(m_concSpeciesSave.data());
        ifunc = SFLUX_INITIALIZE;
        retn = m_surfSolver->solveSurfProb(ifunc, time_scale, TKelvin, PGas,
                                           reltol, atol);
        if (retn != 1) {
            throw CanteraError("ImplicitSurfChem::solvePseudoSteadyStateProblem",
                               "solveSP return an error condition!");
        }
    }
}

void ImplicitSurfChem::getConcSpecies(doublereal* const vecConcSpecies) const
{
    size_t kstart;
    for (size_t ip = 0; ip < m_surf.size(); ip++) {
        ThermoPhase* TP_ptr = m_surf[ip];
        kstart = m_specStartIndex[ip];
        TP_ptr->getConcentrations(vecConcSpecies + kstart);
    }
    kstart = m_nv;
    for (size_t ip = 0; ip < m_bulkPhases.size(); ip++) {
        ThermoPhase* TP_ptr = m_bulkPhases[ip];
        TP_ptr->getConcentrations(vecConcSpecies + kstart);
        kstart += TP_ptr->nSpecies();
    }
}

void ImplicitSurfChem::setConcSpecies(const doublereal* const vecConcSpecies)
{
    size_t kstart;
    for (size_t ip = 0; ip < m_surf.size(); ip++) {
        ThermoPhase* TP_ptr = m_surf[ip];
        kstart = m_specStartIndex[ip];
        TP_ptr->setConcentrations(vecConcSpecies + kstart);
    }
    kstart = m_nv;
    for (size_t ip = 0; ip < m_bulkPhases.size(); ip++) {
        ThermoPhase* TP_ptr = m_bulkPhases[ip];
        TP_ptr->setConcentrations(vecConcSpecies + kstart);
        kstart += TP_ptr->nSpecies();
    }
}

void ImplicitSurfChem::setCommonState_TP(doublereal TKelvin, doublereal PresPa)
{
    for (size_t ip = 0; ip < m_surf.size(); ip++) {
        ThermoPhase* TP_ptr = m_surf[ip];
        TP_ptr->setState_TP(TKelvin, PresPa);
    }
    for (size_t ip = 0; ip < m_bulkPhases.size(); ip++) {
        ThermoPhase* TP_ptr = m_bulkPhases[ip];
        TP_ptr->setState_TP(TKelvin, PresPa);
    }
}

}
