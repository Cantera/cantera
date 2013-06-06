/**
 *  @file Kinetics.cpp Declarations for the base class for kinetics managers
 *      (see \ref  kineticsmgr and class \link Cantera::Kinetics
 *      Kinetics\endlink).
 *
 *  Kinetics managers calculate rates of progress of species due to
 *  homogeneous or heterogeneous kinetics.
 */
// Copyright 2001-2004  California Institute of Technology

#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/StoichManager.h"
#include "cantera/kinetics/RateCoeffMgr.h"

#include "cantera/kinetics/ImplicitSurfChem.h"

using namespace std;


namespace Cantera
{
Kinetics::Kinetics() :
    m_ii(0),
    m_kk(0),
    m_perturb(0),
    m_reactants(0),
    m_products(0),
    m_thermo(0),
    m_start(0),
    m_phaseindex(),
    m_surfphase(npos),
    m_rxnphase(npos),
    m_mindim(4),
    m_dummygroups(0)
{
}

Kinetics::~Kinetics() {}

Kinetics::Kinetics(const Kinetics& right) :
    m_ii(0),
    m_kk(0),
    m_perturb(0),
    m_reactants(0),
    m_products(0),
    m_thermo(0),
    m_start(0),
    m_phaseindex(),
    m_surfphase(npos),
    m_rxnphase(npos),
    m_mindim(4),
    m_dummygroups(0)
{
    /*
     * Call the assignment operator
     */
    *this = right;
}

Kinetics& Kinetics::
operator=(const Kinetics& right)
{
    /*
     * Check for self assignment.
     */
    if (this == &right) {
        return *this;
    }

    m_ii                = right.m_ii;
    m_kk                = right.m_kk;
    m_perturb           = right.m_perturb;
    m_reactants         = right.m_reactants;
    m_products          = right.m_products;

    m_thermo            = right.m_thermo; //  DANGER -> shallow pointer copy

    m_start             = right.m_start;
    m_phaseindex        = right.m_phaseindex;
    m_surfphase         = right.m_surfphase;
    m_rxnphase          = right.m_rxnphase;
    m_mindim            = right.m_mindim;
    m_dummygroups       = right.m_dummygroups;

    return *this;
}

Kinetics* Kinetics::duplMyselfAsKinetics(const std::vector<thermo_t*> & tpVector) const
{
    Kinetics* ko = new Kinetics(*this);

    ko->assignShallowPointers(tpVector);
    return ko;
}

int Kinetics::type() const
{
    return 0;
}

void Kinetics::checkReactionIndex(size_t i) const
{
    if (i >= m_ii) {
        throw IndexError("checkReactionIndex", "reactions", i, m_ii-1);
    }
}

void Kinetics::checkReactionArraySize(size_t ii) const
{
    if (m_ii > ii) {
        throw ArraySizeError("checkReactionArraySize", ii, m_ii);
    }
}

void Kinetics::checkPhaseIndex(size_t m) const
{
    if (m >= nPhases()) {
        throw IndexError("checkPhaseIndex", "phase", m, nPhases()-1);
    }
}

void Kinetics::checkPhaseArraySize(size_t mm) const
{
    if (nPhases() > mm) {
        throw ArraySizeError("checkPhaseArraySize", mm, nPhases());
    }
}

void Kinetics::checkSpeciesIndex(size_t k) const
{
    if (k >= m_kk) {
        throw IndexError("checkSpeciesIndex", "species", k, m_kk-1);
    }
}

void Kinetics::checkSpeciesArraySize(size_t kk) const
{
    if (m_kk > kk) {
        throw ArraySizeError("checkSpeciesArraySize", kk, m_kk);
    }
}

void Kinetics::assignShallowPointers(const std::vector<thermo_t*> & tpVector)
{
    size_t ns = tpVector.size();
    if (ns != m_thermo.size()) {
        throw CanteraError(" Kinetics::assignShallowPointers",
                           " Number of ThermoPhase objects arent't the same");
    }
    for (size_t i = 0; i < ns; i++) {
        ThermoPhase* ntp = tpVector[i];
        ThermoPhase* otp = m_thermo[i];
        if (ntp->id() != otp->id()) {
            throw CanteraError(" Kinetics::assignShallowPointers",
                               " id() of the ThermoPhase objects isn't the same");
        }
        if (ntp->eosType() != otp->eosType()) {
            throw CanteraError(" Kinetics::assignShallowPointers",
                               " eosType() of the ThermoPhase objects isn't the same");
        }
        if (ntp->nSpecies() != otp->nSpecies()) {
            throw CanteraError(" Kinetics::assignShallowPointers",
                               " Number of ThermoPhase objects isn't the same");
        }
        m_thermo[i] = tpVector[i];
    }


}

void Kinetics::selectPhase(const doublereal* data, const thermo_t* phase,
                           doublereal* phase_data)
{
    for (size_t n = 0; n < nPhases(); n++) {
        if (phase == m_thermo[n]) {
            size_t nsp = phase->nSpecies();
            copy(data + m_start[n],
                 data + m_start[n] + nsp, phase_data);
            return;
        }
    }
    throw CanteraError("Kinetics::selectPhase", "Phase not found.");
}

string Kinetics::kineticsSpeciesName(size_t k) const
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return thermo(n).speciesName(k - m_start[n]);
        }
    }
    return "<unknown>";
}

size_t Kinetics::kineticsSpeciesIndex(const std::string& nm) const
{
    for (size_t n = 0; n < m_thermo.size(); n++) {
        string id = thermo(n).id();
        // Check the ThermoPhase object for a match
        size_t k = thermo(n).speciesIndex(nm);
        if (k != npos) {
            return k + m_start[n];
        }
    }
    return npos;
}

size_t Kinetics::kineticsSpeciesIndex(const std::string& nm,
                                      const std::string& ph) const
{
    if (ph == "<any>") {
        return kineticsSpeciesIndex(nm);
    }

    for (size_t n = 0; n < m_thermo.size(); n++) {
        string id = thermo(n).id();
        if (ph == id) {
            size_t k = thermo(n).speciesIndex(nm);
            if (k == npos) {
                return npos;
            }
            return k + m_start[n];
        }
    }
    return npos;
}

thermo_t& Kinetics::speciesPhase(const std::string& nm)
{
    size_t np = m_thermo.size();
    size_t k;
    string id;
    for (size_t n = 0; n < np; n++) {
        k = thermo(n).speciesIndex(nm);
        if (k != npos) {
            return thermo(n);
        }
    }
    throw CanteraError("speciesPhase", "unknown species "+nm);
    return thermo(0);
}

size_t Kinetics::speciesPhaseIndex(size_t k)
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return n;
        }
    }
    throw CanteraError("speciesPhaseIndex", "illegal species index: "+int2str(k));
    return npos;
}

void Kinetics::addPhase(thermo_t& thermo)
{
    // if not the first thermo object, set the start position
    // to that of the last object added + the number of its species
    if (m_thermo.size() > 0) {
        m_start.push_back(m_start.back()
                          + m_thermo.back()->nSpecies());
    }
    // otherwise start at 0
    else {
        m_start.push_back(0);
    }

    // the phase with lowest dimensionality is assumed to be the
    // phase/interface at which reactions take place
    if (thermo.nDim() <= m_mindim) {
        m_mindim = thermo.nDim();
        m_rxnphase = nPhases();
    }

    // there should only be one surface phase
    int ptype = -100;
    if (type() == cEdgeKinetics) {
        ptype = cEdge;
    } else if (type() == cInterfaceKinetics) {
        ptype = cSurf;
    }
    if (thermo.eosType() == ptype) {
        m_surfphase = nPhases();
        m_rxnphase = nPhases();
    }
    m_thermo.push_back(&thermo);
    m_phaseindex[m_thermo.back()->id()] = nPhases();
}

void Kinetics::finalize()
{
    m_kk = 0;
    for (size_t n = 0; n < nPhases(); n++) {
        size_t nsp = m_thermo[n]->nSpecies();
        m_kk += nsp;
    }
}

void Kinetics::err(const std::string& m) const
{
    throw CanteraError("Kinetics::" + m,
                       "The default Base class method was called, when "
                       "the inherited class's method should "
                       "have been called");
}

}
