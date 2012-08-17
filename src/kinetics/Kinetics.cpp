/**
 *  @file Kinetics.cpp
 *      Declarations for the base class for kinetics
 *    managers (see \ref  kineticsmgr and class
 *  \link Cantera::Kinetics Kinetics\endlink).
 *
 *      Kinetics managers calculate rates of progress of species due to homogeneous or heterogeneous kinetics.
 */
// Copyright 2001-2004  California Institute of Technology

// Why InterfaceKinetics.h and not Kinetics.h ??

#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/StoichManager.h"
#include "cantera/kinetics/RateCoeffMgr.h"

#include "ImplicitSurfChem.h"

#include <iostream>
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


//  Copy Constructor for the %Kinetics object.
/*
 * Currently, this is not fully implemented. If called it will
 * throw an exception.
 */
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

// Assignment operator
/*
 *  This is NOT a virtual function.
 *
 * @param right    Reference to %Kinetics object to be copied into the
 *                 current one.
 */
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

//====================================================================================================================
// Duplication routine for objects which inherit from
// Kinetics
/*
 *  This virtual routine can be used to duplicate %Kinetics objects
 *  inherited from %Kinetics even if the application only has
 *  a pointer to %Kinetics to work with.
 *
 *  These routines are basically wrappers around the derived copy
 *  constructor.
 */
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

//====================================================================================================================
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
//====================================================================================================================
/**
 * Takes as input an array of properties for all species in the
 * mechanism and copies those values belonging to a particular
 * phase to the output array.
 * @param data Input data array.
 * @param phase Pointer to one of the phase objects participating
 * in this reaction mechanism
 * @param phase_data Output array where the values for the the
 * specified phase are to be written.
 */
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


/**
 * kineticsSpeciesName():
 *
 * Return the string name of the kth species in the kinetics
 * manager. k is an integer from 0 to ktot - 1, where ktot is
 * the number of species in the kinetics manager, which is the
 * sum of the number of species in all phases participating in
 * the kinetics manager.  If k is out of bounds, the string
 * "<unknown>" is returned.
 */
string Kinetics::kineticsSpeciesName(size_t k) const
{
    for (size_t n = m_start.size()-1; n != npos; n--) {
        if (k >= m_start[n]) {
            return thermo(n).speciesName(k - m_start[n]);
        }
    }
    return "<unknown>";
}

/**
 * This routine will look up a species number based on the input
 * std::string nm. The lookup of species will occur for all phases
 * listed in the kinetics object.
 *
 *  return
 *   - If a match is found, the position in the species list is returned.
 *   - If no match is found, the value -1 is returned.
 *
 * @param nm   Input string name of the species
 */
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

/**
 * This routine will look up a species number based on the input
 * std::string nm. The lookup of species will occur in the specified
 * phase of the object, or all phases if ph is "<any>".
 *
 *  return
 *   - If a match is found, the position in the species list is returned.
 *   - If no match is found, the value npos (-1) is returned.
 *
 * @param nm   Input string name of the species
 * @param ph   Input string name of the phase.
 */
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


/**
 * This function looks up the string name of a species and
 * returns a reference to the ThermoPhase object of the
 * phase where the species resides.
 * Will throw an error if the species string doesn't match.
 */
thermo_t& Kinetics::speciesPhase(std::string nm)
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

//==============================================================================================
/*
 * This function takes as an argument the kineticsSpecies index
 * (i.e., the list index in the list of species in the kinetics
 * manager) and returns the index of the phase owning the
 * species.
 */
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

/*
 * Add a phase to the kinetics manager object. This must
 * be done before the function init() is called or
 * before any reactions are input.
 * The following fields are updated:
 *  m_start -> vector of integers, containing the
 *             starting position of the species for
 *             each phase in the kinetics mechanism.
 *  m_surfphase -> index of the surface phase.
 *  m_thermo -> vector of pointers to ThermoPhase phases
 *              that participate in the kinetics
 *              mechanism.
 *  m_phaseindex -> map containing the string id of each
 *              ThermoPhase phase as a key and the
 *              index of the phase within the kinetics
 *              manager object as the value.
 */
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

// Private function of the class Kinetics, indicating that a function
//  inherited from the base class hasn't had a definition assigned to it
/*
 * @param m String message
 */
void Kinetics::err(std::string m) const
{
    throw CanteraError("Kinetics::" + m,
                       "The default Base class method was called, when "
                       "the inherited class's method should "
                       "have been called");
}

}
