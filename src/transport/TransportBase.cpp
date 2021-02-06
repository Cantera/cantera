/**
 *  @file TransportBase.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/TransportBase.h"

using namespace std;

namespace Cantera
{
Transport::Transport(thermo_t* thermo, size_t ndim) :
    m_thermo(thermo),
    m_ready(false),
    m_nsp(0),
    m_nDim(ndim),
    m_velocityBasis(VB_MASSAVG)
{
}

bool Transport::ready()
{
    return m_ready;
}

void Transport::setNDim(const int ndim)
{
    m_nDim = ndim;
}

void Transport::checkSpeciesIndex(size_t k) const
{
    if (k >= m_nsp) {
        throw IndexError("Transport::checkSpeciesIndex", "species", k, m_nsp-1);
    }
}

void Transport::checkSpeciesArraySize(size_t kk) const
{
    if (m_nsp > kk) {
        throw ArraySizeError("Transport::checkSpeciesArraySize", kk, m_nsp);
    }
}

void Transport::setParameters(const int type, const int k,
                              const doublereal* const p)
{
    throw NotImplementedError("Transport::setParameters");
}

void Transport::setThermo(thermo_t& thermo)
{
    if (!ready()) {
        m_thermo = &thermo;
        m_nsp = m_thermo->nSpecies();
    } else  {
        size_t newNum = thermo.nSpecies();
        size_t oldNum = m_thermo->nSpecies();
        if (newNum != oldNum) {
            throw CanteraError("Transport::setThermo",
                               "base object cannot be changed after "
                               "the transport manager has been constructed because num species isn't the same.");
        }
        for (size_t i = 0; i < newNum; i++) {
            std::string newS0 = thermo.speciesName(i);
            std::string oldS0 = m_thermo->speciesName(i);
            if (newNum != oldNum) {
                throw CanteraError("Transport::setThermo",
                                   "base object cannot be changed after "
                                   "the transport manager has been constructed because species names are not the same");
            }
        }
        m_thermo = &thermo;
    }
}

void Transport::finalize()
{
    if (!ready()) {
        m_ready = true;
    } else {
        throw CanteraError("Transport::finalize",
                           "finalize has already been called.");
    }
}

void Transport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                 size_t ldx, const doublereal* const grad_X,
                                 size_t ldf, doublereal* const fluxes)
{
    throw NotImplementedError("Transport::getSpeciesFluxes");
}
}
