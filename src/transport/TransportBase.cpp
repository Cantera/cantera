/**
 *  @file TransportBase.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/LiquidTransport.h"
#include "cantera/base/ctexceptions.h"

#include "cantera/base/utilities.h"
#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"

#include "cantera/numerics/ctlapack.h"

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

Transport::Transport(const Transport& right)
{
    m_thermo        = right.m_thermo;
    m_ready         = right.m_ready;
    m_nsp          = right.m_nsp;
    m_nDim          = right.m_nDim;
    m_velocityBasis = right.m_velocityBasis;
}


Transport& Transport::operator=(const Transport& right)
{
    if (&right != this) {
        return *this;
    }
    m_thermo        = right.m_thermo;
    m_ready         = right.m_ready;
    m_nsp          = right.m_nsp;
    m_nDim          = right.m_nDim;
    m_velocityBasis = right.m_velocityBasis;
    return *this;
}

Transport* Transport::duplMyselfAsTransport() const
{
    return new Transport(*this);
}

Transport::~Transport()
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
        throw IndexError("checkSpeciesIndex", "species", k, m_nsp-1);
    }
}

void Transport::checkSpeciesArraySize(size_t kk) const
{
    if (m_nsp > kk) {
        throw ArraySizeError("checkSpeciesArraySize", kk, m_nsp);
    }
}

void Transport::setParameters(const int type, const int k,
                              const doublereal* const p)
{
    err("setParameters");
}

void Transport::setThermo(thermo_t& thermo)
{
    if (!ready()) {
        m_thermo = &thermo;
        m_nsp = m_thermo->nSpecies();
    } else  {
        int newNum = thermo.nSpecies();
        int oldNum = m_thermo->nSpecies();
        if (newNum != oldNum) {
            throw CanteraError("Transport::setThermo",
                               "base object cannot be changed after "
                               "the transport manager has been constructed because num species isn't the same.");
        }
        for (int i = 0; i < newNum; i++) {
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

doublereal Transport::err(const std::string& msg) const
{

    throw CanteraError("Transport Base Class",
                       "\n\n\n**** Method "+ msg +" not implemented in model "
                       + int2str(model()) + " ****\n"
                       "(Did you forget to specify a transport model?)\n\n\n");

    return 0.0;
}

void Transport::finalize()
{
    if (!ready()) {
        m_ready = true;
    } else
        throw CanteraError("Transport::finalize",
                           "finalize has already been called.");
}

void Transport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                 size_t ldx, const doublereal* const grad_X,
                                 size_t ldf, doublereal* const fluxes)
{
    err("getSpeciesFluxes");
}
}
