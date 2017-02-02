/**
 *  @file UlnTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/UlnTransport.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
UlnTransport::UlnTransport() :
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{
}

UlnTransport::UlnTransport(const UlnTransport& right) :
    GasTransport(right),
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{
    *this = right;
}

UlnTransport& UlnTransport::operator=(const UlnTransport& right)
{
    if (&right == this) {
        return *this;
    }
    GasTransport::operator=(right);

    m_cond = right.m_cond;
    m_lambda = right.m_lambda;
    m_spcond_ok = right.m_spcond_ok;
    m_condmix_ok = right.m_condmix_ok;
    m_debug = right.m_debug;

    return *this;
}

Transport* UlnTransport::duplMyselfAsTransport() const
{
    return new UlnTransport(*this);
}

void UlnTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    GasTransport::init(thermo, mode, log_level);
    m_cond.resize(m_nsp);

    // set flags all false
    m_spcond_ok = false;
    m_condmix_ok = false;
}

void UlnTransport::getMobilities(doublereal* const mobil)
{
	 getUlnDiffCoeffs(m_spwork.data(),thermalConductivity());
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

doublereal UlnTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_spcond_ok) {
        updateCond_T();
    }
    if (!m_condmix_ok) {
        doublereal sum1 = 0.0, sum2 = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum1 += m_molefracs[k] * m_cond[k];
            sum2 += m_molefracs[k] / m_cond[k];
        }
        m_lambda = 0.5*(sum1 + 1.0/sum2);
        m_condmix_ok = true;
    }
    return m_lambda;
}

void UlnTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void UlnTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{
    update_T();
    update_C();
	 getUlnDiffCoeffs(m_spwork.data(),thermalConductivity());
    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();
    vector_fp sum(ndim,0.0);
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * grad_X[n*ldx + k];
            sum[n] += fluxes[n*ldf + k];
        }
    }
    // add correction flux to enforce sum to zero
    for (size_t n = 0; n < ndim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] -= y[k]*sum[n];
        }
    }
}

void UlnTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp && m_nsp == m_thermo->nSpecies()) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("UlnTransport::update_T",
                           "negative temperature {}", t);
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be redone.
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}

void UlnTransport::update_C()
{
    // signal that concentration-dependent quantities will need to be recomputed
    // before use, and update the local mole fractions.
    m_visc_ok = false;
    m_condmix_ok = false;
    m_thermo->getMoleFractions(m_molefracs.data());

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void UlnTransport::updateCond_T()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = m_sqrt_t * dot5(m_polytempvec, m_condcoeffs[k]);
        }
    }
    m_spcond_ok = true;
    m_condmix_ok = false;
}

}
