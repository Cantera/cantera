/**
 *  @file MixTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */
// copyright 2001 California Institute of Technology

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/MixTransport.h"

#include "cantera/base/utilities.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
MixTransport::MixTransport() :
    m_condcoeffs(0),
    m_cond(0),
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{
}

MixTransport::MixTransport(const MixTransport& right) :
    GasTransport(right),
    m_condcoeffs(0),
    m_cond(0),
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_debug(false)
{
    *this = right;
}

MixTransport&  MixTransport::operator=(const MixTransport& right)
{
    if (&right == this) {
        return *this;
    }
    GasTransport::operator=(right);

    m_condcoeffs = right.m_condcoeffs;
    m_cond = right.m_cond;
    m_lambda = right.m_lambda;
    m_spcond_ok = right.m_spcond_ok;
    m_condmix_ok = right.m_condmix_ok;
    m_debug = right.m_debug;

    return *this;
}

Transport* MixTransport::duplMyselfAsTransport() const
{
    return new MixTransport(*this);
}

bool MixTransport::initGas(GasTransportParams& tr)
{
    GasTransport::initGas(tr);

    m_eps = tr.eps;
    m_sigma = tr.sigma;
    m_alpha = tr.alpha;
    m_dipole = tr.dipole;
    m_zrot = tr.zrot;
    m_crot = tr.crot;

    // copy polynomials and parameters into local storage
    m_condcoeffs = tr.condcoeffs;

    m_cond.resize(m_nsp);

    // set flags all false
    m_spcond_ok = false;
    m_condmix_ok = false;

    return true;
}

void MixTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

doublereal MixTransport::thermalConductivity()
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

void MixTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void MixTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{
    update_T();
    update_C();

    getMixDiffCoeffs(DATA_PTR(m_spwork));

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();
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

void MixTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("MixTransport::update_T",
                           "negative temperature "+fp2str(t));
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be redone.
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}

void MixTransport::update_C()
{
    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_visc_ok = false;
    m_condmix_ok = false;

    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void MixTransport::updateCond_T()
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
