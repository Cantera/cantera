/**
 *  @file UnityLewisNumberTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures with
 *  diffusion coefficients set based on the assumption that Lewis #  = 1.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/UnityLewisNumberTransport.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{
UnityLewisTransport::UnityLewisTransport() :
    m_lambda(0.0),
    m_spcond_ok(false),
    m_condmix_ok(false)
{
}

void UnityLewisTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    GasTransport::init(thermo, mode, log_level);
    m_cond.resize(m_nsp);

    // set flags all false
    m_spcond_ok = false;
    m_condmix_ok = false;
}


void UnityLewisTransport::getMixDiffCoeffs(doublereal* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            double sum2 = 0.0;
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != k) {
                    sum2 += m_molefracs[j] / m_bdiff(j,k);
                }
            }
            if (sum2 <= 0.0) {
                d[k] = m_bdiff(k,k) / p;
            } else {
                d[k] = (mmw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
            }
        }
    }
}


void MixTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                    size_t ldx, const doublereal* const grad_X,
                                    size_t ldf, doublereal* const fluxes)
{
    update_T();
    update_C();
    getMixDiffCoeffs(m_spwork.data());
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

void MixTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp && m_nsp == m_thermo->nSpecies()) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("MixTransport::update_T",
                           "negative temperature {}", t);
    }
    GasTransport::update_T();
    // temperature has changed, so polynomial fits will need to be redone.
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}

void MixTransport::update_C()
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
