/**
 *  @file MixTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/MixTransport.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

void MixTransport::init(shared_ptr<ThermoPhase> thermo, int mode)
{
    GasTransport::init(thermo, mode);
    m_cond.resize(m_nsp);
}

void MixTransport::getMobilities(span<double> mobil)
{
    checkArraySize("MixTransport::getMobilities", mobil.size(), m_nsp);
    getMixDiffCoeffs(m_spwork);
    double c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

double MixTransport::thermalConductivity()
{
    update_T();
    update_C();
    if (!m_spcond_ok) {
        updateCond_T();
    }
    if (!m_condmix_ok) {
        double sum1 = 0.0, sum2 = 0.0;
        for (size_t k = 0; k < m_nsp; k++) {
            sum1 += m_molefracs[k] * m_cond[k];
            sum2 += m_molefracs[k] / m_cond[k];
        }
        m_lambda = 0.5*(sum1 + 1.0/sum2);
        m_condmix_ok = true;
    }
    return m_lambda;
}

void MixTransport::getThermalDiffCoeffs(span<double> dt)
{
    checkArraySize("MixTransport::getThermalDiffCoeffs", dt.size(), m_nsp);
    update_T();
    update_C();
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    if (!m_viscwt_ok) {
        updateViscosity_T();
    }

    auto y = m_thermo->massFractions();

    vector<double>& a = m_spwork;

    for (size_t k=0; k<m_nsp; ++k) {
        dt[k] = 0.;

        if (y[k] < Tiny) {
            a[k] = 0.;
            continue;
        }

        double lambda_mono_k = (15./4.) * m_visc[k] / m_mw[k];

        double sum = 0.;
        for (size_t j=0; j<m_nsp; ++j) {
            if (j != k) {
                sum += m_molefracs[j]*m_phi(k,j);
            }
        };

        a[k] = lambda_mono_k / (1. + 1.065 * sum / m_molefracs[k]);
    }

    double rp = 1./m_thermo->pressure();

    for (size_t k=0; k<m_nsp-1; ++k) {
        for (size_t j=k+1; j<m_nsp; ++j) {

            double log_tstar = std::log(m_kbt/m_epsilon(k,j));

            int ipoly = m_poly[k][j];

            double Cstar = 0.;
            if (m_mode == CK_Mode) {
                Cstar = poly6(log_tstar, m_cstar_poly[ipoly]);
            } else {
                Cstar = poly8(log_tstar, m_cstar_poly[ipoly]);
            }

            double dt_T = ((1.2*Cstar - 1.0)/(m_bdiff(k,j)*rp))
                          / (m_mw[k] + m_mw[j]);

            dt[k] += dt_T * (y[k]*a[j] - y[j]*a[k]);
            dt[j] += dt_T * (y[j]*a[k] - y[k]*a[j]);
        }
    }

    vector<double>& Dm = m_spwork;
    getMixDiffCoeffs(Dm);

    double mmw = m_thermo->meanMolecularWeight();
    double norm = 0.;
    for (size_t k=0; k<m_nsp; ++k) {
        dt[k] *= Dm[k] * m_mw[k] * mmw;
        norm += dt[k];
    }

    // ensure that the sum of all Soret diffusion coefficients is zero
    for (size_t k=0; k<m_nsp; ++k) {
        dt[k] -= y[k]*norm;
    }
}

void MixTransport::getSpeciesFluxes(size_t ndim, span<const double> grad_T,
                                    size_t ldx, span<const double> grad_X,
                                    size_t ldf, span<double> fluxes)
{
    checkArraySize("MixTransport::getSpeciesFluxes: grad_T", grad_T.size(), ndim);
    if (ldx < m_nsp) {
        throw CanteraError("MixTransport::getSpeciesFluxes", "ldx is too small");
    }
    checkArraySize("MixTransport::getSpeciesFluxes: grad_X", grad_X.size(), ldx * ndim);
    if (ldf < m_nsp) {
        throw CanteraError("MixTransport::getSpeciesFluxes", "ldf is too small");
    }
    checkArraySize("MixTransport::getSpeciesFluxes: fluxes", fluxes.size(), ldf * ndim);
    update_T();
    update_C();
    getMixDiffCoeffs(m_spwork);
    auto mw = m_thermo->molecularWeights();
    auto y = m_thermo->massFractions();
    double rhon = m_thermo->molarDensity();
    vector<double> sum(ndim,0.0);
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
    double t = m_thermo->temperature();
    if (t == m_temp && m_nsp == m_thermo->nSpecies()) {
        return;
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
    m_thermo->getMoleFractions(m_molefracs);

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
