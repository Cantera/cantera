//! @file GasTransport.cpp
#include "cantera/transport/GasTransport.h"
#include "cantera/transport/TransportParams.h"

namespace Cantera
{

GasTransport::GasTransport(ThermoPhase* thermo) :
    Transport(thermo),
    m_molefracs(0),
    m_viscmix(0.0),
    m_visc_ok(false),
    m_viscwt_ok(false),
    m_spvisc_ok(false),
    m_bindiff_ok(false),
    m_mode(0),
    m_phi(0,0),
    m_spwork(0),
    m_visc(0),
    m_visccoeffs(0),
    m_mw(0),
    m_wratjk(0,0),
    m_wratkj1(0,0),
    m_sqvisc(0),
    m_polytempvec(5),
    m_temp(-1.0),
    m_kbt(0.0),
    m_sqrt_kbt(0.0),
    m_sqrt_t(0.0),
    m_logt(0.0),
    m_t14(0.0),
    m_t32(0.0),
    m_diffcoeffs(0),
    m_bdiff(0, 0)
{
}

GasTransport::GasTransport(const GasTransport& right) :
    m_molefracs(0),
    m_viscmix(0.0),
    m_visc_ok(false),
    m_viscwt_ok(false),
    m_spvisc_ok(false),
    m_bindiff_ok(false),
    m_mode(0),
    m_phi(0,0),
    m_spwork(0),
    m_visc(0),
    m_visccoeffs(0),
    m_mw(0),
    m_wratjk(0,0),
    m_wratkj1(0,0),
    m_sqvisc(0),
    m_polytempvec(5),
    m_temp(-1.0),
    m_kbt(0.0),
    m_sqrt_kbt(0.0),
    m_sqrt_t(0.0),
    m_logt(0.0),
    m_t14(0.0),
    m_t32(0.0),
    m_diffcoeffs(0),
    m_bdiff(0, 0)
{
}

GasTransport& GasTransport::operator=(const GasTransport& right)
{
    m_molefracs = right.m_molefracs;
    m_viscmix = right.m_viscmix;
    m_visc_ok = right.m_visc_ok;
    m_viscwt_ok = right.m_viscwt_ok;
    m_spvisc_ok = right.m_spvisc_ok;
    m_bindiff_ok = right.m_bindiff_ok;
    m_mode = right.m_mode;
    m_phi = right.m_phi;
    m_spwork = right.m_spwork;
    m_visc = right.m_visc;
    m_mw = right.m_mw;
    m_wratjk = right.m_wratjk;
    m_wratkj1 = right.m_wratkj1;
    m_sqvisc = right.m_sqvisc;
    m_polytempvec = right.m_polytempvec;
    m_temp = right.m_temp;
    m_kbt = right.m_kbt;
    m_sqrt_kbt = right.m_sqrt_kbt;
    m_sqrt_t = right.m_sqrt_t;
    m_logt = right.m_logt;
    m_t14 = right.m_t14;
    m_t32 = right.m_t32;
    m_diffcoeffs = right.m_diffcoeffs;
    m_bdiff = right.m_bdiff;

    return *this;
}

bool GasTransport::initGas(GasTransportParams& tr)
{
    // constant mixture attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();

    // copy polynomials and parameters into local storage
    m_visccoeffs = tr.visccoeffs;
    m_diffcoeffs = tr.diffcoeffs;
    m_mode = tr.mode_;

    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);
    m_visc.resize(m_nsp);
    m_phi.resize(m_nsp, m_nsp, 0.0);
    m_bdiff.resize(m_nsp, m_nsp);

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(),
         m_thermo->molecularWeights().end(), m_mw.begin());

    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t k = j; k < m_nsp; k++) {
            m_wratjk(j,k) = sqrt(m_mw[j]/m_mw[k]);
            m_wratjk(k,j) = sqrt(m_wratjk(j,k));
            m_wratkj1(j,k) = sqrt(1.0 + m_mw[k]/m_mw[j]);
        }
    }

    m_sqvisc.resize(m_nsp);

    // set flags all false
    m_visc_ok = false;
    m_viscwt_ok = false;
    m_spvisc_ok = false;
    m_bindiff_ok = false;

    return true;
}

void GasTransport::update_T(void)
{
    double T = m_thermo->temperature();
    if (T == m_temp) {
        return;
    }

    m_temp = T;
    m_kbt = Boltzmann * m_temp;
    m_sqrt_kbt = sqrt(Boltzmann*m_temp);
    m_logt = log(m_temp);
    m_sqrt_t = sqrt(m_temp);
    m_t14 = sqrt(m_sqrt_t);
    m_t32 = m_temp * m_sqrt_t;

    // compute powers of log(T)
    m_polytempvec[0] = 1.0;
    m_polytempvec[1] = m_logt;
    m_polytempvec[2] = m_logt*m_logt;
    m_polytempvec[3] = m_logt*m_logt*m_logt;
    m_polytempvec[4] = m_logt*m_logt*m_logt*m_logt;

    // temperature has changed, so polynomial fits will need to be redone
    m_visc_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_bindiff_ok = false;
}

doublereal GasTransport::viscosity()
{
    update_T();
    update_C();

    if (m_visc_ok) {
        return m_viscmix;
    }

    doublereal vismix = 0.0;
    // update m_visc and m_phi if necessary
    if (!m_viscwt_ok) {
        updateViscosity_T();
    }

    multiply(m_phi, DATA_PTR(m_molefracs), DATA_PTR(m_spwork));

    for (size_t k = 0; k < m_nsp; k++) {
        vismix += m_molefracs[k] * m_visc[k]/m_spwork[k]; //denom;
    }
    m_viscmix = vismix;
    return vismix;
}

void GasTransport::updateViscosity_T()
{
    doublereal vratiokj, wratiojk, factor1;

    if (!m_spvisc_ok) {
        updateSpeciesViscosities();
    }

    // see Eq. (9-5.15) of Reid, Prausnitz, and Poling
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t k = j; k < m_nsp; k++) {
            vratiokj = m_visc[k]/m_visc[j];
            wratiojk = m_mw[j]/m_mw[k];

            // Note that m_wratjk(k,j) holds the square root of m_wratjk(j,k)!
            factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
            m_phi(k,j) = factor1*factor1 / (SqrtEight * m_wratkj1(j,k));
            m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
        }
    }
    m_viscwt_ok = true;
}

void GasTransport::updateSpeciesViscosities()
{
    update_T();
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_visc[k] = exp(dot4(m_polytempvec, m_visccoeffs[k]));
            m_sqvisc[k] = sqrt(m_visc[k]);
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            // the polynomial fit is done for sqrt(visc/sqrt(T))
            m_sqvisc[k] = m_t14 * dot5(m_polytempvec, m_visccoeffs[k]);
            m_visc[k] = (m_sqvisc[k] * m_sqvisc[k]);
        }
    }
    m_spvisc_ok = true;
}

void GasTransport::updateDiff_T()
{
    update_T();
    // evaluate binary diffusion coefficients at unit pressure
    size_t ic = 0;
    if (m_mode == CK_Mode) {
        for (size_t i = 0; i < m_nsp; i++) {
            for (size_t j = i; j < m_nsp; j++) {
                m_bdiff(i,j) = exp(dot4(m_polytempvec, m_diffcoeffs[ic]));
                m_bdiff(j,i) = m_bdiff(i,j);
                ic++;
            }
        }
    } else {
        for (size_t i = 0; i < m_nsp; i++) {
            for (size_t j = i; j < m_nsp; j++) {
                m_bdiff(i,j) = m_temp * m_sqrt_t*dot5(m_polytempvec,
                                                      m_diffcoeffs[ic]);
                m_bdiff(j,i) = m_bdiff(i,j);
                ic++;
            }
        }
    }
    m_bindiff_ok = true;
}

void GasTransport::getBinaryDiffCoeffs(const size_t ld, doublereal* const d)
{
    update_T();
    // if necessary, evaluate the binary diffusion coefficients from the polynomial fits
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    if (ld < m_nsp) {
        throw CanteraError(" MixTransport::getBinaryDiffCoeffs()", "ld is too small");
    }
    doublereal rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < m_nsp; i++)
        for (size_t j = 0; j < m_nsp; j++) {
            d[ld*j + i] = rp * m_bdiff(i,j);
        }
}

void GasTransport::getMixDiffCoeffs(doublereal* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal sumxw = 0.0;
    doublereal p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            sumxw += m_molefracs[k] * m_mw[k];
        }
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
                d[k] = (sumxw - m_molefracs[k] * m_mw[k])/(p * mmw * sum2);
            }
        }
    }
}

void GasTransport::getMixDiffCoeffsMole(doublereal* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

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
                d[k] = (1 - m_molefracs[k]) / (p * sum2);
            }
        }
    }
}

void GasTransport::getMixDiffCoeffsMass(doublereal* const d)
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
        for (size_t k=0; k<m_nsp; k++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (size_t i=0; i<m_nsp; i++) {
                if (i==k) {
                    continue;
                }
                sum1 += m_molefracs[i] / m_bdiff(k,i);
                sum2 += m_molefracs[i] * m_mw[i] / m_bdiff(k,i);
            }
            sum1 *= p;
            sum2 *= p * m_molefracs[k] / (mmw - m_mw[k]*m_molefracs[k]);
            d[k] = 1.0 / (sum1 +  sum2);
        }
    }
}

}
