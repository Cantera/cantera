/**
 *  @file AqueousTransport.cpp
 *  Transport properties for aqueous systems
 */

#include "cantera/thermo/ThermoPhase.h"
#include "cantera/transport/AqueousTransport.h"

#include "cantera/base/utilities.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/transport/TransportFactory.h"

#include "cantera/numerics/ctlapack.h"

#include "cantera/base/stringUtils.h"

#include <cstdio>

using namespace std;

namespace Cantera
{

AqueousTransport::AqueousTransport() :
    m_iStateMF(-1),
    m_temp(-1.0),
    m_logt(0.0),
    m_sqrt_t(-1.0),
    m_t14(-1.0),
    m_t32(-1.0),
    m_sqrt_kbt(-1.0),
    m_press(-1.0),
    m_lambda(-1.0),
    m_viscmix(-1.0),
    m_viscmix_ok(false),
    m_viscwt_ok(false),
    m_spvisc_ok(false),
    m_diffmix_ok(false),
    m_bindiff_ok(false),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_mode(-1000),
    m_debug(false),
    m_nDim(1)
{
}

bool AqueousTransport::initLiquid(LiquidTransportParams& tr)
{
    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(),
         m_thermo->molecularWeights().end(), m_mw.begin());

    // copy polynomials and parameters into local storage
    //m_visccoeffs = tr.visccoeffs;
    //m_condcoeffs = tr.condcoeffs;
    //m_diffcoeffs = tr.diffcoeffs;
    cout << "In AqueousTransport::initLiquid we need to replace" << endl
         << "LiquidTransportParams polynomial coefficients with" << endl
         <<  "those in LiquidTransportData as in SimpleTransport." << endl;

    m_mode       = tr.mode_;

    m_phi.resize(m_nsp, m_nsp, 0.0);


    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    for (size_t j = 0; j < m_nsp; j++)
        for (size_t k = j; k < m_nsp; k++) {
            m_wratjk(j,k) = sqrt(m_mw[j]/m_mw[k]);
            m_wratjk(k,j) = sqrt(m_wratjk(j,k));
            m_wratkj1(j,k) = sqrt(1.0 + m_mw[k]/m_mw[j]);
        }

    m_polytempvec.resize(5);
    m_visc.resize(m_nsp);
    m_sqvisc.resize(m_nsp);
    m_cond.resize(m_nsp);
    m_bdiff.resize(m_nsp, m_nsp);

    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);

    // resize the internal gradient variables
    m_Grad_X.resize(m_nDim * m_nsp, 0.0);
    m_Grad_T.resize(m_nDim, 0.0);
    m_Grad_V.resize(m_nDim, 0.0);
    m_Grad_mu.resize(m_nDim * m_nsp, 0.0);


    // set all flags to false
    m_viscmix_ok = false;
    m_viscwt_ok  = false;
    m_spvisc_ok  = false;
    m_spcond_ok  = false;
    m_condmix_ok = false;
    m_spcond_ok  = false;
    m_diffmix_ok = false;

    return true;
}

doublereal AqueousTransport::viscosity()
{

    update_T();
    update_C();

    if (m_viscmix_ok) {
        return m_viscmix;
    }

    // update m_visc[] and m_phi[] if necessary
    if (!m_viscwt_ok) {
        updateViscosity_T();
    }

    multiply(m_phi, DATA_PTR(m_molefracs), DATA_PTR(m_spwork));

    m_viscmix = 0.0;
    for (size_t k = 0; k < m_nsp; k++) {
        m_viscmix += m_molefracs[k] * m_visc[k]/m_spwork[k]; //denom;
    }
    return m_viscmix;
}

void AqueousTransport::getSpeciesViscosities(doublereal* const visc)
{
    updateViscosity_T();
    copy(m_visc.begin(), m_visc.end(), visc);
}

void AqueousTransport::getBinaryDiffCoeffs(const size_t ld, doublereal* const d)
{
    update_T();

    // if necessary, evaluate the binary diffusion coefficients
    // from the polynomial fits
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    doublereal pres = m_thermo->pressure();

    doublereal rp = 1.0/pres;
    for (size_t i = 0; i < m_nsp; i++)
        for (size_t j = 0; j < m_nsp; j++) {
            d[ld*j + i] = rp * m_bdiff(i,j);
        }
}

void AqueousTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

void AqueousTransport::getFluidMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = 1.0 / (GasConstant * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}

void AqueousTransport::set_Grad_V(const doublereal* const grad_V)
{
    for (size_t a = 0; a < m_nDim; a++) {
        m_Grad_V[a] = grad_V[a];
    }
}

void AqueousTransport::set_Grad_T(const doublereal* const grad_T)
{
    for (size_t a = 0; a < m_nDim; a++) {
        m_Grad_T[a] = grad_T[a];
    }
}

void AqueousTransport::set_Grad_X(const doublereal* const grad_X)
{
    size_t itop = m_nDim * m_nsp;
    for (size_t i = 0; i < itop; i++) {
        m_Grad_X[i] = grad_X[i];
    }
}

doublereal AqueousTransport::thermalConductivity()
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
    }
    return m_lambda;
}

void AqueousTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}

void AqueousTransport::getSpeciesFluxes(size_t ndim, const doublereal* const grad_T,
                                        size_t ldx, const doublereal* const grad_X,
                                        size_t ldf, doublereal* const fluxes)
{
    set_Grad_T(grad_T);
    set_Grad_X(grad_X);
    getSpeciesFluxesExt(ldf, fluxes);
}

void AqueousTransport::getSpeciesFluxesExt(size_t ldf, doublereal* const fluxes)
{
    update_T();
    update_C();

    getMixDiffCoeffs(DATA_PTR(m_spwork));

    const vector_fp& mw = m_thermo->molecularWeights();
    const doublereal* y  = m_thermo->massFractions();
    doublereal rhon = m_thermo->molarDensity();
    // Unroll wrt ndim
    vector_fp sum(m_nDim,0.0);
    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] = -rhon * mw[k] * m_spwork[k] * m_Grad_X[n*m_nsp + k];
            sum[n] += fluxes[n*ldf + k];
        }
    }
    // add correction flux to enforce sum to zero
    for (size_t n = 0; n < m_nDim; n++) {
        for (size_t k = 0; k < m_nsp; k++) {
            fluxes[n*ldf + k] -= y[k]*sum[n];
        }
    }
}

void AqueousTransport::getMixDiffCoeffs(doublereal* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    size_t k, j;
    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal sumxw = 0.0, sum2;
    doublereal p = m_press;
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (k = 0; k < m_nsp; k++) {
            sumxw += m_molefracs[k] * m_mw[k];
        }
        for (k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
            for (j = 0; j < m_nsp; j++) {
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

void AqueousTransport::update_T()
{
    doublereal t = m_thermo->temperature();
    if (t == m_temp) {
        return;
    }
    if (t < 0.0) {
        throw CanteraError("AqueousTransport::update_T",
                           "negative temperature "+fp2str(t));
    }

    // Compute various functions of temperature
    m_temp = t;
    m_logt = log(m_temp);
    m_kbt = Boltzmann * m_temp;
    m_sqrt_t = sqrt(m_temp);
    m_t14 = sqrt(m_sqrt_t);
    m_t32 = m_temp * m_sqrt_t;
    m_sqrt_kbt = sqrt(Boltzmann*m_temp);

    // compute powers of log(T)
    m_polytempvec[0] = 1.0;
    m_polytempvec[1] = m_logt;
    m_polytempvec[2] = m_logt*m_logt;
    m_polytempvec[3] = m_logt*m_logt*m_logt;
    m_polytempvec[4] = m_logt*m_logt*m_logt*m_logt;

    // temperature has changed, so polynomial temperature
    // interpolations will need to be reevaluated.
    // Set all of these flags to false
    m_viscmix_ok = false;
    m_spvisc_ok  = false;
    m_viscwt_ok  = false;
    m_spcond_ok  = false;
    m_diffmix_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;

    // For now, for a concentration redo also
    m_iStateMF   = -1;
}

void AqueousTransport::update_C()
{

    doublereal pres = m_thermo->pressure();
    // Check for changes in the mole fraction vector.
    //int iStateNew = m_thermo->getIStateMF();
    //if (iStateNew == m_iStateMF) {
    //  if (pres == m_press) {
    //    return;
    //     }
    // } else {
    //  m_iStateMF = iStateNew;
    //}
    m_press = pres;

    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_viscmix_ok = false;
    m_diffmix_ok = false;
    m_condmix_ok = false;

    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // add an offset to avoid a pure species condition or
    // negative mole fractions. *Tiny* is 1.0E-20, a value
    // which is below the additive machine precision of mole fractions.
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(Tiny, m_molefracs[k]);
    }
}

void AqueousTransport::updateCond_T()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = exp(dot4(m_polytempvec, m_condcoeffs[k]));
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            m_cond[k] = m_sqrt_t*dot5(m_polytempvec, m_condcoeffs[k]);
        }
    }
    m_spcond_ok = true;
    m_condmix_ok = false;
}

void AqueousTransport::updateDiff_T()
{
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
    m_diffmix_ok = false;
}

void AqueousTransport::updateSpeciesViscosities()
{
    if (m_mode == CK_Mode) {
        for (size_t k = 0; k < m_nsp; k++) {
            m_visc[k] = exp(dot4(m_polytempvec, m_visccoeffs[k]));
            m_sqvisc[k] = sqrt(m_visc[k]);
        }
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            // the polynomial fit is done for sqrt(visc/sqrt(T))
            m_sqvisc[k] = m_t14*dot5(m_polytempvec, m_visccoeffs[k]);
            m_visc[k] = (m_sqvisc[k]*m_sqvisc[k]);
        }
    }
    m_spvisc_ok = true;
}

void AqueousTransport::updateViscosity_T()
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

            // Note that m_wratjk(k,j) holds the square root of
            // m_wratjk(j,k)!
            factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
            m_phi(k,j) = factor1*factor1 /
                         (SqrtEight * m_wratkj1(j,k));
            m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
        }
    }
    m_viscwt_ok = true;
}

LiquidTransportData AqueousTransport::getLiquidTransportData(int kSpecies)
{
    LiquidTransportData td;
    td.speciesName = m_thermo->speciesName(kSpecies);
    return td;
}

void AqueousTransport::stefan_maxwell_solve()
{
    size_t VIM = 2;
    m_B.resize(m_nsp, VIM);
    // grab a local copy of the molecular weights
    const vector_fp& M =  m_thermo->molecularWeights();


    // get the mean molecular weight of the mixture
    //double M_mix = m_thermo->meanMolecularWeight();


    // get the concentration of the mixture
    //double rho = m_thermo->density();
    //double c = rho/M_mix;


    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    double T = m_thermo->temperature();


    /* electrochemical potential gradient */
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t a = 0; a < VIM; a++) {
            m_Grad_mu[a*m_nsp + i] = m_chargeSpecies[i] * Faraday * m_Grad_V[a]
                                     + (GasConstant*T/m_molefracs[i]) * m_Grad_X[a*m_nsp+i];
        }
    }

    /*
     * Just for Note, m_A(i,j) refers to the ith row and jth column.
     * They are still fortran ordered, so that i varies fastest.
     */
    switch (VIM) {
    case 1:  /* 1-D approximation */
        m_B(0,0) = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            m_A(0,j) = 1.0;
        }
        for (size_t i = 1; i < m_nsp; i++) {
            m_B(i,0) = m_concentrations[i] * m_Grad_mu[i] / (GasConstant * T);
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    m_A(i,j)  = m_molefracs[i] / (M[j] * m_DiffCoeff_StefMax(i,j));
                    m_A(i,i) -= m_molefracs[j] / (M[i] * m_DiffCoeff_StefMax(i,j));
                } else if (j == i)  {
                    m_A(i,i) = 0.0;
                }
            }
        }

        //! invert and solve the system  Ax = b. Answer is in m_B
        solve(m_A, m_B.ptrColumn(0));

        m_flux = m_B;


        break;
    case 2:  /* 2-D approximation */
        m_B(0,0) = 0.0;
        m_B(0,1) = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            m_A(0,j) = 1.0;
        }
        for (size_t i = 1; i < m_nsp; i++) {
            m_B(i,0) = m_concentrations[i] * m_Grad_mu[i] / (GasConstant * T);
            m_B(i,1) = m_concentrations[i] * m_Grad_mu[m_nsp + i] / (GasConstant * T);
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    m_A(i,j)  = m_molefracs[i] / (M[j] * m_DiffCoeff_StefMax(i,j));
                    m_A(i,i) -= m_molefracs[j] / (M[i] * m_DiffCoeff_StefMax(i,j));
                } else if (j == i)  {
                    m_A(i,i) = 0.0;
                }
            }
        }

        //! invert and solve the system  Ax = b. Answer is in m_B
        //solve(m_A, m_B);

        m_flux = m_B;


        break;

    case 3:  /* 3-D approximation */
        m_B(0,0) = 0.0;
        m_B(0,1) = 0.0;
        m_B(0,2) = 0.0;
        for (size_t j = 0; j < m_nsp; j++) {
            m_A(0,j) = 1.0;
        }
        for (size_t i = 1; i < m_nsp; i++) {
            m_B(i,0) = m_concentrations[i] * m_Grad_mu[i] / (GasConstant * T);
            m_B(i,1) = m_concentrations[i] * m_Grad_mu[m_nsp + i] / (GasConstant * T);
            m_B(i,2) = m_concentrations[i] * m_Grad_mu[2*m_nsp + i] / (GasConstant * T);
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    m_A(i,j)  = m_molefracs[i] / (M[j] * m_DiffCoeff_StefMax(i,j));
                    m_A(i,i) -= m_molefracs[j] / (M[i] * m_DiffCoeff_StefMax(i,j));
                } else if (j == i)  {
                    m_A(i,i) = 0.0;
                }
            }
        }

        //! invert and solve the system  Ax = b. Answer is in m_B
        //solve(m_A, m_B);

        m_flux = m_B;


        break;
    default:
        printf("unimplemented\n");
        throw CanteraError("routine", "not done");
        break;
    }
}

}
