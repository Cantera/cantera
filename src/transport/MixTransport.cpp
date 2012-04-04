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

#include <iostream>
using namespace std;

/**
 * Mole fractions below MIN_X will be set to MIN_X when computing
 * transport properties.
 */
#ifndef MIN_X
#define MIN_X 1.e-20
#endif

namespace Cantera
{


//====================================================================================================================
MixTransport::MixTransport() :
    m_tmin(-1.0),
    m_tmax(100000.),
    m_mw(0),
    m_visccoeffs(0),
    m_condcoeffs(0),
    m_diffcoeffs(0),
    m_polytempvec(0),
    m_bdiff(0, 0),
    m_visc(0),
    m_sqvisc(0),
    m_cond(0),
    m_molefracs(0),
    m_phi(0,0),
    m_wratjk(0,0),
    m_wratkj1(0,0),
    m_temp(-1.0),
    m_logt(0.0),
    m_kbt(0.0),
    m_t14(0.0),
    m_t32(0.0),
    m_sqrt_kbt(0.0),
    m_sqrt_t(0.0),
    m_lambda(0.0),
    m_viscmix(0.0),
    m_spwork(0),
    m_viscmix_ok(false),
    m_viscwt_ok(false),
    m_spvisc_ok(false),
    m_bindiff_ok(false),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_mode(0),
    m_eps(0),
    m_diam(0, 0),
    m_dipoleDiag(0),
    m_alpha(0),
    m_crot(0),
    m_zrot(0),
    m_debug(false)
{
}
//====================================================================================================================
MixTransport::MixTransport(const MixTransport& right) :
    m_tmin(-1.0),
    m_tmax(100000.),
    m_mw(0),
    m_visccoeffs(0),
    m_condcoeffs(0),
    m_diffcoeffs(0),
    m_polytempvec(0),
    m_bdiff(0, 0),
    m_visc(0),
    m_sqvisc(0),
    m_cond(0),
    m_molefracs(0),
    m_phi(0,0),
    m_wratjk(0,0),
    m_wratkj1(0,0),
    m_temp(-1.0),
    m_logt(0.0),
    m_kbt(0.0),
    m_t14(0.0),
    m_t32(0.0),
    m_sqrt_kbt(0.0),
    m_sqrt_t(0.0),
    m_lambda(0.0),
    m_viscmix(0.0),
    m_spwork(0),
    m_viscmix_ok(false),
    m_viscwt_ok(false),
    m_spvisc_ok(false),
    m_bindiff_ok(false),
    m_spcond_ok(false),
    m_condmix_ok(false),
    m_mode(0),
    m_eps(0),
    m_diam(0, 0),
    m_dipoleDiag(0),
    m_alpha(0),
    m_crot(0),
    m_zrot(0),
    m_debug(false)
{
    *this = right;
}
//====================================================================================================================
// Assignment operator
/*
 *  This is NOT a virtual function.
 *
 * @param right    Reference to %LiquidTransport object to be copied
 *                 into the current one.
 */
MixTransport&  MixTransport::operator=(const MixTransport& right)
{
    if (&right == this) {
        return *this;
    }
    Transport::operator=(right);

    m_tmin = right.m_tmin;
    m_tmax = right.m_tmax;
    m_mw =right.m_mw;
    m_visccoeffs = right.m_visccoeffs;
    m_condcoeffs = right.m_condcoeffs;
    m_diffcoeffs = right.m_diffcoeffs;
    m_polytempvec = right.m_polytempvec;
    m_bdiff = right.m_bdiff;
    m_visc = right.m_visc;
    m_sqvisc = right.m_sqvisc;
    m_cond = right.m_cond;
    m_molefracs = right.m_molefracs;
    m_phi = right.m_phi;
    m_wratjk = right.m_wratjk;
    m_wratkj1 = right.m_wratkj1;
    m_temp = right.m_temp;
    m_logt = right.m_logt;
    m_kbt = right.m_kbt;
    m_t14 = right.m_t14;
    m_t32 = right.m_t32;
    m_sqrt_kbt = right.m_sqrt_kbt;
    m_sqrt_t = right.m_sqrt_t;
    m_lambda = right.m_lambda;
    m_viscmix = right.m_viscmix;
    m_spwork = right.m_spwork;
    m_viscmix_ok = right.m_viscmix_ok;
    m_viscwt_ok = right.m_viscwt_ok;
    m_spvisc_ok = right.m_spvisc_ok;
    m_bindiff_ok = right.m_bindiff_ok;
    m_spcond_ok = right.m_spcond_ok;
    m_condmix_ok = right.m_condmix_ok;
    m_mode = right.m_mode;
    m_eps = right.m_eps;
    m_diam = right.m_diam;
    m_dipoleDiag = right.m_dipoleDiag;
    m_alpha = right.m_alpha;
    m_crot = right.m_crot;
    m_zrot = right.m_zrot;
    m_debug = right.m_debug;

    return *this;
}
//====================================================================================================================
// Duplication routine for objects which inherit from %Transport
/*
 *  This virtual routine can be used to duplicate %Transport objects
 *  inherited from %Transport even if the application only has
 *  a pointer to %Transport to work with.
 *
 *  These routines are basically wrappers around the derived copy
 *  constructor.
 */
Transport* MixTransport::duplMyselfAsTransport() const
{
    MixTransport* tr = new MixTransport(*this);
    return (dynamic_cast<Transport*>(tr));
}
//====================================================================================================================
MixTransport::~MixTransport()
{
}
//====================================================================================================================
bool MixTransport::initGas(GasTransportParams& tr)
{

    // constant substance attributes
    m_thermo = tr.thermo;
    m_nsp   = m_thermo->nSpecies();
    m_tmin  = m_thermo->minTemp();
    m_tmax  = m_thermo->maxTemp();

    // make a local copy of the molecular weights
    m_mw.resize(m_nsp);
    copy(m_thermo->molecularWeights().begin(),
         m_thermo->molecularWeights().end(), m_mw.begin());

    // copy polynomials and parameters into local storage
    m_visccoeffs = tr.visccoeffs;
    m_condcoeffs = tr.condcoeffs;
    m_diffcoeffs = tr.diffcoeffs;

    m_zrot       = tr.zrot;
    m_crot       = tr.crot;
    m_mode       = tr.mode_;
    m_diam       = tr.diam;
    m_eps        = tr.eps;
    m_alpha      = tr.alpha;
    m_dipoleDiag.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        m_dipoleDiag[i] = tr.dipole(i,i);
    }

    m_phi.resize(m_nsp, m_nsp, 0.0);
    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    size_t j, k;
    for (j = 0; j < m_nsp; j++)
        for (k = j; k < m_nsp; k++) {
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

    // set flags all false
    m_viscmix_ok = false;
    m_viscwt_ok = false;
    m_spvisc_ok = false;
    m_spcond_ok = false;
    m_condmix_ok = false;

    return true;
}
//====================================================================================================================
// Viscosity of the mixture
/*
 *
 * The viscosity is computed using the Wilke mixture rule.
 *
 *    \f[
 *        \mu = \sum_k \frac{\mu_k X_k}{\sum_j \Phi_{k,j} X_j}.
 *    \f]
 *
 *     Here \f$ \mu_k \f$ is the viscosity of pure species \e k, and
 *
 *    \f[
 *        \Phi_{k,j} = \frac{\left[1
 *                     + \sqrt{\left(\frac{\mu_k}{\mu_j}\sqrt{\frac{M_j}{M_k}}\right)}\right]^2}
 *                     {\sqrt{8}\sqrt{1 + M_k/M_j}}
 *    \f]
 *
 * @see updateViscosity_T();
 */
doublereal MixTransport::viscosity()
{
    update_T();
    update_C();

    if (m_viscmix_ok) {
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
//====================================================================================================================
// Returns the matrix of binary diffusion coefficients.
/*
 *
 *        d[ld*j + i] = rp * m_bdiff(i,j);
 *
 *  units of m**2 / s
 *
 * @param ld   offset of rows in the storage
 * @param d    output vector of diffusion coefficients
 */
void MixTransport::getBinaryDiffCoeffs(const size_t ld, doublereal* const d)
{
    update_T();
    // if necessary, evaluate the binary diffusion coefficients from the polynomial fits
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    if (ld < m_nsp) {
        throw CanteraError(" MixTransport::getBinaryDiffCoeffs()", "ld is too small");
    }
    doublereal rp = 1.0/pressure_ig();
    for (size_t i = 0; i < m_nsp; i++)
        for (size_t j = 0; j < m_nsp; j++) {
            d[ld*j + i] = rp * m_bdiff(i,j);
        }
}
//===================================================================================================================
void MixTransport::getMobilities(doublereal* const mobil)
{
    getMixDiffCoeffs(DATA_PTR(m_spwork));
    doublereal c1 = ElectronCharge / (Boltzmann * m_temp);
    for (size_t k = 0; k < m_nsp; k++) {
        mobil[k] = c1 * m_spwork[k];
    }
}
//===================================================================================================================
// Returns the mixture thermal conductivity (W/m /K)
/*
 * The thermal conductivity is computed from the following mixture rule:
 *   \f[
 *          \lambda = 0.5 \left( \sum_k X_k \lambda_k  + \frac{1}{\sum_k X_k/\lambda_k} \right)
 *   \f]
 *
 *  It's used to compute the flux of energy due to a thermal gradient
 *
 *   \f[
 *          j_T =  - \lambda  \nabla T
 *   \f]
 *
 *  The flux of energy has units of energy (kg m2 /s2) per second per area.
 *
 *  The units of lambda are W / m K which is equivalent to kg m / s^3 K.
 *
 * @return Returns the mixture thermal conductivity, with units of W/m/K
 */
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
//===================================================================================================================
// Return the thermal diffusion coefficients
/*
 * For this approximation, these are all zero.
 *
 *  Eqns. (12.168) shows how they are used in an expression for the species flux.
 *
 * @param dt  Vector of thermal diffusion coefficients. Units = kg/m/s
 */
void MixTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    for (size_t k = 0; k < m_nsp; k++) {
        dt[k] = 0.0;
    }
}
//===================================================================================================================
// Get the species diffusive mass fluxes wrt to the mass averaged velocity,
// given the gradients in mole fraction and temperature
/*
 *  Units for the returned fluxes are kg m-2 s-1.
 *
 *
 * The diffusive mass flux of species \e k is computed from
 * \f[
 *          \vec{j}_k = -n M_k D_k \nabla X_k.
 * \f]
 *
 *  @param ndim      Number of dimensions in the flux expressions
 *  @param grad_T    Gradient of the temperature
 *                    (length = ndim)
 * @param ldx        Leading dimension of the grad_X array
 *                   (usually equal to m_nsp but not always)
 * @param grad_X     Gradients of the mole fraction
 *                   Flat vector with the m_nsp in the inner loop.
 *                   length = ldx * ndim
 * @param ldf  Leading dimension of the fluxes array
 *              (usually equal to m_nsp but not always)
 * @param fluxes  Output of the diffusive mass fluxes
 *             Flat vector with the m_nsp in the inner loop.
 *             length = ldx * ndim
 */
void MixTransport::getSpeciesFluxes(size_t ndim,
                                    const doublereal* const grad_T, int ldx, const doublereal* const grad_X,
                                    int ldf, doublereal* const fluxes)
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
//===========================================================================================================
// Mixture-averaged diffusion coefficients [m^2/s].
/*
 * Returns the mixture averaged diffusion coefficients for a gas.
 * Note, for the single species case or the pure fluid case the routine returns the self-diffusion coefficient.
 * This is need to avoid a Nan result in the formula
 * below.
 *
 *  @param d  Output Vector of diffusion coefficients for each species (m^2/s)
 *            length m_nsp
 */
void MixTransport::getMixDiffCoeffs(doublereal* const d)
{
    update_T();
    update_C();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    doublereal mmw = m_thermo->meanMolecularWeight();
    doublereal sumxw = 0.0, sum2;
    doublereal p = pressure_ig();
    if (m_nsp == 1) {
        d[0] = m_bdiff(0,0) / p;
    } else {
        for (size_t k = 0; k < m_nsp; k++) {
            sumxw += m_molefracs[k] * m_mw[k];
        }
        for (size_t k = 0; k < m_nsp; k++) {
            sum2 = 0.0;
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
//===========================================================================================================
/*
 *  @internal This is called whenever a transport property is
 *  requested from ThermoSubstance if the temperature has changed
 *  since the last call to update_T.
 */
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

    // temperature has changed, so polynomial fits will need to be
    // redone.
    m_viscmix_ok = false;
    m_spvisc_ok = false;
    m_viscwt_ok = false;
    m_spcond_ok = false;
    m_bindiff_ok = false;
    m_condmix_ok = false;
}
//====================================================================================================================
/*
 *  @internal This is called the first time any transport property
 *  is requested from Mixture after the concentrations
 *  have changed.
 */
void MixTransport::update_C()
{
    // signal that concentration-dependent quantities will need to
    // be recomputed before use, and update the local mole
    // fractions.

    m_viscmix_ok = false;
    m_condmix_ok = false;

    m_thermo->getMoleFractions(DATA_PTR(m_molefracs));

    // add an offset to avoid a pure species condition
    for (size_t k = 0; k < m_nsp; k++) {
        m_molefracs[k] = std::max(MIN_X, m_molefracs[k]);
    }
}
//====================================================================================================================
/*
 * Update the temperature-dependent parts of the mixture-averaged
 * thermal conductivity.
 */
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
//====================================================================================================================
/*
 * Update the binary diffusion coefficients. These are evaluated
 * from the polynomial fits at unit pressure (1 Pa).
 */
void MixTransport::updateDiff_T()
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
}
//====================================================================================================================
/*
 * Update the pure-species viscosities.
 */
void MixTransport::updateSpeciesViscosities()
{
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
//====================================================================================================================
// Update the temperature-dependent viscosity terms.
/*
 * Updates the array of pure species viscosities, and the weighting functions in the viscosity mixture rule.
 * The flag m_visc_ok is set to true.
 *
 * The formula for the weighting function is from Poling and Prausnitz.
 * See Eq. (9-5.14) of  Poling, Prausnitz, and O'Connell. The equation for the weighting function
 * \f$ \phi_{ij} \f$ is reproduced below.
 *
 *  \f[
 *        \phi_{ij} = \frac{ \left[ 1 + \left( \mu_i / \mu_j \right)^{1/2} \left( M_j / M_i \right)^{1/4} \right]^2 }
 *                    {\left[ 8 \left( 1 + M_i / M_j \right) \right]^{1/2}}
 *  \f]
 */
void MixTransport::updateViscosity_T()
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
            m_phi(k,j) = factor1*factor1 / (SqrtEight * m_wratkj1(j,k));
            m_phi(j,k) = m_phi(k,j)/(vratiokj * wratiojk);
        }
    }
    m_viscwt_ok = true;
}
//====================================================================================================================
/*
 * This function returns a Transport data object for a given species.
 *
 */
struct GasTransportData MixTransport::getGasTransportData(int kSpecies) const {

    struct GasTransportData td;
    td.speciesName = m_thermo->speciesName(kSpecies);

    td.geometry = 2;
    if (m_crot[kSpecies] == 0.0) {
        td.geometry = 0;
    } else if (m_crot[kSpecies] == 1.0) {
        td.geometry = 1;
    }
    td.wellDepth = m_eps[kSpecies] / Boltzmann;
    td.dipoleMoment = m_dipoleDiag[kSpecies] * 1.0E25 / SqrtTen;
    td.diameter = m_diam(kSpecies, kSpecies) * 1.0E10;
    td.polarizability = m_alpha[kSpecies] * 1.0E30;
    td.rotRelaxNumber = m_zrot[kSpecies];

    return td;
}
//====================================================================================================================
}

