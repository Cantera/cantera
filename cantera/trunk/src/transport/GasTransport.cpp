//! @file GasTransport.cpp
#include "cantera/transport/GasTransport.h"
#include "cantera/transport/TransportParams.h"
#include "MMCollisionInt.h"
#include "cantera/base/stringUtils.h"
#include "cantera/numerics/polyfit.h"
#include "cantera/transport/TransportData.h"

namespace Cantera
{

//! polynomial degree used for fitting collision integrals
//! except in CK mode, where the degree is 6.
#define COLL_INT_POLY_DEGREE 8

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
    m_bdiff(0, 0),
    m_verbose(false)
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
    m_bdiff(0, 0),
    m_verbose(false)
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
    m_verbose = right.m_verbose;

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

void GasTransport::init(thermo_t* thermo, int mode, int log_level)
{
    GasTransportParams trParam;
    if (log_level == 0) {
        m_verbose = 0;
    }
    // set up Monchick and Mason collision integrals
    setupMM(thermo, mode, log_level, trParam);
    // do model-specific initialization
    initGas(trParam);
}

void GasTransport::setupMM(thermo_t* thermo, int mode, int log_level,
                           GasTransportParams& tr)
{
    // constant mixture attributes
    tr.thermo = thermo;
    tr.nsp_ = tr.thermo->nSpecies();
    size_t nsp = tr.nsp_;

    tr.tmin = thermo->minTemp();
    tr.tmax = thermo->maxTemp();
    tr.mw.resize(nsp);
    tr.log_level = log_level;

    copy(tr.thermo->molecularWeights().begin(),
         tr.thermo->molecularWeights().end(), tr.mw.begin());

    tr.mode_ = mode;
    tr.epsilon.resize(nsp, nsp, 0.0);
    tr.delta.resize(nsp, nsp, 0.0);
    tr.reducedMass.resize(nsp, nsp, 0.0);
    tr.dipole.resize(nsp, nsp, 0.0);
    tr.diam.resize(nsp, nsp, 0.0);
    tr.crot.resize(nsp);
    tr.zrot.resize(nsp);
    tr.polar.resize(nsp, false);
    tr.alpha.resize(nsp, 0.0);
    tr.poly.resize(nsp);
    tr.sigma.resize(nsp);
    tr.eps.resize(nsp);
    tr.w_ac.resize(nsp);

    getTransportData(*thermo, tr);

    for (size_t i = 0; i < nsp; i++) {
        tr.poly[i].resize(nsp);
    }

    double tstar_min = 1.e8, tstar_max = 0.0;
    double f_eps, f_sigma;

    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = i; j < nsp; j++) {
            // the reduced mass
            tr.reducedMass(i,j) =  tr.mw[i] * tr.mw[j] / (Avogadro * (tr.mw[i] + tr.mw[j]));

            // hard-sphere diameter for (i,j) collisions
            tr.diam(i,j) = 0.5*(tr.sigma[i] + tr.sigma[j]);

            // the effective well depth for (i,j) collisions
            tr.epsilon(i,j) = sqrt(tr.eps[i]*tr.eps[j]);

            //  The polynomial fits of collision integrals vs. T*
            //  will be done for the T* from tstar_min to tstar_max
            tstar_min = std::min(tstar_min, Boltzmann * tr.tmin/tr.epsilon(i,j));
            tstar_max = std::max(tstar_max, Boltzmann * tr.tmax/tr.epsilon(i,j));

            // the effective dipole moment for (i,j) collisions
            tr.dipole(i,j) = sqrt(tr.dipole(i,i)*tr.dipole(j,j));

            // reduced dipole moment delta* (nondimensional)
            double d = tr.diam(i,j);
            tr.delta(i,j) =  0.5 * tr.dipole(i,j)*tr.dipole(i,j)
                             / (4 * Pi * epsilon_0 * tr.epsilon(i,j) * d * d * d);

            makePolarCorrections(i, j, tr, f_eps, f_sigma);
            tr.diam(i,j) *= f_sigma;
            tr.epsilon(i,j) *= f_eps;

            // properties are symmetric
            tr.reducedMass(j,i) = tr.reducedMass(i,j);
            tr.diam(j,i) = tr.diam(i,j);
            tr.epsilon(j,i) = tr.epsilon(i,j);
            tr.dipole(j,i)  = tr.dipole(i,j);
            tr.delta(j,i)   = tr.delta(i,j);
        }
    }

    // Chemkin fits the entire T* range in the Monchick and Mason tables,
    // so modify tstar_min and tstar_max if in Chemkin compatibility mode
    if (mode == CK_Mode) {
        tstar_min = 0.101;
        tstar_max = 99.9;
    }

    // initialize the collision integral calculator for the desired T* range
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** collision_integrals ***\n");
    }
    MMCollisionInt integrals;
    integrals.init(tstar_min, tstar_max, log_level);
    fitCollisionIntegrals(tr, integrals);
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** end of collision_integrals ***\n");
    }
    // make polynomial fits
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** property fits ***\n");
    }
    fitProperties(tr, integrals);
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("*** end of property fits ***\n");
    }
}

void GasTransport::getTransportData(const ThermoPhase& thermo,
                                    GasTransportParams& tr)
{
    for (size_t k = 0; k < thermo.nSpecies(); k++) {
        const Species& s = thermo.species(thermo.speciesName(k));
        const GasTransportData& sptran =
            dynamic_cast<GasTransportData&>(*s.transport.get());
        if (sptran.geometry == "atom") {
            tr.crot[k] = 0.0;
        } else if (sptran.geometry == "linear") {
            tr.crot[k] = 1.0;
        } else if (sptran.geometry == "nonlinear") {
            tr.crot[k] = 1.5;
        }

        tr.sigma[k] = sptran.diameter;
        tr.eps[k] = sptran.well_depth;
        tr.dipole(k,k) = sptran.dipole;
        tr.polar[k] = (sptran.dipole > 0);
        tr.alpha[k] = sptran.polarizability;
        tr.zrot[k] = sptran.rotational_relaxation;
        tr.w_ac[k] = sptran.acentric_factor;
    }
}

void GasTransport::makePolarCorrections(size_t i, size_t j,
        const GasTransportParams& tr, doublereal& f_eps, doublereal& f_sigma)
{
    // no correction if both are nonpolar, or both are polar
    if (tr.polar[i] == tr.polar[j]) {
        f_eps = 1.0;
        f_sigma = 1.0;
        return;
    }

    // corrections to the effective diameter and well depth
    // if one is polar and one is non-polar

    size_t kp = (tr.polar[i] ? i : j);     // the polar one
    size_t knp = (i == kp ? j : i);        // the nonpolar one

    double d3np, d3p, alpha_star, mu_p_star, xi;
    d3np = pow(tr.sigma[knp],3);
    d3p  = pow(tr.sigma[kp],3);
    alpha_star = tr.alpha[knp]/d3np;
    mu_p_star  = tr.dipole(kp,kp)/sqrt(4 * Pi * epsilon_0 * d3p * tr.eps[kp]);
    xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star *
         sqrt(tr.eps[kp]/tr.eps[knp]);
    f_sigma = pow(xi, -1.0/6.0);
    f_eps = xi*xi;
}

void GasTransport::fitCollisionIntegrals(GasTransportParams& tr,
                                         MMCollisionInt& integrals)
{
    vector_fp::iterator dptr;
    double dstar;
    size_t nsp = tr.nsp_;
    int mode = tr.mode_;

    // Chemkin fits to sixth order polynomials
    int degree = (mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("tstar_fits\n"
                 "fits to A*, B*, and C* vs. log(T*).\n"
                 "These are done only for the required dstar(j,k) values.\n\n");
        if (tr.log_level < 3) {
            writelog("*** polynomial coefficients not printed (log_level < 3) ***\n");
        }
    }
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = i; j < nsp; j++)  {
            // Chemkin fits only delta* = 0
            if (mode != CK_Mode) {
                dstar = tr.delta(i,j);
            } else {
                dstar = 0.0;
            }

            // if a fit has already been generated for delta* = tr.delta(i,j),
            // then use it. Otherwise, make a new fit, and add tr.delta(i,j) to
            // the list of delta* values for which fits have been done.

            // 'find' returns a pointer to end() if not found
            dptr = find(tr.fitlist.begin(), tr.fitlist.end(), dstar);
            if (dptr == tr.fitlist.end()) {
                vector_fp ca(degree+1), cb(degree+1), cc(degree+1);
                vector_fp co22(degree+1);
                integrals.fit(degree, dstar,
                              DATA_PTR(ca), DATA_PTR(cb), DATA_PTR(cc));
                integrals.fit_omega22(degree, dstar,
                                      DATA_PTR(co22));
                tr.omega22_poly.push_back(co22);
                tr.astar_poly.push_back(ca);
                tr.bstar_poly.push_back(cb);
                tr.cstar_poly.push_back(cc);
                tr.poly[i][j] = static_cast<int>(tr.astar_poly.size()) - 1;
                tr.fitlist.push_back(dstar);
            }

            // delta* found in fitlist, so just point to this polynomial
            else {
                tr.poly[i][j] = static_cast<int>((dptr - tr.fitlist.begin()));
            }
            tr.poly[j][i] = tr.poly[i][j];
        }
    }
}

void GasTransport::fitProperties(GasTransportParams& tr,
                                 MMCollisionInt& integrals)
{
    int ndeg = 0;
    // number of points to use in generating fit data
    const size_t np = 50;

    int mode = tr.mode_;
    int degree = (mode == CK_Mode ? 3 : 4);

    double dt = (tr.tmax - tr.tmin)/(np-1);
    vector_fp tlog(np), spvisc(np), spcond(np);
    vector_fp w(np), w2(np);

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        double t = tr.tmin + dt*n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector_fp c(degree + 1), c2(degree + 1);

    // fit the pure-species viscosity and thermal conductivity for each species
    if (DEBUG_MODE_ENABLED && tr.log_level < 2 && m_verbose) {
        writelog("*** polynomial coefficients not printed (log_level < 2) ***\n");
    }
    double sqrt_T, visc, err, relerr,
               mxerr = 0.0, mxrelerr = 0.0, mxerr_cond = 0.0, mxrelerr_cond = 0.0;

    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelog("Polynomial fits for viscosity:\n");
        if (mode == CK_Mode) {
            writelog("log(viscosity) fit to cubic polynomial in log(T)\n");
        } else {
            writelogf("viscosity/sqrt(T) fit to polynomial of degree "
                      "%d in log(T)", degree);
        }
    }

    double cp_R, cond, w_RT, f_int, A_factor, B_factor, c1, cv_rot, cv_int,
           f_rot, f_trans, om11, diffcoeff;

    for (size_t k = 0; k < tr.nsp_; k++) {
        for (size_t n = 0; n < np; n++) {
            double t = tr.tmin + dt*n;

            tr.thermo->setTemperature(t);
            vector_fp cp_R_all(tr.thermo->nSpecies());
            tr.thermo->getCp_R_ref(&cp_R_all[0]);
            cp_R = cp_R_all[k];

            double tstar = Boltzmann * t/ tr.eps[k];
            sqrt_T = sqrt(t);
            double om22 = integrals.omega22(tstar, tr.delta(k,k));
            om11 = integrals.omega11(tstar, tr.delta(k,k));

            // self-diffusion coefficient, without polar corrections
            diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/tr.reducedMass(k,k)) *
                        pow((Boltzmann * t), 1.5)/
                        (Pi * tr.sigma[k] * tr.sigma[k] * om11);

            // viscosity
            visc = FiveSixteenths
                   * sqrt(Pi * tr.mw[k] * Boltzmann * t / Avogadro) /
                   (om22 * Pi * tr.sigma[k]*tr.sigma[k]);

            // thermal conductivity
            w_RT = tr.mw[k]/(GasConstant * t);
            f_int = w_RT * diffcoeff/visc;
            cv_rot = tr.crot[k];

            A_factor = 2.5 - f_int;
            B_factor = tr.zrot[k] + 2.0/Pi * (5.0/3.0 * cv_rot + f_int);
            c1 = 2.0/Pi * A_factor/B_factor;
            cv_int = cp_R - 2.5 - cv_rot;

            f_rot = f_int * (1.0 + c1);
            f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);

            cond = (visc/tr.mw[k])*GasConstant*(f_trans * 1.5
                                                + f_rot * cv_rot + f_int * cv_int);

            if (mode == CK_Mode) {
                spvisc[n] = log(visc);
                spcond[n] = log(cond);
                w[n] = -1.0;
                w2[n] = -1.0;
            } else {
                // the viscosity should be proportional approximately to
                // sqrt(T); therefore, visc/sqrt(T) should have only a weak
                // temperature dependence. And since the mixture rule requires
                // the square root of the pure-species viscosity, fit the square
                // root of (visc/sqrt(T)) to avoid having to compute square
                // roots in the mixture rule.
                spvisc[n] = sqrt(visc/sqrt_T);

                // the pure-species conductivity scales approximately with
                // sqrt(T). Unlike the viscosity, there is no reason here to fit
                // the square root, since a different mixture rule is used.
                spcond[n] = cond/sqrt_T;
                w[n] = 1.0/(spvisc[n]*spvisc[n]);
                w2[n] = 1.0/(spcond[n]*spcond[n]);
            }
        }
        polyfit(np, DATA_PTR(tlog), DATA_PTR(spvisc),
                DATA_PTR(w), degree, ndeg, 0.0, DATA_PTR(c));
        polyfit(np, DATA_PTR(tlog), DATA_PTR(spcond),
                DATA_PTR(w), degree, ndeg, 0.0, DATA_PTR(c2));

        // evaluate max fit errors for viscosity
        for (size_t n = 0; n < np; n++) {
            double val, fit;
            if (mode == CK_Mode) {
                val = exp(spvisc[n]);
                fit = exp(poly3(tlog[n], DATA_PTR(c)));
            } else {
                sqrt_T = exp(0.5*tlog[n]);
                val = sqrt_T * pow(spvisc[n],2);
                fit = sqrt_T * pow(poly4(tlog[n], DATA_PTR(c)),2);
            }
            err = fit - val;
            relerr = err/val;
            mxerr = std::max(mxerr, fabs(err));
            mxrelerr = std::max(mxrelerr, fabs(relerr));
        }

        // evaluate max fit errors for conductivity
        for (size_t n = 0; n < np; n++) {
            double val, fit;
            if (mode == CK_Mode) {
                val = exp(spcond[n]);
                fit = exp(poly3(tlog[n], DATA_PTR(c2)));
            } else {
                sqrt_T = exp(0.5*tlog[n]);
                val = sqrt_T * spcond[n];
                fit = sqrt_T * poly4(tlog[n], DATA_PTR(c2));
            }
            err = fit - val;
            relerr = err/val;
            mxerr_cond = std::max(mxerr_cond, fabs(err));
            mxrelerr_cond = std::max(mxrelerr_cond, fabs(relerr));
        }
        tr.visccoeffs.push_back(c);
        tr.condcoeffs.push_back(c2);

        if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
            writelog(tr.thermo->speciesName(k) + ": [" + vec2str(c) + "]\n");
        }
    }
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelogf("Maximum viscosity absolute error:  %12.6g\n", mxerr);
        writelogf("Maximum viscosity relative error:  %12.6g\n", mxrelerr);

        writelog("\nPolynomial fits for conductivity:\n");
        if (mode == CK_Mode)
            writelog("log(conductivity) fit to cubic polynomial in log(T)");
        else {
            writelogf("conductivity/sqrt(T) fit to "
                      "polynomial of degree %d in log(T)", degree);
        }
        if (tr.log_level >= 2)
            for (size_t k = 0; k < tr.nsp_; k++) {
                writelog(tr.thermo->speciesName(k) + ": [" +
                         vec2str(tr.condcoeffs[k]) + "]\n");
            }
        writelogf("Maximum conductivity absolute error:  %12.6g\n", mxerr_cond);
        writelogf("Maximum conductivity relative error:  %12.6g\n", mxrelerr_cond);

        // fit the binary diffusion coefficients for each species pair
        writelogf("\nbinary diffusion coefficients:\n");
        if (mode == CK_Mode)
            writelog("log(D) fit to cubic polynomial in log(T)");
        else {
            writelogf("D/T**(3/2) fit to polynomial of degree %d in log(T)",degree);
        }
    }

    mxerr = 0.0, mxrelerr = 0.0;
    vector_fp diff(np + 1);
    double eps, sigma;
    for (size_t k = 0; k < tr.nsp_; k++)  {
        for (size_t j = k; j < tr.nsp_; j++) {
            for (size_t n = 0; n < np; n++) {
                double t = tr.tmin + dt*n;
                eps = tr.epsilon(j,k);
                double tstar = Boltzmann * t/eps;
                sigma = tr.diam(j,k);
                om11 = integrals.omega11(tstar, tr.delta(j,k));

                diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/tr.reducedMass(k,j)) *
                            pow(Boltzmann * t, 1.5) /
                            (Pi * sigma * sigma * om11);

                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                double fkj, fjk;
                getBinDiffCorrection(t, tr, integrals, k, j, 1.0, 1.0, fkj, fjk);

                if (mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;
                } else {
                    diff[n] = diffcoeff/pow(t, 1.5);
                    w[n] = 1.0/(diff[n]*diff[n]);
                }
            }
            polyfit(np, DATA_PTR(tlog), DATA_PTR(diff),
                    DATA_PTR(w), degree, ndeg, 0.0, DATA_PTR(c));

            for (size_t n = 0; n < np; n++) {
                double val, fit;
                if (mode == CK_Mode) {
                    val = exp(diff[n]);
                    fit = exp(poly3(tlog[n], DATA_PTR(c)));
                } else {
                    double t = exp(tlog[n]);
                    double pre = pow(t, 1.5);
                    val = pre * diff[n];
                    fit = pre * poly4(tlog[n], DATA_PTR(c));
                }
                err = fit - val;
                relerr = err/val;
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }
            tr.diffcoeffs.push_back(c);
            if (DEBUG_MODE_ENABLED && tr.log_level >= 2 && m_verbose) {
                writelog(tr.thermo->speciesName(k) + "__" +
                         tr.thermo->speciesName(j) + ": [" + vec2str(c) + "]\n");
            }
        }
    }
    if (DEBUG_MODE_ENABLED && m_verbose) {
        writelogf("Maximum binary diffusion coefficient absolute error:"
                 "  %12.6g\n", mxerr);
        writelogf("Maximum binary diffusion coefficient relative error:"
                 "%12.6g", mxrelerr);
    }
}

void GasTransport::getBinDiffCorrection(double t,
        const GasTransportParams& tr, MMCollisionInt& integrals,
        size_t k, size_t j, double xk, double xj, double& fkj, double& fjk)
{
    double w1 = tr.mw[k];
    double w2 = tr.mw[j];
    double wsum = w1 + w2;
    double wmwp = (w1 - w2)/wsum;
    double sqw12 = sqrt(w1*w2);

    double sig1 = tr.sigma[k];
    double sig2 = tr.sigma[j];
    double sig12 = 0.5*(tr.sigma[k] + tr.sigma[j]);
    double sigratio = sig1*sig1/(sig2*sig2);
    double sigratio2 = sig1*sig1/(sig12*sig12);
    double sigratio3 = sig2*sig2/(sig12*sig12);

    double tstar1 = Boltzmann * t / tr.eps[k];
    double tstar2 = Boltzmann * t / tr.eps[j];
    double tstar12 = Boltzmann * t / sqrt(tr.eps[k] * tr.eps[j]);

    double om22_1 = integrals.omega22(tstar1, tr.delta(k,k));
    double om22_2 = integrals.omega22(tstar2, tr.delta(j,j));
    double om11_12 = integrals.omega11(tstar12, tr.delta(k,j));
    double astar_12 = integrals.astar(tstar12, tr.delta(k,j));
    double bstar_12 = integrals.bstar(tstar12, tr.delta(k,j));
    double cstar_12 = integrals.cstar(tstar12, tr.delta(k,j));

    double cnst = sigratio * sqrt(2.0*w2/wsum) * 2.0 * w1*w1/(wsum * w2);
    double p1 = cnst * om22_1 / om11_12;

    cnst = (1.0/sigratio) * sqrt(2.0*w1/wsum) * 2.0*w2*w2/(wsum*w1);
    double p2 = cnst * om22_2 / om11_12;

    double p12 = 15.0 * wmwp*wmwp + 8.0*w1*w2*astar_12/(wsum*wsum);

    cnst = (2.0/(w2*wsum))*sqrt(2.0*w2/wsum)*sigratio2;
    double q1 = cnst*((2.5 - 1.2*bstar_12)*w1*w1 + 3.0*w2*w2
               + 1.6*w1*w2*astar_12);

    cnst = (2.0/(w1*wsum))*sqrt(2.0*w1/wsum)*sigratio3;
    double q2 = cnst*((2.5 - 1.2*bstar_12)*w2*w2 + 3.0*w1*w1
               + 1.6*w1*w2*astar_12);

    double q12 = wmwp*wmwp*15.0*(2.5 - 1.2*bstar_12)
          + 4.0*w1*w2*astar_12*(11.0 - 2.4*bstar_12)/(wsum*wsum)
          +  1.6*wsum*om22_1*om22_2/(om11_12*om11_12*sqw12)
          * sigratio2 * sigratio3;

    cnst = 6.0*cstar_12 - 5.0;
    fkj = 1.0 + 0.1*cnst*cnst *
          (p1*xk*xk + p2*xj*xj + p12*xk*xj)/
          (q1*xk*xk + q2*xj*xj + q12*xk*xj);
    fjk = 1.0 + 0.1*cnst*cnst *
          (p2*xk*xk + p1*xj*xj + p12*xk*xj)/
          (q2*xk*xk + q1*xj*xj + q12*xk*xj);
}

}
