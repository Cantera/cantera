//! @file GasTransport.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/GasTransport.h"
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
    m_viscmix(0.0),
    m_visc_ok(false),
    m_viscwt_ok(false),
    m_spvisc_ok(false),
    m_bindiff_ok(false),
    m_mode(0),
    m_polytempvec(5),
    m_temp(-1.0),
    m_kbt(0.0),
    m_sqrt_kbt(0.0),
    m_sqrt_t(0.0),
    m_logt(0.0),
    m_t14(0.0),
    m_t32(0.0),
    m_log_level(0)
{
}

void GasTransport::update_T()
{
    if (m_thermo->nSpecies() != m_nsp) {
        // Rebuild data structures if number of species has changed
        init(m_thermo, m_mode, m_log_level);
    }

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

    multiply(m_phi, m_molefracs.data(), m_spwork.data());

    for (size_t k = 0; k < m_nsp; k++) {
        vismix += m_molefracs[k] * m_visc[k]/m_spwork[k]; //denom;
    }
    m_viscmix = vismix;
    return vismix;
}

void GasTransport::updateViscosity_T()
{
    if (!m_spvisc_ok) {
        updateSpeciesViscosities();
    }

    // see Eq. (9-5.15) of Reid, Prausnitz, and Poling
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t k = j; k < m_nsp; k++) {
            double vratiokj = m_visc[k]/m_visc[j];
            double wratiojk = m_mw[j]/m_mw[k];

            // Note that m_wratjk(k,j) holds the square root of m_wratjk(j,k)!
            double factor1 = 1.0 + (m_sqvisc[k]/m_sqvisc[j]) * m_wratjk(k,j);
            m_phi(k,j) = factor1*factor1 / (sqrt(8.0) * m_wratkj1(j,k));
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
        throw CanteraError("GasTransport::getBinaryDiffCoeffs", "ld is too small");
    }
    doublereal rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            d[ld*j + i] = rp * m_bdiff(i,j);
        }
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
            d[k] = 1.0 / (sum1 + sum2);
        }
    }
}

void GasTransport::init(thermo_t* thermo, int mode, int log_level)
{
    m_thermo = thermo;
    m_nsp = m_thermo->nSpecies();
    m_mode = mode;
    m_log_level = log_level;

    // set up Monchick and Mason collision integrals
    setupCollisionParameters();
    setupCollisionIntegral();

    m_molefracs.resize(m_nsp);
    m_spwork.resize(m_nsp);
    m_visc.resize(m_nsp);
    m_sqvisc.resize(m_nsp);
    m_phi.resize(m_nsp, m_nsp, 0.0);
    m_bdiff.resize(m_nsp, m_nsp);

    // make a local copy of the molecular weights
    m_mw = m_thermo->molecularWeights();

    m_wratjk.resize(m_nsp, m_nsp, 0.0);
    m_wratkj1.resize(m_nsp, m_nsp, 0.0);
    for (size_t j = 0; j < m_nsp; j++) {
        for (size_t k = j; k < m_nsp; k++) {
            m_wratjk(j,k) = sqrt(m_mw[j]/m_mw[k]);
            m_wratjk(k,j) = sqrt(m_wratjk(j,k));
            m_wratkj1(j,k) = sqrt(1.0 + m_mw[k]/m_mw[j]);
        }
    }
}

void GasTransport::setupCollisionParameters()
{
    m_epsilon.resize(m_nsp, m_nsp, 0.0);
    m_delta.resize(m_nsp, m_nsp, 0.0);
    m_reducedMass.resize(m_nsp, m_nsp, 0.0);
    m_dipole.resize(m_nsp, m_nsp, 0.0);
    m_diam.resize(m_nsp, m_nsp, 0.0);
    m_crot.resize(m_nsp);
    m_zrot.resize(m_nsp);
    m_polar.resize(m_nsp, false);
    m_alpha.resize(m_nsp, 0.0);
    m_poly.resize(m_nsp);
    m_sigma.resize(m_nsp);
    m_eps.resize(m_nsp);
    m_w_ac.resize(m_nsp);
    m_disp.resize(m_nsp, 0.0);
    m_quad_polar.resize(m_nsp, 0.0);

    const vector_fp& mw = m_thermo->molecularWeights();
    getTransportData();

    for (size_t i = 0; i < m_nsp; i++) {
        m_poly[i].resize(m_nsp);
    }

    double f_eps, f_sigma;

    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            // the reduced mass
            m_reducedMass(i,j) = mw[i] * mw[j] / (Avogadro * (mw[i] + mw[j]));

            // hard-sphere diameter for (i,j) collisions
            m_diam(i,j) = 0.5*(m_sigma[i] + m_sigma[j]);

            // the effective well depth for (i,j) collisions
            m_epsilon(i,j) = sqrt(m_eps[i]*m_eps[j]);

            // the effective dipole moment for (i,j) collisions
            m_dipole(i,j) = sqrt(m_dipole(i,i)*m_dipole(j,j));

            // reduced dipole moment delta* (nondimensional)
            double d = m_diam(i,j);
            m_delta(i,j) = 0.5 * m_dipole(i,j)*m_dipole(i,j)
                           / (4 * Pi * epsilon_0 * m_epsilon(i,j) * d * d * d);
            makePolarCorrections(i, j, f_eps, f_sigma);
            m_diam(i,j) *= f_sigma;
            m_epsilon(i,j) *= f_eps;

            // properties are symmetric
            m_reducedMass(j,i) = m_reducedMass(i,j);
            m_diam(j,i) = m_diam(i,j);
            m_epsilon(j,i) = m_epsilon(i,j);
            m_dipole(j,i) = m_dipole(i,j);
            m_delta(j,i) = m_delta(i,j);
        }
    }
}

void GasTransport::setupCollisionIntegral()
{
    double tstar_min = 1.e8, tstar_max = 0.0;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            // The polynomial fits of collision integrals vs. T*
            // will be done for the T* from tstar_min to tstar_max
            tstar_min = std::min(tstar_min, Boltzmann * m_thermo->minTemp()/m_epsilon(i,j));
            tstar_max = std::max(tstar_max, Boltzmann * m_thermo->maxTemp()/m_epsilon(i,j));
        }
    }
    // Chemkin fits the entire T* range in the Monchick and Mason tables,
    // so modify tstar_min and tstar_max if in Chemkin compatibility mode
    if (m_mode == CK_Mode) {
        tstar_min = 0.101;
        tstar_max = 99.9;
    }

    // initialize the collision integral calculator for the desired T* range
    debuglog("*** collision_integrals ***\n", m_log_level);
    MMCollisionInt integrals;
    integrals.init(tstar_min, tstar_max, m_log_level);
    fitCollisionIntegrals(integrals);
    debuglog("*** end of collision_integrals ***\n", m_log_level);
    // make polynomial fits
    debuglog("*** property fits ***\n", m_log_level);
    fitProperties(integrals);
    debuglog("*** end of property fits ***\n", m_log_level);
}

void GasTransport::getTransportData()
{
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        shared_ptr<Species> s = m_thermo->species(k);
        const GasTransportData* sptran =
            dynamic_cast<GasTransportData*>(s->transport.get());
        if (!sptran) {
            throw CanteraError("GasTransport::getTransportData",
                "Missing gas-phase transport data for species '{}'.", s->name);
        }

        if (sptran->geometry == "atom") {
            m_crot[k] = 0.0;
        } else if (sptran->geometry == "linear") {
            m_crot[k] = 1.0;
        } else if (sptran->geometry == "nonlinear") {
            m_crot[k] = 1.5;
        }

        m_sigma[k] = sptran->diameter;
        m_eps[k] = sptran->well_depth;
        m_dipole(k,k) = sptran->dipole;
        m_polar[k] = (sptran->dipole > 0);
        m_alpha[k] = sptran->polarizability;
        m_zrot[k] = sptran->rotational_relaxation;
        m_w_ac[k] = sptran->acentric_factor;
        m_disp[k] = sptran->dispersion_coefficient;
        m_quad_polar[k] = sptran->quadrupole_polarizability;
    }
}

void GasTransport::makePolarCorrections(size_t i, size_t j,
        doublereal& f_eps, doublereal& f_sigma)
{
    // no correction if both are nonpolar, or both are polar
    if (m_polar[i] == m_polar[j]) {
        f_eps = 1.0;
        f_sigma = 1.0;
        return;
    }

    // corrections to the effective diameter and well depth
    // if one is polar and one is non-polar
    size_t kp = (m_polar[i] ? i : j); // the polar one
    size_t knp = (i == kp ? j : i); // the nonpolar one
    double d3np, d3p, alpha_star, mu_p_star, xi;
    d3np = pow(m_sigma[knp],3);
    d3p = pow(m_sigma[kp],3);
    alpha_star = m_alpha[knp]/d3np;
    mu_p_star = m_dipole(kp,kp)/sqrt(4 * Pi * epsilon_0 * d3p * m_eps[kp]);
    xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star *
         sqrt(m_eps[kp]/m_eps[knp]);
    f_sigma = pow(xi, -1.0/6.0);
    f_eps = xi*xi;
}

void GasTransport::fitCollisionIntegrals(MMCollisionInt& integrals)
{
    // Chemkin fits to sixth order polynomials
    int degree = (m_mode == CK_Mode ? 6 : COLL_INT_POLY_DEGREE);
    if (m_log_level) {
        writelog("tstar_fits\n"
                 "fits to A*, B*, and C* vs. log(T*).\n"
                 "These are done only for the required dstar(j,k) values.\n\n");
        if (m_log_level < 3) {
            writelog("*** polynomial coefficients not printed (log_level < 3) ***\n");
        }
    }
    vector_fp fitlist;
    m_omega22_poly.clear();
    m_astar_poly.clear();
    m_bstar_poly.clear();
    m_cstar_poly.clear();
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = i; j < m_nsp; j++) {
            // Chemkin fits only delta* = 0
            double dstar = (m_mode != CK_Mode) ? m_delta(i,j) : 0.0;

            // if a fit has already been generated for delta* = m_delta(i,j),
            // then use it. Otherwise, make a new fit, and add m_delta(i,j) to
            // the list of delta* values for which fits have been done.

            // 'find' returns a pointer to end() if not found
            auto dptr = find(fitlist.begin(), fitlist.end(), dstar);
            if (dptr == fitlist.end()) {
                vector_fp ca(degree+1), cb(degree+1), cc(degree+1);
                vector_fp co22(degree+1);
                integrals.fit(degree, dstar, ca.data(), cb.data(), cc.data());
                integrals.fit_omega22(degree, dstar, co22.data());
                m_omega22_poly.push_back(co22);
                m_astar_poly.push_back(ca);
                m_bstar_poly.push_back(cb);
                m_cstar_poly.push_back(cc);
                m_poly[i][j] = static_cast<int>(m_astar_poly.size()) - 1;
                fitlist.push_back(dstar);
            } else {
                // delta* found in fitlist, so just point to this polynomial
                m_poly[i][j] = static_cast<int>((dptr - fitlist.begin()));
            }
            m_poly[j][i] = m_poly[i][j];
        }
    }
}

void GasTransport::fitProperties(MMCollisionInt& integrals)
{
    // number of points to use in generating fit data
    const size_t np = 50;
    int degree = (m_mode == CK_Mode ? 3 : 4);
    double dt = (m_thermo->maxTemp() - m_thermo->minTemp())/(np-1);
    vector_fp tlog(np), spvisc(np), spcond(np);
    vector_fp w(np), w2(np);

    m_visccoeffs.clear();
    m_condcoeffs.clear();

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        double t = m_thermo->minTemp() + dt*n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector_fp c(degree + 1), c2(degree + 1);

    // fit the pure-species viscosity and thermal conductivity for each species
    if (m_log_level && m_log_level < 2) {
        writelog("*** polynomial coefficients not printed (log_level < 2) ***\n");
    }
    double visc, err, relerr,
               mxerr = 0.0, mxrelerr = 0.0, mxerr_cond = 0.0, mxrelerr_cond = 0.0;

    if (m_log_level) {
        writelog("Polynomial fits for viscosity:\n");
        if (m_mode == CK_Mode) {
            writelog("log(viscosity) fit to cubic polynomial in log(T)\n");
        } else {
            writelogf("viscosity/sqrt(T) fit to polynomial of degree "
                      "%d in log(T)", degree);
        }
    }

    double T_save = m_thermo->temperature();
    const vector_fp& mw = m_thermo->molecularWeights();
    for (size_t k = 0; k < m_nsp; k++) {
        double tstar = Boltzmann * 298.0 / m_eps[k];
        // Scaling factor for temperature dependence of z_rot. [Kee2003] Eq.
        // 12.112 or [Kee2017] Eq. 11.115
        double fz_298 = 1.0 + pow(Pi, 1.5) / sqrt(tstar) * (0.5 + 1.0 / tstar) +
            (0.25 * Pi * Pi + 2) / tstar;

        for (size_t n = 0; n < np; n++) {
            double t = m_thermo->minTemp() + dt*n;
            m_thermo->setTemperature(t);
            vector_fp cp_R_all(m_thermo->nSpecies());
            m_thermo->getCp_R_ref(&cp_R_all[0]);
            double cp_R = cp_R_all[k];
            tstar = Boltzmann * t / m_eps[k];
            double sqrt_T = sqrt(t);
            double om22 = integrals.omega22(tstar, m_delta(k,k));
            double om11 = integrals.omega11(tstar, m_delta(k,k));

            // self-diffusion coefficient, without polar corrections
            double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/m_reducedMass(k,k)) *
                               pow((Boltzmann * t), 1.5)/
                               (Pi * m_sigma[k] * m_sigma[k] * om11);

            // viscosity
            visc = 5.0/16.0 * sqrt(Pi * mw[k] * Boltzmann * t / Avogadro) /
                   (om22 * Pi * m_sigma[k]*m_sigma[k]);

            // thermal conductivity
            double f_int = mw[k]/(GasConstant * t) * diffcoeff/visc;
            double cv_rot = m_crot[k];
            double A_factor = 2.5 - f_int;
            double fz_tstar = 1.0 + pow(Pi, 1.5) / sqrt(tstar) * (0.5 + 1.0 / tstar) +
                (0.25 * Pi * Pi + 2) / tstar;
            double B_factor = m_zrot[k] * fz_298 / fz_tstar + 2.0/Pi * (5.0/3.0 * cv_rot + f_int);
            double c1 = 2.0/Pi * A_factor/B_factor;
            double cv_int = cp_R - 2.5 - cv_rot;
            double f_rot = f_int * (1.0 + c1);
            double f_trans = 2.5 * (1.0 - c1 * cv_rot/1.5);
            double cond = (visc/mw[k])*GasConstant*(f_trans * 1.5
                                                + f_rot * cv_rot + f_int * cv_int);

            if (m_mode == CK_Mode) {
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
        polyfit(np, degree, tlog.data(), spvisc.data(), w.data(), c.data());
        polyfit(np, degree, tlog.data(), spcond.data(), w2.data(), c2.data());

        // evaluate max fit errors for viscosity
        for (size_t n = 0; n < np; n++) {
            double val, fit;
            if (m_mode == CK_Mode) {
                val = exp(spvisc[n]);
                fit = exp(poly3(tlog[n], c.data()));
            } else {
                double sqrt_T = exp(0.5*tlog[n]);
                val = sqrt_T * pow(spvisc[n],2);
                fit = sqrt_T * pow(poly4(tlog[n], c.data()),2);
            }
            err = fit - val;
            relerr = err/val;
            mxerr = std::max(mxerr, fabs(err));
            mxrelerr = std::max(mxrelerr, fabs(relerr));
        }

        // evaluate max fit errors for conductivity
        for (size_t n = 0; n < np; n++) {
            double val, fit;
            if (m_mode == CK_Mode) {
                val = exp(spcond[n]);
                fit = exp(poly3(tlog[n], c2.data()));
            } else {
                double sqrt_T = exp(0.5*tlog[n]);
                val = sqrt_T * spcond[n];
                fit = sqrt_T * poly4(tlog[n], c2.data());
            }
            err = fit - val;
            relerr = err/val;
            mxerr_cond = std::max(mxerr_cond, fabs(err));
            mxrelerr_cond = std::max(mxrelerr_cond, fabs(relerr));
        }
        m_visccoeffs.push_back(c);
        m_condcoeffs.push_back(c2);

        if (m_log_level >= 2) {
            writelog(m_thermo->speciesName(k) + ": [" + vec2str(c) + "]\n");
        }
    }
    m_thermo->setTemperature(T_save);

    if (m_log_level) {
        writelogf("Maximum viscosity absolute error:  %12.6g\n", mxerr);
        writelogf("Maximum viscosity relative error:  %12.6g\n", mxrelerr);
        writelog("\nPolynomial fits for conductivity:\n");
        if (m_mode == CK_Mode) {
            writelog("log(conductivity) fit to cubic polynomial in log(T)");
        } else {
            writelogf("conductivity/sqrt(T) fit to "
                      "polynomial of degree %d in log(T)", degree);
        }
        if (m_log_level >= 2) {
            for (size_t k = 0; k < m_nsp; k++) {
                writelog(m_thermo->speciesName(k) + ": [" +
                         vec2str(m_condcoeffs[k]) + "]\n");
            }
        }
        writelogf("Maximum conductivity absolute error:  %12.6g\n", mxerr_cond);
        writelogf("Maximum conductivity relative error:  %12.6g\n", mxrelerr_cond);

        // fit the binary diffusion coefficients for each species pair
        writelogf("\nbinary diffusion coefficients:\n");
        if (m_mode == CK_Mode) {
            writelog("log(D) fit to cubic polynomial in log(T)");
        } else {
            writelogf("D/T**(3/2) fit to polynomial of degree %d in log(T)",degree);
        }
    }

    fitDiffCoeffs(integrals);
}

void GasTransport::fitDiffCoeffs(MMCollisionInt& integrals)
{
    // number of points to use in generating fit data
    const size_t np = 50;
    int degree = (m_mode == CK_Mode ? 3 : 4);
    double dt = (m_thermo->maxTemp() - m_thermo->minTemp())/(np-1);
    vector_fp tlog(np);
    vector_fp w(np), w2(np);

    // generate array of log(t) values
    for (size_t n = 0; n < np; n++) {
        double t = m_thermo->minTemp() + dt*n;
        tlog[n] = log(t);
    }

    // vector of polynomial coefficients
    vector_fp c(degree + 1), c2(degree + 1);
    double err, relerr,
               mxerr = 0.0, mxrelerr = 0.0;

    vector_fp diff(np + 1);
    m_diffcoeffs.clear();
    for (size_t k = 0; k < m_nsp; k++) {
        for (size_t j = k; j < m_nsp; j++) {
            for (size_t n = 0; n < np; n++) {
                double t = m_thermo->minTemp() + dt*n;
                double eps = m_epsilon(j,k);
                double tstar = Boltzmann * t/eps;
                double sigma = m_diam(j,k);
                double om11 = integrals.omega11(tstar, m_delta(j,k));
                double diffcoeff = 3.0/16.0 * sqrt(2.0 * Pi/m_reducedMass(k,j))
                    * pow(Boltzmann * t, 1.5) / (Pi * sigma * sigma * om11);

                // 2nd order correction
                // NOTE: THIS CORRECTION IS NOT APPLIED
                double fkj, fjk;
                getBinDiffCorrection(t, integrals, k, j, 1.0, 1.0, fkj, fjk);

                if (m_mode == CK_Mode) {
                    diff[n] = log(diffcoeff);
                    w[n] = -1.0;
                } else {
                    diff[n] = diffcoeff/pow(t, 1.5);
                    w[n] = 1.0/(diff[n]*diff[n]);
                }
            }
            polyfit(np, degree, tlog.data(), diff.data(), w.data(), c.data());

            for (size_t n = 0; n < np; n++) {
                double val, fit;
                if (m_mode == CK_Mode) {
                    val = exp(diff[n]);
                    fit = exp(poly3(tlog[n], c.data()));
                } else {
                    double t = exp(tlog[n]);
                    double pre = pow(t, 1.5);
                    val = pre * diff[n];
                    fit = pre * poly4(tlog[n], c.data());
                }
                err = fit - val;
                relerr = err/val;
                mxerr = std::max(mxerr, fabs(err));
                mxrelerr = std::max(mxrelerr, fabs(relerr));
            }
            m_diffcoeffs.push_back(c);
            if (m_log_level >= 2) {
                writelog(m_thermo->speciesName(k) + "__" +
                         m_thermo->speciesName(j) + ": [" + vec2str(c) + "]\n");
            }
        }
    }
    if (m_log_level) {
        writelogf("Maximum binary diffusion coefficient absolute error:"
                 "  %12.6g\n", mxerr);
        writelogf("Maximum binary diffusion coefficient relative error:"
                 "%12.6g", mxrelerr);
    }
}

void GasTransport::getBinDiffCorrection(double t, MMCollisionInt& integrals,
        size_t k, size_t j, double xk, double xj, double& fkj, double& fjk)
{
    double w1 = m_thermo->molecularWeight(k);
    double w2 = m_thermo->molecularWeight(j);
    double wsum = w1 + w2;
    double wmwp = (w1 - w2)/wsum;
    double sqw12 = sqrt(w1*w2);
    double sig1 = m_sigma[k];
    double sig2 = m_sigma[j];
    double sig12 = 0.5*(m_sigma[k] + m_sigma[j]);
    double sigratio = sig1*sig1/(sig2*sig2);
    double sigratio2 = sig1*sig1/(sig12*sig12);
    double sigratio3 = sig2*sig2/(sig12*sig12);
    double tstar1 = Boltzmann * t / m_eps[k];
    double tstar2 = Boltzmann * t / m_eps[j];
    double tstar12 = Boltzmann * t / sqrt(m_eps[k] * m_eps[j]);
    double om22_1 = integrals.omega22(tstar1, m_delta(k,k));
    double om22_2 = integrals.omega22(tstar2, m_delta(j,j));
    double om11_12 = integrals.omega11(tstar12, m_delta(k,j));
    double astar_12 = integrals.astar(tstar12, m_delta(k,j));
    double bstar_12 = integrals.bstar(tstar12, m_delta(k,j));
    double cstar_12 = integrals.cstar(tstar12, m_delta(k,j));

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
