/**
 *  @file HighPressureGasTransport.cpp
 *  Implementation file for class HighPressureGasTransport
 *
 *  Transport parameters are calculated using corresponding states models:
 *      Binary diffusion coefficients use the generalized chart described by
 *      Takahashi, et al. and viscosity calculations use the Lucas method.
 *      All methods are described in Reid, Prausnitz, and Polling, "The Properties
 *      of Gases and Liquids, 4th ed., 1987 (viscosity in Ch. 9, Thermal
 *      conductivity in Ch. 10, and Diffusion coefficients in Ch. 11).
 **/

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/HighPressureGasTransport.h"
#include "cantera/numerics/ctlapack.h"
#include "cantera/base/utilities.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/transport/MultiTransport.h"

using namespace std;

namespace Cantera
{

HighPressureGasTransport::HighPressureGasTransport(thermo_t* thermo)
: MultiTransport(thermo)
{
}

double HighPressureGasTransport::thermalConductivity()
{
    //  Method of Ely and Hanley:
    update_T();
    doublereal Lprime_m = 0.0;
    const doublereal c1 = 1./16.04;
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    vector_fp cp_0_R(nsp);
    m_thermo->getCp_R_ref(&cp_0_R[0]);

    vector_fp L_i(nsp);
    vector_fp f_i(nsp);
    vector_fp h_i(nsp);
    vector_fp V_k(nsp);

    m_thermo -> getPartialMolarVolumes(&V_k[0]);
    doublereal L_i_min = BigNumber;

    for (size_t i = 0; i < m_nsp; i++) {
        doublereal Tc_i = Tcrit_i(i);
        doublereal Vc_i = Vcrit_i(i);
        doublereal T_r = m_thermo->temperature()/Tc_i;
        doublereal V_r = V_k[i]/Vc_i;
        doublereal T_p = std::min(T_r,2.0);
        doublereal V_p = std::max(0.5,std::min(V_r,2.0));

        // Calculate variables for density-independent component:
        doublereal theta_p = 1.0 + (m_w_ac[i] - 0.011)*(0.56553
            - 0.86276*log(T_p) - 0.69852/T_p);
        doublereal phi_p = (1.0 + (m_w_ac[i] - 0.011)*(0.38560
            - 1.1617*log(T_p)))*0.288/Zcrit_i(i);
        doublereal f_fac = Tc_i*theta_p/190.4;
        doublereal h_fac = 1000*Vc_i*phi_p/99.2;
        doublereal T_0 = m_temp/f_fac;
        doublereal mu_0 = 1e-7*(2.90774e6/T_0 - 3.31287e6*pow(T_0,-2./3.)
            + 1.60810e6*pow(T_0,-1./3.) - 4.33190e5 + 7.06248e4*pow(T_0,1./3.)
            - 7.11662e3*pow(T_0,2./3.) + 4.32517e2*T_0 - 1.44591e1*pow(T_0,4./3.)
            + 2.03712e-1*pow(T_0,5./3.));
        doublereal H = sqrt(f_fac*16.04/m_mw[i])*pow(h_fac,-2./3.);
        doublereal mu_i = mu_0*H*m_mw[i]*c1;
        L_i[i] = mu_i*1.32*GasConstant*(cp_0_R[i] - 2.5)/m_mw[i];
        L_i_min = min(L_i_min,L_i[i]);
        // Calculate variables for density-dependent component:
        doublereal theta_s = 1 + (m_w_ac[i] - 0.011)*(0.09057 - 0.86276*log(T_p)
            + (0.31664 - 0.46568/T_p)*(V_p - 0.5));
        doublereal phi_s = (1 + (m_w_ac[i] - 0.011)*(0.39490*(V_p - 1.02355)
            - 0.93281*(V_p - 0.75464)*log(T_p)))*0.288/Zcrit_i(i);
        f_i[i] = Tc_i*theta_s/190.4;
        h_i[i] = 1000*Vc_i*phi_s/99.2;
    }

    doublereal h_m = 0;
    doublereal f_m = 0;
    doublereal mw_m = 0;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Density-independent component:
            doublereal L_ij = 2*L_i[i]*L_i[j]/(L_i[i] + L_i[j] + Tiny);
            Lprime_m += molefracs[i]*molefracs[j]*L_ij;
            // Additional variables for density-dependent component:
            doublereal f_ij = sqrt(f_i[i]*f_i[j]);
            doublereal h_ij = 0.125*pow(pow(h_i[i],1./3.) + pow(h_i[j],1./3.),3.);
            doublereal mw_ij_inv = (m_mw[i] + m_mw[j])/(2*m_mw[i]*m_mw[j]);
            f_m += molefracs[i]*molefracs[j]*f_ij*h_ij;
            h_m += molefracs[i]*molefracs[j]*h_ij;
            mw_m += molefracs[i]*molefracs[j]*sqrt(mw_ij_inv*f_ij)*pow(h_ij,-4./3.);
        }
    }

    f_m = f_m/h_m;
    mw_m = pow(mw_m,-2.)*f_m*pow(h_m,-8./3.);

    doublereal rho_0 = 16.04*h_m/(1000*m_thermo->molarVolume());
    doublereal T_0 = m_temp/f_m;
    doublereal mu_0 = 1e-7*(2.90774e6/T_0 - 3.31287e6*pow(T_0,-2./3.)
                + 1.60810e6*pow(T_0,-1./3.) - 4.33190e5 + 7.06248e4
                *pow(T_0,1./3.) - 7.11662e3*pow(T_0,2./3.) + 4.32517e2*T_0
                - 1.44591e1*pow(T_0,4./3.) + 2.03712e-1*pow(T_0,5./3.));
    doublereal L_1m = 1944*mu_0;
    doublereal L_2m = (-2.5276e-4 + 3.3433e-4*pow(1.12 - log(T_0/1.680e2),2))*rho_0;
    doublereal L_3m = exp(-7.19771 + 85.67822/T_0)*(exp((12.47183
                - 984.6252*pow(T_0,-1.5))*pow(rho_0,0.1) + (rho_0/0.1617 - 1)
                *sqrt(rho_0)*(0.3594685 + 69.79841/T_0 - 872.8833*pow(T_0,-2))) - 1.)*1e-3;
    doublereal H_m = sqrt(f_m*16.04/mw_m)*pow(h_m,-2./3.);
    doublereal Lstar_m = H_m*(L_1m + L_2m + L_3m);
    return Lprime_m + Lstar_m;
}

void HighPressureGasTransport::getThermalDiffCoeffs(doublereal* const dt)
{
    // Method for MultiTransport class:
    // solveLMatrixEquation();
    // const doublereal c = 1.6/GasConstant;
    // for (size_t k = 0; k < m_nsp; k++) {
    // dt[k] = c * m_mw[k] * m_molefracs[k] * m_a[k];
    // }
    throw NotImplementedError("HighPressureGasTransport::getThermalDiffCoeffs");
}

void HighPressureGasTransport::getBinaryDiffCoeffs(const size_t ld, doublereal* const d)
{
    vector_fp PcP(5);
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    update_T();
    // Evaluate the binary diffusion coefficients from the polynomial fits.
    // This should perhaps be preceded by a check to see whether any of T, P, or
    //   C have changed.
    //if (!m_bindiff_ok) {
    updateDiff_T();
    //}
    if (ld < nsp) {
        throw CanteraError("HighPressureGasTransport::getBinaryDiffCoeffs",
                           "ld is too small");
    }
    doublereal rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = 0; j < nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            // zero (this would lead to Pr_ij = Inf):
            doublereal x_i = std::max(Tiny, molefracs[i]);
            doublereal x_j = std::max(Tiny, molefracs[j]);

            // Weight mole fractions of i and j so that X_i + X_j = 1.0:
            x_i = x_i/(x_i + x_j);
            x_j = x_j/(x_i + x_j);

            //Calculate Tr and Pr based on mole-fraction-weighted crit constants:
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            double P_corr_ij;
            if (Pr_ij < 0.1) {
                // If pressure is low enough, no correction is needed:
                P_corr_ij = 1;
            }else {
                // Otherwise, calculate the parameters for Takahashi correlation
                // by interpolating on Pr_ij:
                P_corr_ij = setPcorr(Pr_ij, Tr_ij);

                // If the reduced temperature is too low, the correction factor
                // P_corr_ij will be < 0:
                if (P_corr_ij<0) {
                    P_corr_ij = Tiny;
                }
            }

            // Multiply the standard low-pressure binary diffusion coefficient
            // (m_bdiff) by the Takahashi correction factor P_corr_ij:
            d[ld*j + i] = P_corr_ij*rp * m_bdiff(i,j);
        }
    }
}

void HighPressureGasTransport::getMultiDiffCoeffs(const size_t ld, doublereal* const d)
{
    // Not currently implemented.  m_Lmatrix inversion returns NaN.  Needs to be
    //   fixed.  --SCD - 2-28-2014
    throw NotImplementedError("HighPressureGasTransport:getMultiDiffCoeffs");
    // Calculate the multi-component Stefan-Maxwell diffusion coefficients,
    // based on the Takahashi-correlation-corrected binary diffusion coefficients.

    // update the mole fractions
    update_C();

    // update the binary diffusion coefficients
    update_T();
    updateThermal_T();

    // Correct the binary diffusion coefficients for high-pressure effects; this
    // is basically the same routine used in 'getBinaryDiffCoeffs,' above:
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    update_T();
    // Evaluate the binary diffusion coefficients from the polynomial fits -
    // this should perhaps be preceded by a check for changes in T, P, or C.
    updateDiff_T();

    if (ld < m_nsp) {
        throw CanteraError("HighPressureGasTransport::getMultiDiffCoeffs",
                           "ld is too small");
    }
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            //   zero (this would lead to Pr_ij = Inf):
            doublereal x_i = std::max(Tiny, molefracs[i]);
            doublereal x_j = std::max(Tiny, molefracs[j]);
            x_i = x_i/(x_i+x_j);
            x_j = x_j/(x_i+x_j);
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            double P_corr_ij;
            if (Pr_ij < 0.1) {
                P_corr_ij = 1;
            }else {
                P_corr_ij = setPcorr(Pr_ij, Tr_ij);
                if (P_corr_ij<0) {
                    P_corr_ij = Tiny;
                }
            }

            m_bdiff(i,j) *= P_corr_ij;
        }
    }
    m_bindiff_ok = false; // m_bdiff is overwritten by the above routine.

    // Having corrected m_bdiff for pressure and concentration effects, the
    //    routine now proceeds the same as in the low-pressure case:

    // evaluate L0000 if the temperature or concentrations have
    // changed since it was last evaluated.
    if (!m_l0000_ok) {
        eval_L0000(molefracs.data());
    }

    // invert L00,00
    int ierr = invert(m_Lmatrix, m_nsp);
    if (ierr != 0) {
        throw CanteraError("HighPressureGasTransport::getMultiDiffCoeffs",
                           "invert returned ierr = {}", ierr);
    }
    m_l0000_ok = false; // matrix is overwritten by inverse
    m_lmatrix_soln_ok = false;

    doublereal prefactor = 16.0 * m_temp
        *m_thermo->meanMolecularWeight()/(25.0*m_thermo->pressure());

    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            double c = prefactor/m_mw[j];
            d[ld*j + i] = c*molefracs[i]*(m_Lmatrix(i,j) - m_Lmatrix(i,i));
        }
    }
}

doublereal HighPressureGasTransport::viscosity()
{
    // Calculate the high-pressure mixture viscosity, based on the Lucas method.
    double Tc_mix = 0.;
    double Pc_mix_n = 0.;
    double Pc_mix_d = 0.;
    double MW_mix = m_thermo->meanMolecularWeight();
    double MW_H = m_mw[0];
    double MW_L = m_mw[0];
    doublereal FP_mix_o = 0;
    doublereal FQ_mix_o = 0;
    doublereal tKelvin = m_thermo->temperature();
    double Pvp_mix = m_thermo->satPressure(tKelvin);
    size_t nsp = m_thermo->nSpecies();
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    double x_H = molefracs[0];
    for (size_t i = 0; i < m_nsp; i++) {
        // Calculate pure-species critical constants and add their contribution
        // to the mole-fraction-weighted mixture averages:
        double Tc = Tcrit_i(i);
        double Tr = tKelvin/Tc;
        double Zc = Zcrit_i(i);
        Tc_mix += Tc*molefracs[i];
        Pc_mix_n += molefracs[i]*Zc; //numerator
        Pc_mix_d += molefracs[i]*Vcrit_i(i); //denominator

        // Need to calculate ratio of heaviest to lightest species:
        if (m_mw[i] > MW_H) {
            MW_H = m_mw[i];
            x_H = molefracs[i];
        } else if (m_mw[i] < MW_L) {
            MW_L = m_mw[i];        }

        // Calculate reduced dipole moment for polar correction term:
        doublereal mu_ri = 52.46*100000*m_dipole(i,i)*m_dipole(i,i)
            *Pcrit_i(i)/(Tc*Tc);
        if (mu_ri < 0.022) {
            FP_mix_o += molefracs[i];
        } else if (mu_ri < 0.075) {
            FP_mix_o += molefracs[i]*(1. + 30.55*pow(0.292 - Zc, 1.72));
        } else { FP_mix_o += molefracs[i]*(1. + 30.55*pow(0.292 - Zc, 1.72)
                                    *fabs(0.96 + 0.1*(Tr - 0.7)));
        }

        // Calculate contribution to quantum correction term.
        // SCD Note:  This assumes the species of interest (He, H2, and D2) have
        //   been named in this specific way.  They are perhaps the most obvious
        //   names, but it would of course be preferred to have a more general
        //   approach, here.
        std::vector<std::string> spnames = m_thermo->speciesNames();
        if (spnames[i] == "He") {
            FQ_mix_o += molefracs[i]*FQ_i(1.38,Tr,m_mw[i]);
        } else if (spnames[i] == "H2") {
            FQ_mix_o += molefracs[i]*(FQ_i(0.76,Tr,m_mw[i]));
        } else if (spnames[i] == "D2") {
            FQ_mix_o += molefracs[i]*(FQ_i(0.52,Tr,m_mw[i]));
        } else {
            FQ_mix_o += molefracs[i];
        }
    }

    double Tr_mix = tKelvin/Tc_mix;
    double Pc_mix = GasConstant*Tc_mix*Pc_mix_n/Pc_mix_d;
    double Pr_mix = m_thermo->pressure()/Pc_mix;
    double ratio = MW_H/MW_L;
    double ksi = pow(GasConstant*Tc_mix*3.6277*pow(10.0,53.0)/(pow(MW_mix,3)
                        *pow(Pc_mix,4)),1.0/6.0);

    if (ratio > 9 && x_H > 0.05 && x_H < 0.7) {
        FQ_mix_o *= 1 - 0.01*pow(ratio,0.87);
    }

    // Calculate Z1m
    double Z1m = (0.807*pow(Tr_mix,0.618) - 0.357*exp(-0.449*Tr_mix)
                 + 0.340*exp(-4.058*Tr_mix)+0.018)*FP_mix_o*FQ_mix_o;

    // Calculate Z2m:
    double Z2m;
    if (Tr_mix <= 1.0) {
        if (Pr_mix < Pvp_mix/Pc_mix) {
            doublereal alpha = 3.262 + 14.98*pow(Pr_mix,5.508);
            doublereal beta = 1.390 + 5.746*Pr_mix;
            Z2m = 0.600 + 0.760*pow(Pr_mix,alpha) + (0.6990*pow(Pr_mix,beta) -
                0.60)*(1- Tr_mix);
        } else {
            throw CanteraError("HighPressureGasTransport::viscosity",
                               "State is outside the limits of the Lucas model, Tr <= 1");
        }
    } else if ((Tr_mix > 1.0) && (Tr_mix < 40.0)) {
        if ((Pr_mix > 0.0) && (Pr_mix <= 100.0)) {
            doublereal a_fac = 0.001245*exp(5.1726*pow(Tr_mix,-0.3286))/Tr_mix;
            doublereal b_fac = a_fac*(1.6553*Tr_mix - 1.2723);
            doublereal c_fac = 0.4489*exp(3.0578*pow(Tr_mix,-37.7332))/Tr_mix;
            doublereal d_fac = 1.7368*exp(2.2310*pow(Tr_mix,-7.6351))/Tr_mix;
            doublereal f_fac = 0.9425*exp(-0.1853*pow(Tr_mix,0.4489));

            Z2m = Z1m*(1 + a_fac*pow(Pr_mix,1.3088)/(b_fac*pow(Pr_mix,f_fac)
                        + pow(1+c_fac*pow(Pr_mix,d_fac),-1)));
        } else {
            throw CanteraError("HighPressureGasTransport::viscosity",
                           "State is outside the limits of the Lucas model, 1.0 < Tr < 40");
        }
    } else {
        throw CanteraError("HighPressureGasTransport::viscosity",
                           "State is outside the limits of the Lucas model, Tr > 40");
    }

    // Calculate Y:
    doublereal Y = Z2m/Z1m;

    // Return the viscosity:
    return Z2m*(1 + (FP_mix_o - 1)*pow(Y,-3))*(1 + (FQ_mix_o - 1)
            *(1/Y - 0.007*pow(log(Y),4)))/(ksi*FP_mix_o*FQ_mix_o);
}

// Pure species critical properties - Tc, Pc, Vc, Zc:
doublereal HighPressureGasTransport::Tcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector_fp molefracs = store(i, m_thermo->nSpecies());

    double tc = m_thermo->critTemperature();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return tc;
}

doublereal HighPressureGasTransport::Pcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector_fp molefracs = store(i, m_thermo->nSpecies());

    double pc = m_thermo->critPressure();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return pc;
}

doublereal HighPressureGasTransport::Vcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector_fp molefracs = store(i, m_thermo->nSpecies());

    double vc = m_thermo->critVolume();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return vc;
}

doublereal HighPressureGasTransport::Zcrit_i(size_t i)
{
    // Store current molefracs and set temp molefrac of species i to 1.0:
    vector_fp molefracs = store(i, m_thermo->nSpecies());

    double zc = m_thermo->critCompressibility();
    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return zc;
}

vector_fp HighPressureGasTransport::store(size_t i, size_t nsp)
{
    vector_fp molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    vector_fp mf_temp(nsp, 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);
    return molefracs;
}

// Calculates quantum correction term for a species based on Tr and MW, used in
//   viscosity calculation:
doublereal HighPressureGasTransport::FQ_i(doublereal Q, doublereal Tr, doublereal MW)
{
    return 1.22*pow(Q,0.15)*(1 + 0.00385*pow(pow(Tr - 12.,2.),1./MW)
                             *fabs(Tr-12)/(Tr-12));
}

// Set value of parameter values for Takahashi correlation, by interpolating
//   table of constants vs. Pr:
doublereal HighPressureGasTransport::setPcorr(doublereal Pr, doublereal Tr)
{
    const static double Pr_lookup[17] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0,
        1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0};
    const static double DP_Rt_lookup[17] = {1.01, 1.01, 1.01, 1.01, 1.01, 1.01,
        1.01, 1.02, 1.02, 1.02, 1.02, 1.03, 1.03, 1.04, 1.05, 1.06, 1.07};
    const static double A_ij_lookup[17] = {0.038042, 0.067433, 0.098317,
        0.137610, 0.175081, 0.216376, 0.314051, 0.385736, 0.514553, 0.599184,
        0.557725, 0.593007, 0.696001, 0.790770, 0.502100, 0.837452, 0.890390};
    const static double B_ij_lookup[17] = {1.52267, 2.16794, 2.42910, 2.77605,
        2.98256, 3.11384, 3.50264, 3.07773, 3.54744, 3.61216, 3.41882, 3.18415,
        3.37660, 3.27984, 3.39031, 3.23513, 3.13001};
    const static double C_ij_lookup[17] = {0., 0., 0., 0., 0., 0., 0., 0.141211,
        0.278407, 0.372683, 0.504894, 0.678469, 0.665702, 0., 0.602907, 0., 0.};
    const static double E_ij_lookup[17] = {1., 1., 1., 1., 1., 1., 1., 13.45454,
        14., 10.00900, 8.57519, 10.37483, 11.21674, 1., 6.19043, 1., 1.};

    // Interpolate Pr vs. those used in Takahashi table:
    int Pr_i = 0;
    double frac = 0.;

    if (Pr < 0.1) {
        frac = (Pr - Pr_lookup[0])/(Pr_lookup[1] - Pr_lookup[0]);
    } else {
        for (int j = 1; j < 17; j++) {
            if (Pr_lookup[j] > Pr) {
                frac = (Pr - Pr_lookup[j-1])/(Pr_lookup[j] - Pr_lookup[j-1]);
                break;
            }
            Pr_i++;
        }
    }
    // If Pr is greater than the greatest value used by Takahashi (5.0), use the
    //   final table value.  Should eventually add in an extrapolation:
    if (Pr_i == 17) {
        frac = 1.0;
    }

    doublereal P_corr_1 = DP_Rt_lookup[Pr_i]*(1.0 - A_ij_lookup[Pr_i]
        *pow(Tr,-B_ij_lookup[Pr_i]))*(1-C_ij_lookup[Pr_i]
        *pow(Tr,-E_ij_lookup[Pr_i]));
    doublereal P_corr_2 = DP_Rt_lookup[Pr_i+1]*(1.0 - A_ij_lookup[Pr_i+1]
        *pow(Tr,-B_ij_lookup[Pr_i+1]))*(1-C_ij_lookup[Pr_i+1]
        *pow(Tr,-E_ij_lookup[Pr_i+1]));
    return P_corr_1*(1.0-frac) + P_corr_2*frac;
}

}
