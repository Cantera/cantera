/**
 *  @file HighPressureGasTransport.cpp
 *  Implementation file for class HighPressureGasTransport
 **/

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/transport/HighPressureGasTransport.h"
#include "cantera/base/utilities.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/transport/MultiTransport.h"

using namespace std;

namespace Cantera
{

/**
 * @brief Returns interpolated value of (DP)_R obtained from the data
 * in Table 2 of the Takahashi 1975 paper, given a value of the reduced
 * pressure (Pr) and reduced temperature (Tr).
 *
 * @param Pr  Reduced pressure
 * @param Tr  Reduced temperature
 */
double takahashi_correction_factor(double Pr, double Tr)
{
    // In the low pressure limit, no correction is needed. Interpolate
    // the value towards 1 as pressure drops below the 0.1 threshold.
    if (Pr < 0.1) {
        return 1.0;
    }

    // Data from Table 2 of Takahashi 1975 paper:
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

    // Interpolate to obtain the value of the constants (DP)_R, A, B, C, E at
    // the provided value of the reduced pressure (Pr).
    int Pr_i = 0;
    double frac = 0.0;
    double A, B, C,E, DP_Rt;

    for (int j = 1; j < 17; j++){
        if (Pr_lookup[j] > Pr) {
            frac = (Pr - Pr_lookup[j-1])/(Pr_lookup[j] - Pr_lookup[j-1]);
            break;
        }
        Pr_i++;
    }
    // If this loop completes without finding a bounding value of Pr, use
    // the final table value.
    frac = 1.0;


    DP_Rt = DP_Rt_lookup[Pr_i]*(1.0 - frac) + DP_Rt_lookup[Pr_i+1]*frac;
    A = A_ij_lookup[Pr_i]*(1.0 - frac) + A_ij_lookup[Pr_i+1]*frac;
    B = B_ij_lookup[Pr_i]*(1.0 - frac) + B_ij_lookup[Pr_i+1]*frac;
    C = C_ij_lookup[Pr_i]*(1.0 - frac) + C_ij_lookup[Pr_i+1]*frac;
    E = E_ij_lookup[Pr_i]*(1.0 - frac) + E_ij_lookup[Pr_i+1]*frac;

    double DP_R = DP_Rt*(1.0 - A*pow(Tr,-B))*(1.0 - C*pow(Tr,-E));
    return DP_R;
}



double HighPressureGasTransport::thermalConductivity()
{
    //  Method of Ely and Hanley:
    update_T();
    double Lprime_m = 0.0;
    const double c1 = 1./16.04;
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> cp_0_R(nsp);
    m_thermo->getCp_R_ref(&cp_0_R[0]);

    vector<double> L_i(nsp);
    vector<double> f_i(nsp);
    vector<double> h_i(nsp);
    vector<double> V_k(nsp);

    m_thermo -> getPartialMolarVolumes(&V_k[0]);
    double L_i_min = BigNumber;

    for (size_t i = 0; i < m_nsp; i++) {
        double Tc_i = Tcrit_i(i);
        double Vc_i = Vcrit_i(i);
        double T_r = m_thermo->temperature()/Tc_i;
        double V_r = V_k[i]/Vc_i;
        double T_p = std::min(T_r,2.0);
        double V_p = std::max(0.5,std::min(V_r,2.0));

        // Calculate variables for density-independent component:
        double theta_p = 1.0 + (m_w_ac[i] - 0.011)*(0.56553
            - 0.86276*log(T_p) - 0.69852/T_p);
        double phi_p = (1.0 + (m_w_ac[i] - 0.011)*(0.38560
            - 1.1617*log(T_p)))*0.288/Zcrit_i(i);
        double f_fac = Tc_i*theta_p/190.4;
        double h_fac = 1000*Vc_i*phi_p/99.2;
        double T_0 = m_temp/f_fac;
        double mu_0 = 1e-7*(2.90774e6/T_0 - 3.31287e6*pow(T_0,-2./3.)
            + 1.60810e6*pow(T_0,-1./3.) - 4.33190e5 + 7.06248e4*pow(T_0,1./3.)
            - 7.11662e3*pow(T_0,2./3.) + 4.32517e2*T_0 - 1.44591e1*pow(T_0,4./3.)
            + 2.03712e-1*pow(T_0,5./3.));
        double H = sqrt(f_fac*16.04/m_mw[i])*pow(h_fac,-2./3.);
        double mu_i = mu_0*H*m_mw[i]*c1;
        L_i[i] = mu_i*1.32*GasConstant*(cp_0_R[i] - 2.5)/m_mw[i];
        L_i_min = min(L_i_min,L_i[i]);
        // Calculate variables for density-dependent component:
        double theta_s = 1 + (m_w_ac[i] - 0.011)*(0.09057 - 0.86276*log(T_p)
            + (0.31664 - 0.46568/T_p)*(V_p - 0.5));
        double phi_s = (1 + (m_w_ac[i] - 0.011)*(0.39490*(V_p - 1.02355)
            - 0.93281*(V_p - 0.75464)*log(T_p)))*0.288/Zcrit_i(i);
        f_i[i] = Tc_i*theta_s/190.4;
        h_i[i] = 1000*Vc_i*phi_s/99.2;
    }

    double h_m = 0;
    double f_m = 0;
    double mw_m = 0;
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Density-independent component:
            double L_ij = 2*L_i[i]*L_i[j]/(L_i[i] + L_i[j] + Tiny);
            Lprime_m += molefracs[i]*molefracs[j]*L_ij;
            // Additional variables for density-dependent component:
            double f_ij = sqrt(f_i[i]*f_i[j]);
            double h_ij = 0.125*pow(pow(h_i[i],1./3.) + pow(h_i[j],1./3.),3.);
            double mw_ij_inv = (m_mw[i] + m_mw[j])/(2*m_mw[i]*m_mw[j]);
            f_m += molefracs[i]*molefracs[j]*f_ij*h_ij;
            h_m += molefracs[i]*molefracs[j]*h_ij;
            mw_m += molefracs[i]*molefracs[j]*sqrt(mw_ij_inv*f_ij)*pow(h_ij,-4./3.);
        }
    }

    f_m = f_m/h_m;
    mw_m = pow(mw_m,-2.)*f_m*pow(h_m,-8./3.);

    double rho_0 = 16.04*h_m/(1000*m_thermo->molarVolume());
    double T_0 = m_temp/f_m;
    double mu_0 = 1e-7*(2.90774e6/T_0 - 3.31287e6*pow(T_0,-2./3.)
                + 1.60810e6*pow(T_0,-1./3.) - 4.33190e5 + 7.06248e4
                *pow(T_0,1./3.) - 7.11662e3*pow(T_0,2./3.) + 4.32517e2*T_0
                - 1.44591e1*pow(T_0,4./3.) + 2.03712e-1*pow(T_0,5./3.));
    double L_1m = 1944*mu_0;
    double L_2m = (-2.5276e-4 + 3.3433e-4*pow(1.12 - log(T_0/1.680e2),2))*rho_0;
    double L_3m = exp(-7.19771 + 85.67822/T_0)*(exp((12.47183
                - 984.6252*pow(T_0,-1.5))*pow(rho_0,0.1) + (rho_0/0.1617 - 1)
                *sqrt(rho_0)*(0.3594685 + 69.79841/T_0 - 872.8833*pow(T_0,-2))) - 1.)*1e-3;
    double H_m = sqrt(f_m*16.04/mw_m)*pow(h_m,-2./3.);
    double Lstar_m = H_m*(L_1m + L_2m + L_3m);
    return Lprime_m + Lstar_m;
}

void HighPressureGasTransport::getBinaryDiffCoeffs(const size_t ld, double* const d)
{
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    update_T();
    // If necessary, evaluate the binary diffusion coefficients from the polynomial fits
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    if (ld < m_nsp) {
        throw CanteraError("HighPressureGasTransport::getBinaryDiffCoeffs", "ld is too small");
    }

    double rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = 0; j < nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            // zero (this would lead to Pr_ij = Inf).
            double x_i = std::max(Tiny, molefracs[i]);
            double x_j = std::max(Tiny, molefracs[j]);

            // Weight mole fractions of i and j so that X_i + X_j = 1.0.
            x_i = x_i/(x_i + x_j);
            x_j = x_j/(x_i + x_j);

            // Calculate Tr and Pr based on mole-fraction-weighted critical constants.
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            // Calculate the parameters for Takahashi correlation
            double P_corr_ij;
            P_corr_ij = takahashi_correction_factor(Pr_ij, Tr_ij);

            // If the reduced temperature is too low, the correction factor
            // P_corr_ij will be < 0.
            if (P_corr_ij<0) {
                P_corr_ij = Tiny;
            }
            // Multiply the standard low-pressure binary diffusion coefficient
            // (m_bdiff) by the Takahashi correction factor P_corr_ij.
            d[ld*j + i] = P_corr_ij*(rp * m_bdiff(i,j));
        }
    }
}


// Calculate the high-pressure mixture viscosity, based on the Lucas method
// described in chapter 9-7 of Poling et al. (2001).
//
// The mixture pseudo-critical temperature and pressure are calculated using
// equation 9-5.18 and 9-5.19 in Poling et al. (2001).
//
// The mixture molecular weight is computed using equation 9-5.20 in Poling et al. (2001).
//
// The mixture values of the low-pressure polarity and quantum correction factors are
// computed using equations 9-5.21 and 9-5.22 in Poling et al. (2001).
double HighPressureGasTransport::viscosity()
{
    // Most of this function consists of computing mixture values of various critical
    // properties and other parameters that are used in the viscosity calculation.
    double Tc_mix = 0.0;
    double Pc_mix_n = 0.0; // Numerator in equation 9-5.18 in Poling et al. (2001)
    double Pc_mix_d = 0.0; // Denominator in equation 9-5.18 in Poling et al. (2001)
    double MW_mix = m_thermo->meanMolecularWeight();

    double FP_mix_o = 0; // The mole-fraction-weighted mixture average of the low-pressure polarity correction factor
    double FQ_mix_o = 0; // The mole-fraction-weighted mixture average of the low-pressure quantum correction factor
    double MW_H = m_mw[0]; // Holds the molecular weight of the heaviest species
    double MW_L = m_mw[0]; // Holds the molecular weight of the lightest species

    double tKelvin = m_thermo->temperature();
    double P_vap_mix = m_thermo->satPressure(tKelvin);
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    double x_H = molefracs[0]; // Holds the mole fraction of the heaviest species
    for (size_t i = 0; i < m_nsp; i++) {
        // Calculate pure-species critical constants and add their contribution
        // to the mole-fraction-weighted mixture averages:
        double Tc = Tcrit_i(i);
        double Tr = tKelvin/Tc;
        double Zc = Zcrit_i(i);
        Tc_mix += Tc*molefracs[i];
        Pc_mix_n += molefracs[i]*Zc; //numerator
        Pc_mix_d += molefracs[i]*Vcrit_i(i); //denominator (Units of Vcrit_i are m^3/kmol, which are fine because they cancel out with the Cantera gas constant's units used later)

        // Calculate ratio of heaviest to lightest species
        if (m_mw[i] > MW_H) {
            MW_H = m_mw[i];
            x_H = molefracs[i];
        } else if (m_mw[i] < MW_L) {
            MW_L = m_mw[i];
        }

        // Calculate pure-species reduced dipole moment for pure-species
        // polar correction term.
        // Equation 9-4.17 in Poling et al. (2001) requires the pressure
        // to be in units of bar, so we convert from Pa to bar.
        // The dipole moment is stored in SI units, and it needs to be in
        // units of Debye for the Lucas method.
        double pascals_to_bar = 1e-5;
        double SI_to_Debye = 1.0 / 3.335e-30; // Conversion factor from C*m to Debye
        double dipole_ii = m_dipole(i,i)*SI_to_Debye;
        double mu_ri = 52.46*dipole_ii*dipole_ii*(Pcrit_i(i)*pascals_to_bar)/(Tc*Tc);
        FP_mix_o += molefracs[i] * FP_i(mu_ri, Tr, Zc); // mole-fraction weighting of pure-species polar correction term


        // Calculate contribution to quantum correction term.
        // Note:  This assumes the species of interest (He, H2, and D2) have
        //        been named in this specific way.
        vector<string> spnames = m_thermo->speciesNames();
        if (spnames[i] == "He") {
            FQ_mix_o += molefracs[i]*FQ_i(1.38, Tr, m_mw[i]);
        } else if (spnames[i] == "H2") {
            FQ_mix_o += molefracs[i]*(FQ_i(0.76, Tr, m_mw[i]));
        } else if (spnames[i] == "D2") {
            FQ_mix_o += molefracs[i]*(FQ_i(0.52, Tr, m_mw[i]));
        } else {
            FQ_mix_o += molefracs[i];
        }
    }

    double Tr_mix = tKelvin/Tc_mix;
    double Pc_mix = GasConstant*Tc_mix*Pc_mix_n/Pc_mix_d;
    double Pr_mix = m_thermo->pressure()/Pc_mix;

    // Compute the mixture value of the low-pressure quantum correction factor
    double ratio = MW_H/MW_L;
    double A = 1.0;
    if (ratio > 9 && x_H > 0.05 && x_H < 0.7) {
        A = 1 - 0.01*pow(ratio,0.87);
    }
    FQ_mix_o *= A;


    // This is η*ξ
    double nondimensional_viscosity = high_pressure_nondimensional_viscosity(Tr_mix, Pr_mix, FP_mix_o, FQ_mix_o, P_vap_mix, Pc_mix);

    // Using equation 9-4.14 in Poling et al. (2001), with units of 1/(Pa*s)
    double numerator = GasConstant*Tc_mix*pow(Avogadro,2.0);
    double denominator = pow(MW_mix,3.0)*pow(Pc_mix,4.0);
    double ksi = pow(numerator / denominator, 1.0/6.0);

    // Return the viscosity in kg/m/s
    return nondimensional_viscosity / ksi;
}

// Pure species critical properties - Tc, Pc, Vc, Zc:
double HighPressureGasTransport::Tcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double tc = m_thermo->critTemperature();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return tc;
}

double HighPressureGasTransport::Pcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double pc = m_thermo->critPressure();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return pc;
}

double HighPressureGasTransport::Vcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double vc = m_thermo->critVolume();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return vc;
}

double HighPressureGasTransport::Zcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double zc = m_thermo->critCompressibility();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return zc;
}

// The low-pressure nondimensional viscosity equation 9-4.16 in Poling et al. (2001).
// This relation is used for pure species and mixtures at low pressure. The only
// difference is the the values that are passed to the function (pure values versus mixture values).
double HighPressureGasTransport::low_pressure_nondimensional_viscosity(double Tr, double FP, double FQ) {
    double first_term = 0.807*pow(Tr,0.618) - 0.357*exp(-0.449*Tr);
    double second_term = 0.340*exp(-4.058*Tr) + 0.018;
    return (first_term + second_term)*FP*FQ;
}

// The high-pressure nondimensional viscosity equation 9-6.12 in Poling et al. (2001).
// This relation is used for pure species and mixtures at high pressure. The only
// difference is the the values that are passed to the function (pure values versus mixture values).
// This returns the value of η*ξ (by multiplying both sides of 9-6.12 by ξ and simply returning the RHS of the
// resulting equation)
double HighPressureGasTransport::high_pressure_nondimensional_viscosity(double Tr, double Pr, double FP_low, double FQ_low, double P_vap, double P_crit){

    //cout << "Input values are:" << endl;
    //cout << "Tr: " << Tr << endl;
    //cout << "Pr: " << Pr << endl;
    //cout << "FP_low: " << FP_low << endl;
    //cout << "FQ_low: " << FQ_low << endl;
    //cout << "P_vap: " << P_vap << endl;
    //cout << "P_crit: " << P_crit << endl;

    double Z_1 = low_pressure_nondimensional_viscosity(Tr, FP_low, FQ_low); // This is η_0*ξ

    double Z_2;
    if (Tr <= 1.0) {
        if (Pr < P_vap/P_crit) {
            double alpha = 3.262 + 14.98*pow(Pr, 5.508);
            double beta = 1.390 + 5.746*Pr;
            Z_2 = 0.600 + 0.760*pow(Pr, alpha) + (0.6990*pow(Pr, beta) - 0.60) * (1-Tr);
        } else {
            throw CanteraError("HighPressureGasTransport::viscosity",
                               "State is outside the limits of the Lucas model, Tr <= 1");
        }
    } else if ((Tr > 1.0) && (Tr < 40.0)) {
        if ((Pr > 0.0) && (Pr <= 100.0)) {
            // The following expressions are given in page 9.36 of Poling et al. (2001)
            // and correspond to parameters in equation 9-6.8.
            double a_1 = 1.245e-3;
            double a_2 = 5.1726;
            double gamma = -0.3286;
            double a = a_1*exp(a_2*pow(Tr,gamma))/Tr;

            double b_1 = 1.6553;
            double b_2 = 1.2723;
            double b = a*(b_1*Tr - b_2);

            double c_1 = 0.4489;
            double c_2 = 3.0578;
            double delta = -37.7332;
            double c = c_1*exp(c_2*pow(Tr, delta))/Tr;

            double d_1 = 1.7368;
            double d_2 = 2.2310;
            double epsilon = -7.6351;
            double d = d_1*exp(d_2*pow(Tr, epsilon))/Tr;

            double e = 1.3088;

            double f_1 = 0.9425;
            double f_2 = -0.1853;
            double zeta = 0.4489;
            double f = f_1*exp(f_2*pow(Tr, zeta));

            Z_2 = Z_1*(1 + (a*pow(Pr,e)) / (b*pow(Pr,f) + pow(1+c*pow(Pr,d),-1)));
        } else {
            throw CanteraError("HighPressureGasTransport::viscosity",
                           "State is outside the limits of the Lucas model, 1.0 < Tr < 40");
        }
    } else {
        throw CanteraError("HighPressureGasTransport::viscosity",
                           "State is outside the limits of the Lucas model, Tr > 40");
    }

    double Y = Z_2 / Z_1;
    double FP = (1 + (FP_low - 1)*pow(Y,-3.0)) / FP_low;
    double FQ = (1 + (FQ_low - 1)*(1.0/Y - 0.007*pow(log(Y),4.0))) / FQ_low;

    // Return the non-dimensional viscosity η*ξ
    return Z_2 * FP * FQ;
}

// Calculates quantum correction term of the Lucas method for a species based
// on Tr and MW, used in viscosity calculation from equation 9-4.19 in
// Poling et al. (2001)
double HighPressureGasTransport::FQ_i(double Q, double Tr, double MW)
{
    return 1.22*pow(Q,0.15)*(1 + 0.00385*pow(pow(Tr - 12.0, 2.0), 1.0/MW)
                             *sign(Tr - 12.0));
}

// Calculates the pressure correction factor for the Lucas method that
// is used to calculate high pressure pure fluid viscosity. This is equation
// 9-4.18 in Poling et al. (2001).
//
// NOTE: The original description in the book neglects to mention what happens
// when the quantity raised to the 1.72 power goes negative. That is an undefined
// operation that generates real+imaginary numbers. For now, I simply am
// taking the absolute value of the argument and raising it to the power.
double HighPressureGasTransport::FP_i(double mu_r, double Tr, double Z_crit)
{
    //cout << "HighPressureGasTransport::FP_i Inputs:" << endl;
    //cout << "mu_r: " << mu_r << endl;
    //cout << "Tr: " << Tr << endl;
    //cout << "Z_crit: " << Z_crit << endl;

    if (mu_r < 0.022) {
        return 1;
    } else if (mu_r < 0.075) {
        return 1 + 30.55*pow(fabs(0.292 - Z_crit), 1.72);
    } else {
        return 1 + 30.55*pow(fabs(0.292 - Z_crit), 1.72)*fabs(0.96 + 0.1*(Tr - 0.7));
    }
}

}

//Chung Implementation
namespace Cantera
{

void ChungHighPressureGasTransport::getBinaryDiffCoeffs(const size_t ld, double* const d)
{
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    update_T();
    // If necessary, evaluate the binary diffusion coefficients from the polynomial fits
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    if (ld < m_nsp) {
        throw CanteraError("ChungHighPressureGasTransport::getBinaryDiffCoeffs", "ld is too small");
    }

    double rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < nsp; i++) {
        for (size_t j = 0; j < nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            // zero (this would lead to Pr_ij = Inf).
            double x_i = std::max(Tiny, molefracs[i]);
            double x_j = std::max(Tiny, molefracs[j]);

            // Weight mole fractions of i and j so that X_i + X_j = 1.0.
            x_i = x_i/(x_i + x_j);
            x_j = x_j/(x_i + x_j);

            // Calculate Tr and Pr based on mole-fraction-weighted critical constants.
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            // Calculate the parameters for Takahashi correlation
            double P_corr_ij;
            P_corr_ij = takahashi_correction_factor(Pr_ij, Tr_ij);

            // If the reduced temperature is too low, the correction factor
            // P_corr_ij will be < 0.
            if (P_corr_ij<0) {
                P_corr_ij = Tiny;
            }
            // Multiply the standard low-pressure binary diffusion coefficient
            // (m_bdiff) by the Takahashi correction factor P_corr_ij.
            d[ld*j + i] = P_corr_ij*(rp * m_bdiff(i,j));
        }
    }
}

// Calculate the high-pressure mixture thermal conductivity, based on the Chung method
// described in on page 10.23 of Poling et al. (2001).
//
// The mixture pseudo-critical temperature and pressure are calculated using
// the equations defined on page 9.25 of Poling et al. (2001).
//
// The mixture value of the specific heat is computed using equation 10-6.6.
double ChungHighPressureGasTransport::thermalConductivity()
{
    ChungMixtureParameters params;
    compute_mixture_parameters(params);

    // Compute T_star using equation 9-5.26, using the mixture parameters
    double tKelvin = m_thermo->temperature();
    double T_star = tKelvin / params.epsilon_over_k_mix;

    // The density is required for high-pressure gases.
    // The Chung method requires density to be units of mol/cm^3
    // We use the mixture molecular weight (units of kg/kmol) here.
    double kg_per_m3_to_mol_per_cm3 = (1.0 / params.MW_mix)*1e-3; // 1 kmol/m^3 = 1e-3 mol/cm^3
    double density = m_thermo->density()*kg_per_m3_to_mol_per_cm3;

    // The value of Cv is already a mole-weighted average of the pure species values
    double Cv_mix = m_thermo->cv_mole(); // Specific heat at constant volume, units are J/kmol/K

    // This result is in units of W/m/K
    double thermal_conductivity = high_pressure_thermal_conductivity(tKelvin, T_star, params.MW_mix, density, Cv_mix, params.Vc_mix, params.Tc_mix, params.sigma_mix, params.acentric_factor_mix, params.mu_r_mix, params.kappa_mix);

    // Return the thermal conductivity in W/m/K
    return thermal_conductivity;
}



// Computes the high-pressure thermal conductivity using the Chung method (Equation 10-5.5).
// This function is structured such that it can be used for pure species or mixtures, with the
// only difference being the values that are passed to the function (pure values versus mixture values).
//
// This method utilizes the low-pressure Chung viscosity as that is a required parameter in the model, and
// thus makes a call to the low pressure viscosity implementation.
//
// M_prime (M' in the model) has units of kg/mol, and is just the molecular weight (kg/kmol) divided by 1000.
double ChungHighPressureGasTransport::high_pressure_thermal_conductivity(double T, double T_star, double MW, double rho, double Cv, double Vc, double Tc, double sigma, double acentric_factor, double mu_r, double kappa)
{
    // Calculate the low-pressure viscosity using the Chung method
    // This method returns viscosity in micropoise, but the thermal
    // conductivity model needs the low-pressure viscosity to be in units of Pa*s
    double micropoise_to_pascals_second = 1e-7;
    double viscosity = low_pressure_viscosity(T, T_star, MW, acentric_factor, mu_r, sigma, kappa)*micropoise_to_pascals_second;
    //cout << "Low-pressure viscosity: " << viscosity << " Pa*s" << endl;

    double M_prime = MW / 1000.0; // Converting to kg/mol

    // Definition of tabulated coefficients for the Chung method, as shown in Table 10-3 on page 10.23 of Poling et al. (2001)
    vector<double> a = {2.44166, -5.0924e-1, 6.6107, 1.4543e1, 7.9274e-1, -5.8634, 9.1089e1};
    vector<double> b = {7.4824e-1, -1.5094, 5.6207, -8.9139, 8.2019e-1, 1.2801e1, 1.2811e2};
    vector<double> c = {-9.1858e-1, -4.9991e1, 6.4760e1, -5.6379, -6.9369e-1, 9.5893, -5.4217e1};
    vector<double> d ={1.2172e2, 6.9983e1, 2.7039e1, 7.4344e1, 6.3173, 6.5529e1, 5.2381e2};

    // This is slightly pedantic, but this is done to have the naming convention in the
    // equations used match the variable names in the code.
    double B_1, B_2, B_3, B_4, B_5, B_6, B_7;
    B_1 = a[0] + b[0]*acentric_factor + c[0]*pow(mu_r, 4.0) + d[0]*kappa;
    B_2 = a[1] + b[1]*acentric_factor + c[1]*pow(mu_r, 4.0) + d[1]*kappa;
    B_3 = a[2] + b[2]*acentric_factor + c[2]*pow(mu_r, 4.0) + d[2]*kappa;
    B_4 = a[3] + b[3]*acentric_factor + c[3]*pow(mu_r, 4.0) + d[3]*kappa;
    B_5 = a[4] + b[4]*acentric_factor + c[4]*pow(mu_r, 4.0) + d[4]*kappa;
    B_6 = a[5] + b[5]*acentric_factor + c[5]*pow(mu_r, 4.0) + d[5]*kappa;
    B_7 = a[6] + b[6]*acentric_factor + c[6]*pow(mu_r, 4.0) + d[6]*kappa;

    double y = rho*Vc/6.0; // Equation 10-5.6 (with rho = 1/V)

    double G_1 = (1.0 - 0.5*y)/(pow(1.0-y, 3.0)); // Equation  10-5.7

    // Equation 10-5.8
    double G_2 = (B_1*((1.0-exp(-B_4*y)) / y) + B_2*G_1*exp(B_5*y) + B_3*G_1)/ (B_1*B_4 + B_2 + B_3);

    double q = 3.586e-3*sqrt(Tc/M_prime) / pow(Vc, 2.0/3.0); // Shown below equation 10-5.5

    double Tr = T/Tc; // Reduced temperature
    double alpha = (Cv / GasConstant) - 3.0/2.0; // Shown below equation 10-3.14, using Cantera's R (J/kmol/K )
    double beta = 0.7862 - 0.7109*acentric_factor + 1.3168*acentric_factor*acentric_factor; // Shown below Equation 10-3.14
    double Z = 2.0 + 10.5*Tr*Tr; // Shown below Equation 10-3.14
    double psi = 1.0 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z) / (0.6366 + beta*Z + 1.061*alpha*beta)); // Shown below Equation 10-3.14

    //cout << "Parameters: " << endl;
    //cout << "y: " << y << endl;
    //cout << "G_1: " << G_1 << endl;
    //cout << "G_2: " << G_2 << endl;
    //cout << "q: " << q << endl;
    //cout << "Tr: " << Tr << endl;
    //cout << "alpha: " << alpha << endl;
    //cout << "beta: " << beta << endl;
    //cout << "Z: " << Z << endl;
    //cout << "psi: " << psi << endl;

    double lambda = (31.2*viscosity*psi/M_prime)*(1.0/G_2 + B_6*y) + q*B_7*y*y*sqrt(Tr)*G_2; // Equation 10-5.5

    // Units are W/m/K
    return lambda;
}




// This implements the high-pressure gas mixture viscosity model of Chung
// described in chapter 9-7 of Poling et al. (2001).
double ChungHighPressureGasTransport::viscosity()
{
    ChungMixtureParameters params;
    compute_mixture_parameters(params);

    // Compute T_star using equation 9-5.26, using the mixture parameters
    double tKelvin = m_thermo->temperature();
    double T_star = tKelvin / params.epsilon_over_k_mix;

    // The density is required for high-pressure gases.
    // The Chung method requires density to be units of mol/cm^3
    // We use the mixture molecular weight (units of kg/kmol) here.
    double kg_per_m3_to_mol_per_cm3 = (1.0 / params.MW_mix)*1e-3; // 1 kmol/m^3 = 1e-3 mol/cm^3
    double density = m_thermo->density()*kg_per_m3_to_mol_per_cm3;

    // This result is in units of micropoise
    double viscosity = high_pressure_viscosity(T_star, params.MW_mix, density, params.Vc_mix, params.Tc_mix, params.acentric_factor_mix, params.mu_r_mix, params.kappa_mix);

    double micropoise_to_pascals_second = 1e-7;
    return viscosity*micropoise_to_pascals_second;

}



// Implementation of equation the mixing rules defined on page 9.25 of Poling et al. (2001)
// for the Chung method's composition dependent parameters.
void ChungHighPressureGasTransport::compute_mixture_parameters(ChungMixtureParameters& params)
{
    size_t nsp = m_thermo->nSpecies();
    vector<double> molefracs(nsp);
    m_thermo->getMoleFractions(&molefracs[0]);

    // First fill the species-specific values that will then later be used in the combining
    // rules for the Chung method.
    vector<double> sigma_i(nsp), epsilon_over_k_i(nsp), acentric_factor_i(nsp), MW_i(nsp), kappa_i(nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        // From equation 9-5.32 in Poling et al. (2001)
        double m3_per_kmol_to_cm3_per_mol = 1e3; // Convert from m^3/kmol to cm^3/mol
        double Vc = Vcrit_i(i) * m3_per_kmol_to_cm3_per_mol;
        sigma_i[i] = 0.809*pow(Vc, 1.0/3.0);

        // From equation 9-5.34 in Poling et al. (2001)
        epsilon_over_k_i[i] = Tcrit_i(i)/1.2593;

        // NOTE: The association parameter is assumed to be zero for all species, but
        // is left here for completeness or future revision.
        kappa_i[i] = 0.0;

        // These values are available from the base class
        acentric_factor_i[i] = m_w_ac[i];
        MW_i[i] = m_mw[i];
    }

    // Here we use the combining rules defined on page 9.25 of Poling et al. (2001)
    // We have ASSUMED that the binary interaction parameters are unity for all species
    // as was done in the Chung method.
    //
    // The sigma & kappa relations can be fully computed in the loop. The other ones require a final
    // division by the final mixture values, and so the quantities in the loop are only the
    // numerators of the equations that are referenced in the comments. After the loop, the
    // final mixture values are computed and stored.
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j <m_nsp; j++){
            double sigma_ij = sqrt(sigma_i[i]*sigma_i[j]); // Equation 9-5.33
            params.sigma_mix += molefracs[i]*molefracs[j]*pow(sigma_ij,3.0); // Equation 9-5.25

            double epsilon_over_k_ij = sqrt(epsilon_over_k_i[i]*epsilon_over_k_i[j]); // Equation 9-5.35
            params.epsilon_over_k_mix += molefracs[i]*molefracs[j]*epsilon_over_k_ij*pow(sigma_ij,3.0); // Equation 9-5.27 (numerator only)

            // This equation is raised to the power 2, so what we do here is first store only the
            // double summation into this variable, and then later square the result.
            double MW_ij =  (2*MW_i[i]*MW_i[j]) / (MW_i[i] + MW_i[j]); // Equation 9-5.40
            params.MW_mix += molefracs[i]*molefracs[j]*epsilon_over_k_ij*pow(sigma_ij,2.0)*sqrt(MW_ij); // Equation 9-5.28 (double summation only)

            double acentric_factor_ij = 0.5*(acentric_factor_i[i] + acentric_factor_i[j]); // Equation 9-5.36
            params.acentric_factor_mix += molefracs[i]*molefracs[j]*acentric_factor_ij*pow(sigma_ij,3.0); // Equation 9-5.29

            // Using the base class' dipole moment values. These are in the SI units, so we need to convert to
            // Debye units.
            double SI_to_Debye = 1.0 / 3.335e-30; // Conversion factor from C*m to Debye
            double dipole_ii = m_dipole(i,i)*SI_to_Debye;
            double dipole_jj = m_dipole(j,j)*SI_to_Debye;
            params.mu_mix += molefracs[i]*molefracs[j]*pow(dipole_ii,2.0)*pow(dipole_jj,2.0)/pow(sigma_ij,3.0); // Equation 9-5.30

            // Using equation 9-5.31
            double kappa_ij = sqrt(kappa_i[i]*kappa_i[j]); // Equation 9-5.38
            params.kappa_mix += molefracs[i]*molefracs[j]*kappa_ij; // Equation 9-5.31
        }
    }

    // Finalize the expressions for the mixture values of the parameters
    params.sigma_mix = pow(params.sigma_mix, 1.0/3.0); // Equation 9-5.25 computed the cube of sigma_mix
    params.epsilon_over_k_mix /= pow(params.sigma_mix,3.0);

    // The MW_mix was only the numerator inside the brackets of equation 9-5.28
    params.MW_mix = pow(params.MW_mix/(params.epsilon_over_k_mix*pow(params.sigma_mix,2.0)), 2.0);

    params.acentric_factor_mix /= pow(params.sigma_mix,3.0);

    params.mu_mix = pow(pow(params.sigma_mix,3.0)*params.mu_mix, 1.0/4.0); // Equation 9-5.30 computed the 4th power of mu_mix

    // Tc_mix is computed using equation 9-5.44 in Poling et al. (2001)
    params.Tc_mix = 1.2593*params.epsilon_over_k_mix;

    // Vc_mix is computed using equation 9-5.43 in Poling et al. (2001)
    params.Vc_mix = pow(params.sigma_mix/0.809, 3.0);

    // mu_r_mix is computed using equation 9-5.42 in Poling et al. (2001)
    params.mu_r_mix = 131.3*params.mu_mix/sqrt(params.Vc_mix*params.Tc_mix);

    //cout << "Mixture parameters are:" << endl;
    //cout << "sigma_mix: " << params.sigma_mix << endl;
    //cout << "epsilon_over_k_mix: " << params.epsilon_over_k_mix << endl;
    //cout << "MW_mix: " << params.MW_mix << endl;
    //cout << "acentric_factor_mix: " << params.acentric_factor_mix << endl;
    //cout << "mu_mix: " << params.mu_mix << endl;
    //cout << "Tc_mix: " << params.Tc_mix << endl;
    //cout << "Vc_mix: " << params.Vc_mix << endl;
    //cout << "mu_r_mix: " << params.mu_r_mix << endl;
    //cout << "kappa_mix: " << params.kappa_mix << endl;

}

// Implementation of equation 9-4.3 in Poling et al. (2001).
// Applicable over the range of 0.3 <= T_star <= 100.
double neufeld_collision_integral(double T_star)
{
    double A = 1.16145;
    double B = 0.14874;
    double C = 0.52487;
    double D = 0.77320;
    double E = 2.16178;
    double F = 2.43787;

    double omega = A / pow(T_star, B) + C / exp(D*T_star) + E / exp(F*T_star);
    return omega;
}

// This function is structured such that it can be used for pure species or mixtures, with the
// only difference being the values that are passed to the function (pure values versus mixture values).
double ChungHighPressureGasTransport::low_pressure_viscosity(double T, double T_star, double MW, double acentric_factor, double mu_r, double sigma, double kappa)
{
    //cout << "Input values are:" << endl;
    //cout << "T: " << T << endl;
    //cout << "T_star: " << T_star << endl;
    //cout << "MW: " << MW << endl;
    //cout << "acentric_factor: " << acentric_factor << endl;
    //cout << "mu_r: " << mu_r << endl;
    //cout << "sigma: " << sigma << endl;
    //cout << "kappa: " << kappa << endl;

    // Equation 9-4.3 in Poling et al. (2001)
    double omega = neufeld_collision_integral(T_star);

    // Molecular shapes and polarities factor, equation 9-4.11
    double Fc = 1 -0.2756*acentric_factor + 0.059035*pow(mu_r, 4.0) + kappa;
    //cout << "Fc: " << Fc << ", Omega: " << omega << endl;

    // Equation 9-3.9 in Poling et al. (2001), multiplied by the Chung factor, Fc
    // (another way of writing 9-4.10 that avoids explicit use of the critical volume in this method)
    double viscosity = Fc* (26.69*sqrt(MW*T)/(sigma*sigma*omega));

    return viscosity;
}

// Computes the high-pressure viscosity using the Chung method (Equation 9-6.18).
// This function is structured such that it can be used for pure species or mixtures, with the
// only difference being the values that are passed to the function (pure values versus mixture values).
//
// Renamed eta_star and eta_star_star from the Poling description to eta_1 and eta_2 for
// naming simplicity.
double ChungHighPressureGasTransport::high_pressure_viscosity(double T_star, double MW, double rho, double Vc, double Tc, double acentric_factor, double mu_r, double kappa)
{
    //cout << "Input values are:" << endl;
    //cout << "T_star: " << T_star << endl;
    //cout << "MW: " << MW << endl;
    //cout << "rho: " << rho << endl;
    //cout << "Vc: " << Vc << endl;
    //cout << "Tc: " << Tc << endl;
    //cout << "acentric_factor: " << acentric_factor << endl;
    //cout << "mu_r: " << mu_r << endl;
    //cout << "kappa: " << kappa << endl;

    // Definition of tabulated coefficients for the Chung method, as shown in Table 9-6 on page 9.40 of Poling et al. (2001)
    vector<double> a = {6.324, 1.210e-3, 5.283, 6.623, 19.745, -1.900, 24.275, 0.7972, -0.2382, 0.06863};
    vector<double> b = {50.412, -1.154e-3, 254.209, 38.096, 7.630, -12.537, 3.450, 1.117, 0.06770, 0.3479};
    vector<double> c = {-51.680, -6.257e-3, -168.48, -8.464, -14.354, 4.958, -11.291, 0.01235, -0.8163, 0.5926};
    vector<double> d ={1189.0, 0.03728, 3898.0, 31.42, 31.53, -18.15, 69.35, -4.117, 4.025, -0.727};

    // This is slightly pedantic, but this is done to have the naming convention in the
    // equations used match the variable names in the code.
    double E_1, E_2, E_3, E_4, E_5, E_6, E_7, E_8, E_9, E_10;
    E_1 = a[0] + b[0]*acentric_factor + c[0]*pow(mu_r, 4.0) + d[0]*kappa;
    E_2 = a[1] + b[1]*acentric_factor + c[1]*pow(mu_r, 4.0) + d[1]*kappa;
    E_3 = a[2] + b[2]*acentric_factor + c[2]*pow(mu_r, 4.0) + d[2]*kappa;
    E_4 = a[3] + b[3]*acentric_factor + c[3]*pow(mu_r, 4.0) + d[3]*kappa;
    E_5 = a[4] + b[4]*acentric_factor + c[4]*pow(mu_r, 4.0) + d[4]*kappa;
    E_6 = a[5] + b[5]*acentric_factor + c[5]*pow(mu_r, 4.0) + d[5]*kappa;
    E_7 = a[6] + b[6]*acentric_factor + c[6]*pow(mu_r, 4.0) + d[6]*kappa;
    E_8 = a[7] + b[7]*acentric_factor + c[7]*pow(mu_r, 4.0) + d[7]*kappa;
    E_9 = a[8] + b[8]*acentric_factor + c[8]*pow(mu_r, 4.0) + d[8]*kappa;
    E_10 = a[9] + b[9]*acentric_factor + c[9]*pow(mu_r, 4.0) + d[9]*kappa;

    double y = rho*Vc/6.0; // Equation 9-6.20

    double G_1 = (1.0 - 0.5*y)/(pow(1.0-y, 3.0)); // Equation 9-6.21

    // Equation 9-6.22
    double G_2 = (E_1*((1.0-exp(-E_4*y)) / y) + E_2*G_1*exp(E_5*y) + E_3*G_1)/ (E_1*E_4 + E_2 + E_3);

    double eta_2 = E_7*y*y*G_2*exp(E_8 + E_9/T_star + E_10/(T_star*T_star)); // Equation 9-6.23

    // Equation 9-4.3 in Poling et al. (2001)
    double omega = neufeld_collision_integral(T_star);

    // Molecular shapes and polarities factor, equation 9-4.11
    double Fc = 1 -0.2756*acentric_factor + 0.059035*pow(mu_r, 4.0) + kappa;

    double eta_1 = (sqrt(T_star)/omega) * (Fc*(1.0/G_2 + E_6*y)) + eta_2; // Equation 9-6.19

    double eta = eta_1 * 36.344 * sqrt(MW*Tc) / pow(Vc, 2.0/3.0); // Equation 9-6.18

    return eta;
}


// Pure species critical properties - Tc, Pc, Vc, Zc:
double ChungHighPressureGasTransport::Tcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double tc = m_thermo->critTemperature();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return tc;
}

double ChungHighPressureGasTransport::Pcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double pc = m_thermo->critPressure();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return pc;
}

double ChungHighPressureGasTransport::Vcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double vc = m_thermo->critVolume();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return vc;
}

double ChungHighPressureGasTransport::Zcrit_i(size_t i)
{
    vector<double> molefracs(m_thermo->nSpecies());
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> mf_temp(m_thermo->nSpecies(), 0.0);
    mf_temp[i] = 1;
    m_thermo->setMoleFractions(&mf_temp[0]);

    double zc = m_thermo->critCompressibility();

    // Restore actual molefracs:
    m_thermo->setMoleFractions(&molefracs[0]);
    return zc;
}


}
