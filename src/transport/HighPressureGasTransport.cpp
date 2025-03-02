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
#include "cantera/thermo/Species.h"
#include <boost/algorithm/string.hpp>

namespace Cantera
{

/**
 * @brief Returns interpolated value of (DP)_R obtained from the data in Table 2 of
 * the Takahashi 1975 paper, given a value of the reduced pressure (Pr) and reduced
 * temperature (Tr).
 *
 * @param Pr  Reduced pressure
 * @param Tr  Reduced temperature
 */
double takahashiCorrectionFactor(double Pr, double Tr)
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

    // Interpolate to obtain the value of (DP)_R at
    // the provided value of the reduced pressure (Pr).
    int Pr_lower = 0; // Index of the lower bounding value of Pr
    int Pr_upper = 0; // Index of the upper bounding value of Pr
    double frac = 0.0;

    bool found = false;
    for (int j = 1; j < 17; j++){
        if (Pr_lookup[j] > Pr) {
            frac = (Pr - Pr_lookup[j-1])/(Pr_lookup[j] - Pr_lookup[j-1]);
            found = true;
            Pr_lower = j-1;
            Pr_upper = j;
            break;
        }
    }
    // If this loop completes without finding a bounding value of Pr, use
    // the final table value.
    if (!found) {
        Pr_lower = 16;
        Pr_upper = 16;
        frac = 1.0;
    }

    // Compute the value of (DP)_R at the given Pr value by interpolating the
    // bounding values of (DP)_R.
    double A, B, C, E, DP_Rt, DP_R_lower, DP_R_upper;
    DP_Rt = DP_Rt_lookup[Pr_lower];
    A = A_ij_lookup[Pr_lower];
    B = B_ij_lookup[Pr_lower];
    C = C_ij_lookup[Pr_lower];
    E = E_ij_lookup[Pr_lower];

    DP_R_lower = DP_Rt*(1.0 - A*pow(Tr,-B))*(1.0 - C*pow(Tr,-E));

    DP_Rt = DP_Rt_lookup[Pr_upper];
    A = A_ij_lookup[Pr_upper];
    B = B_ij_lookup[Pr_upper];
    C = C_ij_lookup[Pr_upper];
    E = E_ij_lookup[Pr_upper];

    DP_R_upper = DP_Rt*(1.0 - A*pow(Tr,-B))*(1.0 - C*pow(Tr,-E));

    // Linear interpolation of the two bounding values of (DP)_R.
    return DP_R_lower*(1.0 - frac) + DP_R_upper*frac;
}

void HighPressureGasTransportBase::getTransportData()
{
    // Call the base class's method to fill the properties
    GasTransport::getTransportData();

    // Contents of 'critical-properties.yaml', loaded later if needed
    AnyMap critPropsDb;
    std::unordered_map<string, AnyMap*> dbSpecies;

    // If a species has a zero acentric factor, check the critical-properties.yaml
    // database to see if it has a value specified for the species.
    for (size_t k = 0; k < m_thermo->nSpecies(); k++) {
        if (m_w_ac[k] == 0.0) {
            // Load 'crit-properties.yaml' file if not already loaded
            if (critPropsDb.empty()) {
                critPropsDb = AnyMap::fromYamlFile("critical-properties.yaml");
                dbSpecies = critPropsDb["species"].asMap("name");
            }

            // All names in critical-properties.yaml are upper case
            auto ucName = boost::algorithm::to_upper_copy(m_thermo->species(k)->name);
            if (dbSpecies.count(ucName)) {
                auto& spec = *dbSpecies.at(ucName);
                auto& critProps = spec["critical-parameters"].as<AnyMap>();
                if (critProps.hasKey("acentric-factor")) {
                    m_w_ac[k] = critProps.convert("acentric-factor", "1");
                }
            }
        }
    }
}

void HighPressureGasTransportBase::initializeCriticalProperties()
{
    size_t nSpecies = m_thermo->nSpecies();
    m_Tcrit.resize(nSpecies);
    m_Pcrit.resize(nSpecies);
    m_Vcrit.resize(nSpecies);
    m_Zcrit.resize(nSpecies);

    std::vector<double> molefracs(nSpecies);
    m_thermo->getMoleFractions(&molefracs[0]);

    std::vector<double> mf_temp(nSpecies, 0.0);

    for (size_t i = 0; i < nSpecies; ++i) {
        mf_temp[i] = 1.0;
        m_thermo->setMoleFractions(&mf_temp[0]);

        if (m_thermo->critTemperature() > 1e4) {
            throw CanteraError(
                "HighPressureGasTransport::initializeCriticalProperties",
                "Species '{}' must have critical properties defined or non-zero "
                "cubic parameters. Check the species definition or the thermo "
                "data file.",
                m_thermo->species(i)->name);
        }
        m_Tcrit[i] = m_thermo->critTemperature();
        m_Pcrit[i] = m_thermo->critPressure();
        m_Vcrit[i] = m_thermo->critVolume();
        m_Zcrit[i] = m_thermo->critCompressibility();

        mf_temp[i] = 0.0;  // Reset for the next iteration
    }

    // Restore actual mole fractions
    m_thermo->setMoleFractions(&molefracs[0]);
}

// Pure species critical properties - Tc, Pc, Vc, Zc:
double HighPressureGasTransportBase::Tcrit_i(size_t i)
{
    return m_Tcrit[i];
}

double HighPressureGasTransportBase::Pcrit_i(size_t i)
{
    return m_Pcrit[i];
}

double HighPressureGasTransportBase::Vcrit_i(size_t i)
{
    return m_Vcrit[i];
}

double HighPressureGasTransportBase::Zcrit_i(size_t i)
{
    return m_Zcrit[i];
}

void HighPressureGasTransportBase::updateCorrectionFactors() {
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Add an offset to avoid a condition where x_i and x_j both equal
            // zero (this would lead to Pr_ij = Inf).
            double x_i = std::max(Tiny, m_molefracs[i]);
            double x_j = std::max(Tiny, m_molefracs[j]);

            // Weight mole fractions of i and j so that X_i + X_j = 1.0.
            double sum_x_ij = x_i + x_j;
            x_i = x_i/(sum_x_ij);
            x_j = x_j/(sum_x_ij);

            // Calculate Tr and Pr based on mole-fraction-weighted critical constants.
            double Tr_ij = m_temp/(x_i*Tcrit_i(i) + x_j*Tcrit_i(j));
            double Pr_ij = m_thermo->pressure()/(x_i*Pcrit_i(i) + x_j*Pcrit_i(j));

            // Calculate the parameters for Takahashi correlation
            double P_corr_ij;
            P_corr_ij = takahashiCorrectionFactor(Pr_ij, Tr_ij);

            // If the reduced temperature is too low, the correction factor
            // P_corr_ij will be < 0.
            if (P_corr_ij<0) {
                P_corr_ij = Tiny;
            }
            m_P_corr_ij(i, j) = P_corr_ij;
        }
    }
}

void HighPressureGasTransportBase::getBinaryDiffCoeffs(const size_t ld, double* const d)
{
    update_C();
    update_T();
    updateCorrectionFactors();
    // If necessary, evaluate the binary diffusion coefficients from the polynomial
    // fits
    if (!m_bindiff_ok) {
        updateDiff_T();
    }
    if (ld < m_nsp) {
        throw CanteraError("HighPressureGasTransport::getBinaryDiffCoeffs",
                           "ld is too small");
    }

    double rp = 1.0/m_thermo->pressure();
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Multiply the standard low-pressure binary diffusion coefficient
            // (m_bdiff) by the Takahashi correction factor P_corr_ij.
            d[ld*j + i] = m_P_corr_ij(i,j)*(rp * m_bdiff(i,j));
        }
    }
}

void HighPressureGasTransportBase::getMixDiffCoeffs(double* const d)
{
    update_T();
    update_C();
    updateCorrectionFactors();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    double mmw = m_thermo->meanMolecularWeight();
    double p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_P_corr_ij(0,0)*m_bdiff(0,0) / p;
    } else {
        for (size_t i = 0; i < m_nsp; i++) {
            double sum2 = 0.0;
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    sum2 += m_molefracs[j] / (m_P_corr_ij(i,j)*m_bdiff(j,i));
                }
            }
            if (sum2 <= 0.0) {
                d[i] = m_P_corr_ij(i,i)*m_bdiff(i,i) / p;
            } else {
                d[i] = (mmw - m_molefracs[i] * m_mw[i])/(p * mmw * sum2);
            }
        }
    }
}

void HighPressureGasTransportBase::getMixDiffCoeffsMole(double* const d)
{
    update_T();
    update_C();
    updateCorrectionFactors();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    double p = m_thermo->pressure();
    if (m_nsp == 1) {
        d[0] = m_P_corr_ij(0,0)*m_bdiff(0,0) / p;
    } else {
        for (size_t i = 0; i < m_nsp; i++) {
            double sum2 = 0.0;
            for (size_t j = 0; j < m_nsp; j++) {
                if (j != i) {
                    sum2 += m_molefracs[j] / (m_P_corr_ij(i,j)*m_bdiff(j,i));
                }
            }
            if (sum2 <= 0.0) {
                d[i] = m_P_corr_ij(i,i)*m_bdiff(i,i) / p;
            } else {
                d[i] = (1 - m_molefracs[i]) / (p * sum2);
            }
        }
    }
}

void HighPressureGasTransportBase::getMixDiffCoeffsMass(double* const d)
{
    update_T();
    update_C();
    updateCorrectionFactors();

    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    double mmw = m_thermo->meanMolecularWeight();
    double p = m_thermo->pressure();

    if (m_nsp == 1) {
        d[0] = m_P_corr_ij(0,0)*m_bdiff(0,0) / p;
    } else {
        for (size_t i=0; i<m_nsp; i++) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (size_t j=0; j<m_nsp; j++) {
                if (j==i) {
                    continue;
                }
                sum1 += m_molefracs[j] / (m_P_corr_ij(i,j)*m_bdiff(i,j));
                sum2 += m_molefracs[j] * m_mw[j] / (m_P_corr_ij(i,j)*m_bdiff(i,j));
            }
            sum1 *= p;
            sum2 *= p * m_molefracs[i] / (mmw - m_mw[i]*m_molefracs[i]);
            d[i] = 1.0 / (sum1 + sum2);
        }
    }
}

// HighPressureGasTransport Implementation
// ---------------------------------------
void HighPressureGasTransport::init(ThermoPhase* thermo, int mode)
{
    MixTransport::init(thermo, mode);
    initializeCriticalProperties();
    m_P_corr_ij.resize(m_nsp, m_nsp);
}

double HighPressureGasTransport::thermalConductivity()
{
    update_T();
    vector<double> molefracs(m_nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    vector<double> cp_0_R(m_nsp);
    m_thermo->getCp_R_ref(&cp_0_R[0]); // Cp/R

    // A model constant from the Euken correlation for polyatomic gases, described
    // below Equation 1 in ely-hanley1981 .
    const double f_int = 1.32;

    // Pure-species model parameters
    // Internal contribution to thermal conductivity (lamba'')
    vector<double> Lambda_1_i(m_nsp);
    vector<double> f_i(m_nsp);
    vector<double> h_i(m_nsp);
    vector<double> V_k(m_nsp);

    m_thermo -> getPartialMolarVolumes(&V_k[0]);
    for (size_t i = 0; i < m_nsp; i++) {
        // Calculate variables for density-independent component, Equation 1,
        // the equation requires the pure-species viscosity estimate from
        // Ely and Hanley.
        double mu_i = elyHanleyDilutePureSpeciesViscosity(V_k[i], Tcrit_i(i),
                                                          Vcrit_i(i), Zcrit_i(i),
                                                          m_w_ac[i], m_mw[i]);

        // This is the internal contribution to the thermal conductivity of
        // pure-species component, i, from Equation 1 in ely-hanley1983
        Lambda_1_i[i] = (mu_i / m_mw[i])*f_int*GasConstant*(cp_0_R[i] - 2.5);

        // Calculate variables for density-dependent component (lambda')
        double Tr = m_thermo->temperature() / Tcrit_i(i);
        double Vr = V_k[i] / Vcrit_i(i);
        double theta_i = thetaShapeFactor(Tr, Vr, m_w_ac[i]);
        double phi_i = phiShapeFactor(Tr, Vr, Zcrit_i(i), m_w_ac[i]);

        f_i[i] = (Tcrit_i(i) / m_ref_Tc)*theta_i; // Equation 12 ely-hanley1983
        h_i[i] = (Vcrit_i(i) / m_ref_Vc)*phi_i; // Equation 13 ely-hanley1983
    }

    double h_m = 0; // Corresponding states parameter, h_x,0 from ely-hanley1983
    double f_m = 0; // Corresponding states parameter, f_x,0 from ely-hanley1983
    double mw_m = 0; // Mixture molecular weight
    double Lambda_1_m = 0.0; // Internal component of mixture thermal conductivity
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j < m_nsp; j++) {
            // Compute the internal contribution to the thermal conductivity of the
            // mixture

            // Equation 3 ely-hanley1983
            double Lambda_1_ij = 2*Lambda_1_i[i]*Lambda_1_i[j] /
                                 (Lambda_1_i[i] + Lambda_1_i[j] + Tiny);
            // Equation 2, ely-hanley1983
            Lambda_1_m += molefracs[i]*molefracs[j]*Lambda_1_ij;

            // Variables for density-dependent translational/collisional component of
            // the mixture.

            // Equation 10, ely-hanley1983
            double f_ij = sqrt(f_i[i]*f_i[j]);

            // Equation 11, ely-hanley1983
            double h_ij = 0.125*pow(pow(h_i[i],1.0/3.0) + pow(h_i[j],1.0/3.0), 3.0);

            // Equation 15, ely-hanley1983
            double mw_ij_inv = 0.5*(m_mw[i] + m_mw[j])/(m_mw[i]*m_mw[j]);

            // Equation 8, ely-hanley1983
            f_m += molefracs[i]*molefracs[j]*f_ij*h_ij;

            // Equation 9, ely-hanley1983
            h_m += molefracs[i]*molefracs[j]*h_ij;

            // Equation, 14 ely-hanley1983
            mw_m += molefracs[i]*molefracs[j]*sqrt(mw_ij_inv*f_ij)*pow(h_ij,-4.0/3.0);
        }
    }

    // The following two equations are the final steps for Equations 8 and 14 in
    // ely-hanley1983 . The calculations in the loop above computed the
    // right-hand-side of the equations, but the left-hand-side of the equations
    // contain other variables that must be moved to the right-hand-side in order to
    // get the values of the variables of interest.
    f_m = f_m/h_m;
    mw_m = pow(mw_m,-2.0)*f_m*pow(h_m,-8.0/3.0);

    // The two relations below are from Equation 7, ely-hanley1983 . This must
    // be in units of g/cm^3 for use with the empirical correlation.
    const double kg_m3_to_g_cm3 = 1e-3; // Conversion factor from kg/m^3 to g/cm^3
    double rho_0 = m_thermo->density()*h_m*kg_m3_to_g_cm3;
    double T_0 = m_temp/f_m;

    // Equation 18, ely-hanley1983
    double Lambda_2_ref = elyHanleyReferenceThermalConductivity(rho_0, T_0);

    // Equation 6, ely-hanley1983
    double F_m = sqrt(f_m*m_ref_MW/mw_m)*pow(h_m,-2.0/3.0);

    // Equation 5, ely-hanley1983
    double Lambda_2_m = F_m*Lambda_2_ref;

    return Lambda_1_m + Lambda_2_m;
}

double HighPressureGasTransport::elyHanleyDilutePureSpeciesViscosity(double V,
    double Tc, double Vc, double Zc, double acentric_factor, double mw)
{
    double Tr = m_thermo->temperature() / Tc;
    double Vr = V / Vc;
    double theta_i = thetaShapeFactor(Tr, Vr, acentric_factor);
    double phi_i = phiShapeFactor(Tr, Vr, Zc, acentric_factor);

    double f_fac = (Tc / m_ref_Tc)*theta_i; // Equation 7 ely-hanley1981
    double h_fac = (Vc / m_ref_Vc)*phi_i;   // Equation 8 ely-hanley1981
    double T_0 = m_temp/f_fac; // Equation 3, ely-hanley1981

    // Dilute reference fluid viscosity correlation, from Table III in
    // ely-hanley1981
    double mu_0 = elyHanleyDiluteReferenceViscosity(T_0);

    // Equation 2, ely-hanley1981
    double F = sqrt(f_fac*(mw/m_ref_MW))*pow(h_fac,-2.0/3.0);

    return mu_0*F;
}

double HighPressureGasTransport::thetaShapeFactor(double Tr, double Vr,
                                                  double acentric_factor)
{
    double T_p = std::min(std::max(Tr,0.5), 2.0);
    double V_p = std::min(std::max(Vr,0.5), 2.0);

    return 1 + (acentric_factor - m_ref_acentric_factor)*(0.090569 - 0.862762*log(T_p)
           + (0.316636 - 0.465684/T_p)*(V_p - 0.5));

}

double HighPressureGasTransport::phiShapeFactor(double Tr, double Vr, double Zc,
                                                double acentric_factor)
{
    double T_p = std::min(std::max(Tr,0.5), 2.0);
    double V_p = std::min(std::max(Vr,0.5), 2.0);

    return (1 + (acentric_factor - m_ref_acentric_factor)*(0.394901*(V_p - 1.023545)
            - 0.932813*(V_p - 0.754639)*log(T_p)))*(m_ref_Zc/Zc);

}

double HighPressureGasTransport::elyHanleyDiluteReferenceViscosity(double T0)
{
    // Conversion factor from the correlation in micrograms/cm/s to Pa*s.
    double const correlation_viscosity_conversion = 1e-7;

    if (T0 > 10000) {
        T0 = 10000; // Limit the temperature to 10000 K
    }

    // Coefficients for the correlation from Table III of ely-hanley1981
    const std::vector<double> c = {2.907741307e6, -3.312874033e6, 1.608101838e6,
                                   -4.331904871e5, 7.062481330e4, -7.116620750e3,
                                   4.325174400e2, -1.445911210e1, 2.037119479e-1};

    double mu_0 = 0.0;
    for (size_t i = 0; i < 9; i++) {
        mu_0 += c[i]*pow(T0,(i+1.0-4.0)/3.0);
    }
    return correlation_viscosity_conversion*mu_0;
}

double HighPressureGasTransport::elyHanleyReferenceThermalConductivity(double rho0,
                                                                       double T0)
{
    // Computing the individual terms of Equation 18, ely-hanley1983 . This is an
    // evaluation of the expressions shown in Table I of ely-hanley1983 . The
    // correlation returns values of thermal conductivity in mW/m/K,
    // so a conversion is needed.
    const double correlation_factor = 1e-3; // mW/m/K to W/m/K

    // This is the reference gas, dilute gas viscosity (eta_0 in Table III of
    // ely-hanley1981)
    double mu_0 = elyHanleyDiluteReferenceViscosity(T0);

    // First term in Equation 18. This expression has the correct units because
    // it does not use any empirical correlation, so it is excluded at the end from
    // the unit conversion.
    double Lambda_ref_star = (15*GasConstant / (4*m_ref_MW))*mu_0;

    // Second term in Equation 18
    const vector<double> b = {-2.52762920e-1, 3.34328590e-1, 1.12, 1.680e2};
    double Lambda_ref_1 = (b[0] + b[1]*pow(b[2] - log(T0/b[3]), 2))*rho0;

    // Third term in Equation 18
    const vector<double> a = {-7.197708227, 8.5678222640e1, 1.2471834689e1,
                              -9.8462522975e2, 3.5946850007e-1, 6.9798412538e1,
                              -8.7288332851e2};
    double delta_lambda_ref = exp(a[0] + a[1]/T0)
                              * (exp((a[2] + a[3]*pow(T0,-1.5))*pow(rho0,0.1)
                              + (rho0/m_ref_rhoc - 1)*sqrt(rho0)*(a[4] + a[5]/T0
                              + a[6]*pow(T0,-2))) - 1.0);

    return Lambda_ref_star + (Lambda_ref_1 + delta_lambda_ref)*correlation_factor;
}

double HighPressureGasTransport::viscosity()
{
    computeMixtureParameters();

    // This is η*ξ
    double nondimensional_viscosity = highPressureNondimensionalViscosity(
        m_Tr_mix, m_Pr_mix, m_FP_mix_o, m_FQ_mix_o,
        m_P_vap_mix, m_Pc_mix);

    // Using equation 9-4.14, with units of 1/(Pa*s)
    double numerator = GasConstant*m_Tc_mix*pow(Avogadro,2.0);
    double denominator = pow(m_MW_mix,3.0)*pow(m_Pc_mix,4.0);
    double xi = pow(numerator / denominator, 1.0/6.0);

    // Return the viscosity in kg/m/s
    return nondimensional_viscosity / xi;
}

void HighPressureGasTransport::computeMixtureParameters()
{
    double Tc_mix = 0.0;
    double Pc_mix_n = 0.0; // Numerator in equation 9-5.19
    double Pc_mix_d = 0.0; // Denominator in equation 9-5.19

    // Equation 9.5.20, Cantera already mole-weights the molecular weights
    double MW_mix = m_thermo->meanMolecularWeight();

    // Mole-fraction-weighted mixture average of the low-pressure polarity correction
    // factor
    double FP_mix_o = 0;

    // Mole-fraction-weighted mixture average of the low-pressure quantum correction
    // factor
    double FQ_mix_o = 0;

    double MW_H = m_mw[0]; // Molecular weight of the heaviest species
    double MW_L = m_mw[0]; // Molecular weight of the lightest species

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

        Tc_mix += Tc*molefracs[i]; // Equation 9-5.18

        Pc_mix_n += molefracs[i]*Zc; // Numerator of 9-5.19
        // (Units of Vcrit_i are m^3/kmol, which are fine because they cancel out
        // with the Cantera gas constant's units used later)
        Pc_mix_d += molefracs[i]*Vcrit_i(i); // Denominator of 9-5.19

        // Calculate ratio of heaviest to lightest species, used in
        // Equation 9-5.23.
        if (m_mw[i] > MW_H) {
            MW_H = m_mw[i];
            x_H = molefracs[i];
        } else if (m_mw[i] < MW_L) {
            MW_L = m_mw[i];
        }

        // Calculate pure-species reduced dipole moment for pure-species
        // polar correction term.
        // Equation 9-4.17 requires the pressure
        // to be in units of bar, so we convert from Pa to bar.
        // The dipole moment is stored in SI units, and it needs to be in
        // units of Debye for the Lucas method.
        double pascals_to_bar = 1e-5;
        double SI_to_Debye = lightSpeed / 1e-21; // Conversion factor from C*m to Debye
        double dipole_ii = m_dipole(i,i)*SI_to_Debye;
        double mu_ri = 52.46*dipole_ii*dipole_ii*(Pcrit_i(i)*pascals_to_bar)/(Tc*Tc);

        // mole-fraction weighting of pure-species polar correction term
        FP_mix_o += molefracs[i] * polarityCorrectionFactor(mu_ri, Tr, Zc);

        // Calculate contribution to quantum correction term.
        // Note:  This assumes the species of interest (He, H2, and D2) have
        //        been named in this specific way.
        vector<string> spnames = m_thermo->speciesNames();
        if (spnames[i] == "He") {
            FQ_mix_o += molefracs[i]*quantumCorrectionFactor(1.38, Tr, m_mw[i]);
        } else if (spnames[i] == "H2") {
            FQ_mix_o += molefracs[i]*(quantumCorrectionFactor(0.76, Tr, m_mw[i]));
        } else if (spnames[i] == "D2") {
            FQ_mix_o += molefracs[i]*(quantumCorrectionFactor(0.52, Tr, m_mw[i]));
        } else {
            FQ_mix_o += molefracs[i];
        }
    }

    double Tr_mix = tKelvin/Tc_mix;
    double Pc_mix = GasConstant*Tc_mix*Pc_mix_n/Pc_mix_d;
    double Pr_mix = m_thermo->pressure()/Pc_mix;

    // Compute the mixture value of the low-pressure quantum correction factor
    // Equation 9-5.23.
    double ratio = MW_H/MW_L;
    double A = 1.0;
    if (ratio > 9 && x_H > 0.05 && x_H < 0.7) {
        A = 1 - 0.01*pow(ratio,0.87);
    }
    FQ_mix_o *= A;

    m_FQ_mix_o = FQ_mix_o;
    m_FP_mix_o = FP_mix_o;
    m_Tc_mix = Tc_mix;
    m_Tr_mix = Tr_mix;
    m_Pc_mix = Pc_mix;
    m_Pr_mix = Pr_mix;
    m_MW_mix = MW_mix;
    m_P_vap_mix = P_vap_mix;
}

double HighPressureGasTransport::lowPressureNondimensionalViscosity(
    double Tr, double FP, double FQ)
{
    double first_term = 0.807*pow(Tr,0.618) - 0.357*exp(-0.449*Tr);
    double second_term = 0.340*exp(-4.058*Tr) + 0.018;
    return (first_term + second_term)*FP*FQ;
}

double HighPressureGasTransport::highPressureNondimensionalViscosity(
    double Tr, double Pr, double FP_low, double FQ_low, double P_vap, double P_crit)
{
    // This is η_0*ξ
    double Z_1 = lowPressureNondimensionalViscosity(Tr, FP_low, FQ_low);

    double Z_2;
    if (Tr <= 1.0) {
        if (Pr < P_vap/P_crit) {
            double alpha = 3.262 + 14.98*pow(Pr, 5.508);
            double beta = 1.390 + 5.746*Pr;
            Z_2 = 0.600 + 0.760*pow(Pr,alpha) + (0.6990*pow(Pr,beta) - 0.60) * (1-Tr);
        } else {
            throw CanteraError(
                "HighPressureGasTransport::highPressureNondimensionalViscosity",
                "State is outside the limits of the Lucas model, Pr ({}) >= "
                "P_vap / P_crit ({}) when Tr ({}) <= 1.0", Pr, P_vap / P_crit, Tr);
        }
    } else if (Tr > 1.0 && Tr < 40.0) {
        if (Pr > 0.0 && Pr <= 100.0) {
            // The following expressions are given in page 9.36 of Poling and
            // correspond to parameters in equation 9-6.8.
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
            throw CanteraError(
                "HighPressureGasTransport::highPressureNondimensionalViscosity",
                "The value of Pr ({}) is outside the limits of the Lucas model, "
                "valid values of Pr are: 0.0 < Pr <= 100", Pr);
        }
    } else {
        throw CanteraError(
            "HighPressureGasTransport::highPressureNondimensionalViscosity",
            "The value of Tr is outside the limits of the Lucas model, "
            "valid  values of Tr are: 1.0 < Tr < 40", Tr);
    }

    double Y = Z_2 / Z_1;
    double FP = (1 + (FP_low - 1)*pow(Y,-3.0)) / FP_low;
    double FQ = (1 + (FQ_low - 1)*(1.0/Y - 0.007*pow(log(Y),4.0))) / FQ_low;

    // Return the non-dimensional viscosity η*ξ
    return Z_2 * FP * FQ;
}

double HighPressureGasTransport::quantumCorrectionFactor(double Q, double Tr,
                                                         double MW)
{
    return 1.22*pow(Q,0.15)*(1 + 0.00385*pow(pow(Tr - 12.0, 2.0), 1.0/MW)
                             *sign(Tr - 12.0));
}

double HighPressureGasTransport::polarityCorrectionFactor(double mu_r, double Tr,
                                                          double Z_crit)
{
    if (mu_r < 0.022) {
        return 1;
    } else if (mu_r < 0.075) {
        return 1 + 30.55*pow(std::max(0.292 - Z_crit,0.0), 1.72);
    } else {
        return 1 + 30.55*pow(std::max(0.292 - Z_crit, 0.0), 1.72)
               *fabs(0.96 + 0.1*(Tr - 0.7));
    }
}


// ChungHighPressureGasTransport Implementation
// --------------------------------------------
void ChungHighPressureGasTransport::init(ThermoPhase* thermo, int mode)
{
    MixTransport::init(thermo, mode);
    initializeCriticalProperties();
    initializePureFluidProperties();
    m_P_corr_ij.resize(m_nsp, m_nsp);
}

void ChungHighPressureGasTransport::initializePureFluidProperties()
{
    // First fill the species-specific values that will then be used in the
    // combining rules for the Chung method.
    m_sigma_i.resize(m_nsp);
    m_epsilon_over_k_i.resize(m_nsp);
    m_acentric_factor_i.resize(m_nsp);
    m_MW_i.resize(m_nsp);
    m_kappa_i.resize(m_nsp);
    for (size_t i = 0; i < m_nsp; i++) {
        // From equation 9-5.32.
        double m3_per_kmol_to_cm3_per_mol = 1e3; // Convert from m^3/kmol to cm^3/mol
        double Vc = Vcrit_i(i) * m3_per_kmol_to_cm3_per_mol;
        m_sigma_i[i] = 0.809*pow(Vc, 1.0/3.0);

        // From equation 9-5.34.
        m_epsilon_over_k_i[i] = Tcrit_i(i)/1.2593;

        // NOTE: The association parameter is assumed to be zero for all species, but
        // is left here for completeness or future revision.
        m_kappa_i[i] = 0.0;

        // These values are available from the base class
        m_acentric_factor_i[i] = m_w_ac[i];
        m_MW_i[i] = m_mw[i];
    }
}

double ChungHighPressureGasTransport::thermalConductivity()
{
    computeMixtureParameters();

    // Compute T_star using equation 9-5.26, using the mixture parameters
    double tKelvin = m_thermo->temperature();
    double T_star = tKelvin / m_epsilon_over_k_mix;

    // The density is required for high-pressure gases.
    // The Chung method requires density to be units of mol/cm^3
    // Use the mixture molecular weight (units of kg/kmol).
    // 1 kmol/m^3 = 1e-3 mol/cm^3
    double kg_per_m3_to_mol_per_cm3 = (1.0 / m_MW_mix)*1e-3;
    double density = m_thermo->density()*kg_per_m3_to_mol_per_cm3;

    // The value of Cv is already a mole-weighted average of the pure species values
    double Cv_mix = m_thermo->cv_mole(); // Units are J/kmol/K

    // This result is in units of W/m/K
    double thermal_conductivity = highPressureThermalConductivity(
        tKelvin, T_star, m_MW_mix, density, Cv_mix, m_Vc_mix,
        m_Tc_mix, m_sigma_mix, m_acentric_factor_mix,
        m_mu_r_mix, m_kappa_mix);

    // Return the thermal conductivity in W/m/K
    return thermal_conductivity;
}

double ChungHighPressureGasTransport::highPressureThermalConductivity(
    double T, double T_star, double MW, double rho, double Cv, double Vc,
    double Tc, double sigma, double acentric_factor, double mu_r,
    double kappa)
{
    // Calculate the low-pressure viscosity using the Chung method (units of Pa*s)
    double viscosity = lowPressureViscosity(T, T_star, MW, acentric_factor, mu_r,
                                            sigma, kappa);


    double M_prime = MW / 1000.0; // Convert kg/kmol to kg/mol

    // Definition of tabulated coefficients for the Chung method, as shown in
    // Table 10-3 on page 10.23.
    static const vector<double> a = {2.44166, -5.0924e-1, 6.6107, 1.4543e1, 7.9274e-1,
                                     -5.8634, 9.1089e1};
    static const vector<double> b = {7.4824e-1, -1.5094, 5.6207, -8.9139, 8.2019e-1,
                                     1.2801e1, 1.2811e2};
    static const vector<double> c = {-9.1858e-1, -4.9991e1, 6.4760e1, -5.6379,
                                     -6.9369e-1, 9.5893, -5.4217e1};
    static const vector<double> d = {1.2172e2, 6.9983e1, 2.7039e1, 7.4344e1, 6.3173,
                                     6.5529e1, 5.2381e2};

    // This is slightly pedantic, but this is done to have the naming convention in the
    // equations used match the variable names in the code. This is equation 10-5.9.
    double B_1 = a[0] + b[0]*acentric_factor + c[0]*pow(mu_r, 4.0) + d[0]*kappa;
    double B_2 = a[1] + b[1]*acentric_factor + c[1]*pow(mu_r, 4.0) + d[1]*kappa;
    double B_3 = a[2] + b[2]*acentric_factor + c[2]*pow(mu_r, 4.0) + d[2]*kappa;
    double B_4 = a[3] + b[3]*acentric_factor + c[3]*pow(mu_r, 4.0) + d[3]*kappa;
    double B_5 = a[4] + b[4]*acentric_factor + c[4]*pow(mu_r, 4.0) + d[4]*kappa;
    double B_6 = a[5] + b[5]*acentric_factor + c[5]*pow(mu_r, 4.0) + d[5]*kappa;
    double B_7 = a[6] + b[6]*acentric_factor + c[6]*pow(mu_r, 4.0) + d[6]*kappa;

    double y = rho*Vc/6.0; // Equation 10-5.6 (with rho = 1/V)

    double G_1 = (1.0 - 0.5*y)/(pow(1.0-y, 3.0)); // Equation  10-5.7

    // Equation 10-5.8
    double G_2 = (B_1*((1.0-exp(-B_4*y)) / y) + B_2*G_1*exp(B_5*y) + B_3*G_1)
                 / (B_1*B_4 + B_2 + B_3);

    double q = 3.586e-3*sqrt(Tc/M_prime) / pow(Vc, 2.0/3.0); // Below equation 10-5.5

    double Tr = T/Tc; // Reduced temperature
    // The following 4 equations are defined below Equation 10-3.14
    double alpha = (Cv / GasConstant) - 3.0/2.0; // GasConstant(R) has units (J/kmol/K)
    double beta = 0.7862 - 0.7109*acentric_factor
                  + 1.3168*acentric_factor*acentric_factor;
    double Z = 2.0 + 10.5*Tr*Tr;
    double psi = 1.0 + alpha*((0.215 + 0.28288*alpha - 1.061*beta + 0.26665*Z)
                              / (0.6366 + beta*Z + 1.061*alpha*beta));

    // Equation 10-5.5
    double lambda = (31.2*viscosity*psi/M_prime)*(1.0/G_2 + B_6*y)
                    + q*B_7*y*y*sqrt(Tr)*G_2;

    // Units are W/m/K
    return lambda;
}

double ChungHighPressureGasTransport::viscosity()
{
    computeMixtureParameters();

    // Compute T_star using equation 9-5.26, using the mixture parameters
    double tKelvin = m_thermo->temperature();
    double T_star = tKelvin / m_epsilon_over_k_mix;

    // The density is required for high-pressure gases.
    // The Chung method requires density to be units of mol/cm^3
    // Use the mixture molecular weight (units of kg/kmol) here.
    // 1 kmol/m^3 = 1e-3 mol/cm^3
    double kg_per_m3_to_mol_per_cm3 = (1.0 / m_MW_mix)*1e-3;
    double molar_density = m_thermo->density()*kg_per_m3_to_mol_per_cm3;

    // This result is in units of micropoise
    double viscosity = highPressureViscosity(T_star, m_MW_mix, molar_density,
                                             m_Vc_mix, m_Tc_mix,
                                             m_acentric_factor_mix,
                                             m_mu_r_mix, m_kappa_mix);

    double micropoise_to_pascals_second = 1e-7;
    return viscosity*micropoise_to_pascals_second;
}

void ChungHighPressureGasTransport::computeMixtureParameters()
{
    // Here we use the combining rules defined on page 9.25.
    // We have ASSUMED that the binary interaction parameters are unity for all species
    // as was done in the Chung method.
    //
    // The sigma & kappa relations can be fully computed in the loop. The other ones
    // require a final division by the final mixture values, and so the quantities in
    // the loop are only the numerators of the equations that are referenced in the
    // comments. After the loop, the final mixture values are computed and stored.

    // Zero out the mixture values of the parameters
    m_sigma_mix = 0.0;
    m_epsilon_over_k_mix = 0.0;
    m_MW_mix = 0.0;
    m_acentric_factor_mix = 0.0;
    m_mu_mix = 0.0;
    m_kappa_mix = 0.0;
    m_Tc_mix = 0.0;
    m_Vc_mix = 0.0;
    m_mu_r_mix = 0.0;

    vector<double> molefracs(m_nsp);
    m_thermo->getMoleFractions(&molefracs[0]);
    for (size_t i = 0; i < m_nsp; i++) {
        for (size_t j = 0; j <m_nsp; j++){
            double sigma_ij = sqrt(m_sigma_i[i]*m_sigma_i[j]); // Equation 9-5.33
            // Equation 9-5.25
            m_sigma_mix += molefracs[i]*molefracs[j]*pow(sigma_ij,3.0);

            // Equation 9-5.35
            double epsilon_over_k_ij = sqrt(m_epsilon_over_k_i[i]*m_epsilon_over_k_i[j]);

            // Equation 9-5.27 (numerator only)
            m_epsilon_over_k_mix += molefracs[i]*molefracs[j]*epsilon_over_k_ij
                                    *pow(sigma_ij,3.0);

            // This equation is raised to the power 2, so what we do here is first
            // store only the double summation into this variable, and then later
            // square the result.
            // Equation 9-5.40
            double MW_ij =  (2*m_MW_i[i]*m_MW_i[j]) / (m_MW_i[i] + m_MW_i[j]);

            // Equation 9-5.28 (double summation only)
            m_MW_mix += molefracs[i]*molefracs[j]*epsilon_over_k_ij*pow(sigma_ij,2.0)
                        *sqrt(MW_ij);

            // Equation 9-5.36
            double acentric_factor_ij = 0.5*(m_acentric_factor_i[i]
                                        + m_acentric_factor_i[j]);

            // Equation 9-5.29
            m_acentric_factor_mix += molefracs[i]*molefracs[j]*acentric_factor_ij
                                     *pow(sigma_ij,3.0);

            // The base class' dipole moment values are in the SI units, so we need to
            // convert to Debye units.
            double SI_to_Debye = lightSpeed / 1e-21; // Conversion from C*m to Debye
            double dipole_ii = m_dipole(i,i)*SI_to_Debye;
            double dipole_jj = m_dipole(j,j)*SI_to_Debye;

            // Equation 9-5.30
            m_mu_mix += molefracs[i]*molefracs[j]*pow(dipole_ii*dipole_jj,2.0)
                        /pow(sigma_ij,3.0);

            // Using equation 9-5.31
            double kappa_ij = sqrt(m_kappa_i[i]*m_kappa_i[j]); // Equation 9-5.39
            m_kappa_mix += molefracs[i]*molefracs[j]*kappa_ij; // Equation 9-5.31
        }
    }

    // Finalize the expressions for the mixture values of the parameters

    // Equation 9-5.25 computed the cube of sigma_mix
    m_sigma_mix = pow(m_sigma_mix, 1.0/3.0);
    m_epsilon_over_k_mix /= pow(m_sigma_mix,3.0);

    // The MW_mix was only the numerator inside the brackets of equation 9-5.28
    m_MW_mix = pow(m_MW_mix/(m_epsilon_over_k_mix*m_sigma_mix*m_sigma_mix), 2.0);

    m_acentric_factor_mix /= pow(m_sigma_mix, 3.0);

    // Equation 9-5.30 computed the 4th power of mu_mix
    m_mu_mix = pow(pow(m_sigma_mix, 3.0)*m_mu_mix, 1.0/4.0);

    // Tc_mix is computed using equation 9-5.44.
    m_Tc_mix = 1.2593*m_epsilon_over_k_mix;

    // Vc_mix is computed using equation 9-5.43.
    m_Vc_mix = pow(m_sigma_mix/0.809, 3.0);

    // mu_r_mix is computed using equation 9-5.42.
    m_mu_r_mix = 131.3*m_mu_mix/sqrt(m_Vc_mix*m_Tc_mix);
}

/**
 * Returns the value of the Neufeld collision integral for a given dimensionless
 * temperature. Implementation of equation 9-4.3.
 * Applicable over the range of 0.3 <= T_star <= 100.
 *
 * @param T_star  Dimensionless temperature (Defined in Equation 9-4.1)
 */
double neufeldCollisionIntegral(double T_star)
{
    double A = 1.16145;
    double B = 0.14874;
    double C = 0.52487;
    double D = 0.77320;
    double E = 2.16178;
    double F = 2.43787;

    return A / pow(T_star, B) + C / exp(D*T_star) + E / exp(F*T_star);
}

double ChungHighPressureGasTransport::lowPressureViscosity(double T, double T_star,
    double MW, double acentric_factor, double mu_r, double sigma, double kappa)
{
    double omega = neufeldCollisionIntegral(T_star); // Equation 9-4.3.

    // Molecular shapes and polarities factor, Equation 9-4.11
    double Fc = 1 - 0.2756*acentric_factor + 0.059035*pow(mu_r, 4.0) + kappa;

    // Equation 9-3.9, multiplied by the Chung factor, Fc
    // (another way of writing 9-4.10 that avoids explicit use of the critical volume
    // in this method)
    double viscosity = Fc*(26.69*sqrt(MW*T)/(sigma*sigma*omega));

    double micropoise_to_pascals_second = 1e-7;
    return micropoise_to_pascals_second*viscosity;
}

double ChungHighPressureGasTransport::highPressureViscosity(double T_star, double MW,
    double rho, double Vc, double Tc, double acentric_factor, double mu_r,
    double kappa)
{
    // Definition of tabulated coefficients for the Chung method, as shown in Table 9-6
    // on page 9.40 of Poling.
    static const vector<double> a = {6.324, 1.210e-3, 5.283, 6.623, 19.745, -1.900,
                                     24.275, 0.7972, -0.2382, 0.06863};
    static const vector<double> b = {50.412, -1.154e-3, 254.209, 38.096, 7.630,
                                     -12.537, 3.450, 1.117, 0.06770, 0.3479};
    static const vector<double> c = {-51.680, -6.257e-3, -168.48, -8.464, -14.354,
                                     4.958, -11.291, 0.01235, -0.8163, 0.5926};
    static const vector<double> d ={1189.0, 0.03728, 3898.0, 31.42, 31.53, -18.15,
                                    69.35, -4.117, 4.025, -0.727};

    // This is slightly pedantic, but this is done to have the naming convention in the
    // equations used match the variable names in the code.
    double E_1 = a[0] + b[0]*acentric_factor + c[0]*pow(mu_r, 4.0) + d[0]*kappa;
    double E_2 = a[1] + b[1]*acentric_factor + c[1]*pow(mu_r, 4.0) + d[1]*kappa;
    double E_3 = a[2] + b[2]*acentric_factor + c[2]*pow(mu_r, 4.0) + d[2]*kappa;
    double E_4 = a[3] + b[3]*acentric_factor + c[3]*pow(mu_r, 4.0) + d[3]*kappa;
    double E_5 = a[4] + b[4]*acentric_factor + c[4]*pow(mu_r, 4.0) + d[4]*kappa;
    double E_6 = a[5] + b[5]*acentric_factor + c[5]*pow(mu_r, 4.0) + d[5]*kappa;
    double E_7 = a[6] + b[6]*acentric_factor + c[6]*pow(mu_r, 4.0) + d[6]*kappa;
    double E_8 = a[7] + b[7]*acentric_factor + c[7]*pow(mu_r, 4.0) + d[7]*kappa;
    double E_9 = a[8] + b[8]*acentric_factor + c[8]*pow(mu_r, 4.0) + d[8]*kappa;
    double E_10 = a[9] + b[9]*acentric_factor + c[9]*pow(mu_r, 4.0) + d[9]*kappa;

    double y = rho*Vc/6.0; // Equation 9-6.20

    double G_1 = (1.0 - 0.5*y)/(pow(1.0-y, 3.0)); // Equation 9-6.21

    // Equation 9-6.22
    double G_2 = (E_1*((1.0-exp(-E_4*y)) / y) + E_2*G_1*exp(E_5*y) + E_3*G_1)
                 / (E_1*E_4 + E_2 + E_3);

    // Equation 9-6.23
    double eta_2 = E_7*y*y*G_2*exp(E_8 + E_9/T_star + E_10/(T_star*T_star));

    double omega = neufeldCollisionIntegral(T_star); // Equation 9-4.3

    // Molecular shapes and polarities factor, Equation 9-4.11
    double Fc = 1 -0.2756*acentric_factor + 0.059035*pow(mu_r, 4.0) + kappa;

    // Renamed eta_star and eta_star_star from the Poling description to eta_1 and
    // eta_2 for naming simplicity.
    // Equation 9-6.19
    double eta_1 = (sqrt(T_star)/omega) * (Fc*(1.0/G_2 + E_6*y)) + eta_2;

    double eta = eta_1 * 36.344 * sqrt(MW*Tc) / pow(Vc, 2.0/3.0); // Equation 9-6.18

    return eta;
}

}
