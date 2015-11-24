#include "gtest/gtest.h"
#include "cantera/base/ct_defs.h"
#include "cantera/thermo/WaterPropsIAPWSphi.h"
#include "cantera/thermo/WaterPropsIAPWS.h"

using namespace Cantera;

const double T_c = 647.096;
const double Rho_c = 322.0;

class WaterPropsIAPWSphi_Test : public testing::Test, public WaterPropsIAPWSphi
{
};

// Test agreement with values given in Table 6.6 of: W. Wagner, A. Pruss, "The
// IAPWS Formulation 1995 for the Thermodynamic Properties of Ordinary Water
// Substance for General and Scientific Use," J. Phys. Chem. Ref. Dat, 31, 387,
// 2002.

TEST_F(WaterPropsIAPWSphi_Test, check1) {
    tdpolycalc(T_c / 500.0, 838.025 / Rho_c);
    EXPECT_NEAR(phi0(), 2.047977334796e+00, 1e-11);
    EXPECT_NEAR(phiR(), -3.426932056816e+00, 1e-11);
    EXPECT_NEAR(phi0_d(), 3.842367471137e-01, 1e-11);
    EXPECT_NEAR(phiR_d(), -3.643666503639e-01, 1e-11);
    EXPECT_NEAR(phi0_dd(), -1.476378778326e-01, 1e-11);
    EXPECT_NEAR(phiR_dd(), 8.560637009746e-01, 1e-11);
    EXPECT_NEAR(phi0_t(), 9.046111061752e+00, 1e-11);
    EXPECT_NEAR(phiR_t(), -5.814034352384e+00, 1e-11);
    EXPECT_NEAR(phi0_tt(), -1.932491850131e+00, 1e-11);
    EXPECT_NEAR(phiR_tt(), -2.234407368843e+00, 1e-11);
    EXPECT_NEAR(phi0_dt(), 0.000000000000e+00, 1e-11);
    EXPECT_NEAR(phiR_dt(), -1.121769146703e+00, 1e-11);
}

TEST_F(WaterPropsIAPWSphi_Test, check2) {
    tdpolycalc(T_c / 647.0, 358.0 / Rho_c);
    EXPECT_NEAR(phi0(), -1.563196050525e+00, 1e-11);
    EXPECT_NEAR(phiR(), -1.212026565041e+00, 1e-11);
    EXPECT_NEAR(phi0_d(), 8.994413407821e-01, 1e-11);
    EXPECT_NEAR(phiR_d(), -7.140120243713e-01, 1e-11);
    EXPECT_NEAR(phi0_dd(), -8.089947255079e-01, 1e-11);
    EXPECT_NEAR(phiR_dd(), 4.757306956457e-01, 1e-11);
    EXPECT_NEAR(phi0_t(), 9.803439179390e+00, 1e-11);
    EXPECT_NEAR(phiR_t(), -3.217225007752e+00, 1e-11);
    EXPECT_NEAR(phi0_tt(), -3.433163341431e+00, 1e-11);
    EXPECT_NEAR(phiR_tt(), -9.960295065593e+00, 1e-11);
    EXPECT_NEAR(phi0_dt(), 0.000000000000e+00, 1e-11);
    EXPECT_NEAR(phiR_dt(), -1.332147204361e+00, 1e-11);
}

class WaterPropsIAPWS_Test : public testing::Test
{
public:
    double dPdT(double T, double P) {
        double rho = water.density(T, P);
        water.setState_TR(T, rho);
        double P1 = water.pressure();
        double T2 = T + 0.001;
        water.setState_TR(T2, rho);
        double P2 = water.pressure();
        return (P2 - P1) / 0.001;
    }

    WaterPropsIAPWS water;
};

// See values on p. 395 of Wagner & Pruss.
TEST_F(WaterPropsIAPWS_Test, triple_point_liquid)
{
    double T = 273.16;
    double pres = water.psat(T);
    EXPECT_NEAR(pres, 611.655, 2e-3);
    EXPECT_NEAR(water.density(T, pres, WATER_LIQUID), 999.793, 2e-3);
    EXPECT_NEAR(water.intEnergy(), 0.0, 5e-7);
    EXPECT_NEAR(water.entropy(), 0.0, 5e-9);
    EXPECT_NEAR(water.enthalpy(), 11.0214, 2e-4);
    EXPECT_NEAR(water.Gibbs(), 11.0214, 2e-4);
    EXPECT_NEAR(water.cv(), 75978.2, 2e-1);
    EXPECT_NEAR(water.cp(), 76022.8, 2e-1);
}

TEST_F(WaterPropsIAPWS_Test, triple_point_gas)
{
    double T = 273.16;
    double pres = water.psat(T);
    EXPECT_NEAR(water.density(T, pres, WATER_GAS), 4.85458e-3, 2e-8);
    EXPECT_NEAR(water.intEnergy(), 4.27848e7, 2e2);
    EXPECT_NEAR(water.entropy(), 164939., 2e0);
    EXPECT_NEAR(water.enthalpy(), 4.50547e7, 2e2);
    EXPECT_NEAR(water.Gibbs(), 11.0214, 2e-4);
    EXPECT_NEAR(water.cv(), 25552.6, 2e-1);
    EXPECT_NEAR(water.cp(), 33947.1, 2e-1);
}

TEST_F(WaterPropsIAPWS_Test, normal_boiling_point)
{
    double T = 373.124;
    double P = water.psat(T);
    EXPECT_NEAR(P, 101324., 1e0);
    double rho = water.density(T, P, WATER_LIQUID);
    EXPECT_NEAR(rho, 958.368, 2e-3);
    EXPECT_NEAR(water.isothermalCompressibility(), 4.901779037782e-10, 2e-21);

    water.density(T, 1.001 * P, WATER_LIQUID);
    EXPECT_NEAR(water.isothermalCompressibility(), 4.901777340771e-10, 2e-21);

    rho = water.density(T, P, WATER_GAS);
    EXPECT_NEAR(rho, 0.597651, 2e-6);
    EXPECT_NEAR(water.isothermalCompressibility(), 1.003322591472e-05, 2e-17);

    rho = water.density(T, P * 0.999, WATER_GAS);
    EXPECT_NEAR(rho, 0.597043, 2e-6);
    EXPECT_NEAR(water.isothermalCompressibility(), 1.004308000545e-05, 2e-17);
}

TEST_F(WaterPropsIAPWS_Test, saturation_pressure_estimate)
{
    vector_fp TT{273.15, 313.9999, 314.0001, 373.15, 647.25};
    vector_fp psat{611.212, 7722.3, 7675.46, 101007, 2.2093e+07};

    for (size_t i = 0; i < TT.size(); i++) {
        double P = water.psat_est(TT[i]);
        EXPECT_NEAR(P, psat[i], 2e-6 * psat[i]);
    }
}

TEST_F(WaterPropsIAPWS_Test, expansion_coeffs)
{
    vector_fp TT{300.0, 300.0, 700.0};
    vector_fp PP{10.0, 10.0e6, 10.0e6};
    vector_fp alpha{0.003333433139236, -0.02277763412159, 0.002346416555069};
    vector_fp beta{1.000020308917, 1265.572840683, 1.240519813089};
    vector_fp beta_num{1.0000203087, 1265.46651311, 1.240519294};
    for (size_t i = 0; i < TT.size(); i++) {
        double rho = water.density(TT[i], PP[i], WATER_GAS);
        water.setState_TR(TT[i], rho);
        EXPECT_NEAR(water.coeffThermExp(), alpha[i], 2e-14);
        EXPECT_NEAR(water.coeffPresExp(), beta[i], beta[i] * 2e-12);
        EXPECT_NEAR(dPdT(TT[i], PP[i]) * 18.015268 / (8.314371E3 * rho),
                    beta_num[i], 2e-10 * beta_num[i]);
    }
}
