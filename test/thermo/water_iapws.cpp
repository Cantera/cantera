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
        water.setState_TD(T, rho);
        double P1 = water.pressure();
        double T2 = T + 0.001;
        water.setState_TD(T2, rho);
        double P2 = water.pressure();
        return (P2 - P1) / 0.001;
    }

    WaterPropsIAPWS water;
};

// See values on p. 395 and Table 13.1 (p. 486) of Wagner & Pruss
TEST_F(WaterPropsIAPWS_Test, triple_point_liquid)
{
    double T = 273.16;
    double pres = water.psat(T);
    EXPECT_NEAR(pres, 611.654771, 5e-6);
    EXPECT_NEAR(water.density(T, pres, WATER_LIQUID), 999.793, 2e-3);
    EXPECT_NEAR(water.intEnergy_mass(), 0.0, 5e-8);
    EXPECT_NEAR(water.entropy_mass(), 0.0, 5e-10);
    EXPECT_NEAR(water.enthalpy_mass(), 0.611782, 2e-5);
    EXPECT_NEAR(water.cv_mass(), 4217.4, 1e-1);
    EXPECT_NEAR(water.cp_mass(), 4219.9, 1e-1);
}

TEST_F(WaterPropsIAPWS_Test, triple_point_gas)
{
    double T = 273.16;
    double pres = water.psat(T);
    EXPECT_NEAR(water.density(T, pres, WATER_GAS), 4.85458e-3, 2e-8);
    EXPECT_NEAR(water.entropy_mass(), 9155.5, 2e0);
    EXPECT_NEAR(water.enthalpy_mass(), 2500920., 2e+1);
    EXPECT_NEAR(water.cv_mass(), 1418.4, 2e-1);
    EXPECT_NEAR(water.cp_mass(), 1884.4, 2e-1);
}

TEST_F(WaterPropsIAPWS_Test, normal_boiling_point)
{
    double T = 373.124;
    double P = water.psat(T);
    EXPECT_NEAR(P, 101324., 1e0);
    double rho = water.density(T, P, WATER_LIQUID);
    EXPECT_NEAR(rho, 958.368, 2e-3);
    EXPECT_NEAR(water.isothermalCompressibility(), 4.901778826964e-10, 2e-21);

    water.density(T, 1.001 * P, WATER_LIQUID);
    EXPECT_NEAR(water.isothermalCompressibility(), 4.901777129953e-10, 2e-21);

    rho = water.density(T, P, WATER_GAS);
    EXPECT_NEAR(rho, 0.597651, 2e-6);
    EXPECT_NEAR(water.isothermalCompressibility(), 1.003322546019e-05, 2e-17);

    rho = water.density(T, P * 0.999, WATER_GAS);
    EXPECT_NEAR(rho, 0.597043, 2e-6);
    EXPECT_NEAR(water.isothermalCompressibility(), 1.004307955057e-05, 2e-17);
}

TEST_F(WaterPropsIAPWS_Test, saturation_pressure_estimate)
{
    vector<double> TT{273.15, 313.9999, 314.0001, 373.15, 647.25};
    vector<double> psat{611.212, 7722.3, 7675.46, 101007, 2.2093e+07};

    for (size_t i = 0; i < TT.size(); i++) {
        double P = water.psat_est(TT[i]);
        EXPECT_NEAR(P, psat[i], 2e-6 * psat[i]);
    }
}

TEST_F(WaterPropsIAPWS_Test, expansion_coeffs)
{
    vector<double> TT{300.0, 300.0, 700.0};
    vector<double> PP{10.0, 10.0e6, 10.0e6};
    vector<double> alpha{0.003333433139236, -0.02277763412159, 0.002346416506367};
    vector<double> beta{1.000020308917, 1265.572840683, 1.240519803181};
    vector<double> beta_num{1.0000203087, 1265.46651311, 1.2405192843};
    for (size_t i = 0; i < TT.size(); i++) {
        double rho = water.density(TT[i], PP[i], WATER_GAS);
        water.setState_TD(TT[i], rho);
        EXPECT_NEAR(water.coeffThermExp(), alpha[i], 2e-14);
        EXPECT_NEAR(water.coeffPresExp(), beta[i], beta[i] * 2e-12);
        EXPECT_NEAR(dPdT(TT[i], PP[i]) / (461.51805 * rho),
                    beta_num[i], 4e-10 * beta_num[i]);
    }
}
