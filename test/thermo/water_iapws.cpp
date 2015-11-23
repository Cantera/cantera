#include "gtest/gtest.h"
#include "cantera/thermo/WaterPropsIAPWSphi.h"

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
