#include "gtest/gtest.h"

#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/transport/MixTransport.h"

using namespace Cantera;

class TransportPolynomialsTest : public testing::Test
{
public:
    TransportPolynomialsTest() {
        phase = newPhase("gri30.cti");
    }

    void check_viscosity_poly(const std::string& speciek, const vector_fp& visc_coeff_expected, int cmode) {
        MixTransport tran;
        tran.init(phase, cmode);
        size_t k = phase->speciesIndex(speciek);
        vector_fp coeffs (cmode == CK_Mode ? 4 : 5);
        tran.getViscosityPolynomials(k, &coeffs[0]);
        for (size_t i = 0; i < visc_coeff_expected.size(); i++) {
            EXPECT_NEAR(coeffs[i], visc_coeff_expected[i], 1e-5);
        }
    }

    void check_cond_poly(const std::string& speciek, const vector_fp& cond_coeff_expected, int cmode) {
        MixTransport tran;
        tran.init(phase, cmode);
        size_t k = phase->speciesIndex(speciek);
        vector_fp coeffs (cmode == CK_Mode ? 4 : 5);
        tran.getConductivityPolynomials(k, &coeffs[0]);
        for (size_t i = 0; i < cond_coeff_expected.size(); i++) {
            EXPECT_NEAR(coeffs[i], cond_coeff_expected[i], 1e-5);
        }
    }

    void check_bindiff_poly(const std::string& speciek, const std::string& speciej, const vector_fp& bindiff_coeff_expected, int cmode) {
        MixTransport tran;
        tran.init(phase, cmode);
        size_t k = phase->speciesIndex(speciek);
        size_t j = phase->speciesIndex(speciej);
        vector_fp coeffs (cmode == CK_Mode ? 4 : 5);
        tran.getBinDiffusivityPolynomials(k, j, &coeffs[0]);
        for (size_t i = 0; i < bindiff_coeff_expected.size(); i++) {
            EXPECT_NEAR(coeffs[i], bindiff_coeff_expected[i], 1e-5);
        }
    }

    ThermoPhase* phase;
};



TEST_F(TransportPolynomialsTest, viscosityPolynomials)
{
    check_viscosity_poly("H2", {-0.00044135741,0.00054143229,-0.00010356811,9.67388e-06,-3.32304e-07}, 0);
    check_viscosity_poly("O2", {-0.00628086,0.0036755258,-0.00069891557,6.0422681e-05,-1.9517881e-06}, 0);
    check_viscosity_poly("H2O", {0.01283557,-0.0070486222,0.0014476475,-0.00012428231,3.86979e-06}, 0);
    check_viscosity_poly("H2", {-15.80391,0.853653,-0.028234165,0.00126934}, CK_Mode);
    check_viscosity_poly("O2", {-19.450466,2.6764033,-0.27227259,0.012162949}, CK_Mode);
    check_viscosity_poly("H2O", {-14.902224,-0.528783,0.3032122,-0.018641952}, CK_Mode);
}

TEST_F(TransportPolynomialsTest, conductivityPolynomials)
{
    check_cond_poly("H2", vector_fp({-1.0486523,0.622758,-0.13650987,0.01318645,-0.00047094451}), 0);
    check_cond_poly("O2", vector_fp({0.1259865,-0.0751831,0.016766921,-0.0016427659,6.02148e-05}), 0);
    check_cond_poly("H2O", vector_fp({ -0.33407864,0.20877084,-0.0484958,0.0049524531,-0.00018565441}), 0);
    check_cond_poly("H2", vector_fp({-1.86328,-0.650749,0.14256814,-0.003901962}), CK_Mode);
    check_cond_poly("O2", vector_fp({-13.457014,2.8963447,-0.272143,0.011597346}), CK_Mode);
    check_cond_poly("H2O", vector_fp({8.676,-7.62436,1.33959,-0.066989854}), CK_Mode);
}

TEST_F(TransportPolynomialsTest, binDiffusivityPolynomials)
{
    check_bindiff_poly("H2", "H2O", vector_fp({-0.00796698, 0.0029745048, -0.00023517416, -3.66864e-06, 9.44696e-07}), 0);
    check_bindiff_poly("O2", "O2", vector_fp({-0.0031017126, 0.0016181327, -0.000286911, 2.347239e-05, -7.01879e-07}), 0);
    check_bindiff_poly("H2O", "O2",  vector_fp({ -18.805519, 5.5532975, -0.48501682, 0.020188987}), CK_Mode);
    check_bindiff_poly("H2", "O2", vector_fp({-9.45132, 2.5185864, -0.11599084, 0.0051866062}), CK_Mode);
}
