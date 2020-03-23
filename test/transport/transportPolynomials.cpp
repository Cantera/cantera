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

    void check_viscosity_poly(const std::string speciek, vector_fp visc_coeff_expected, int mode) {
        MixTransport tran;
        tran.init(phase, mode);
        size_t k = phase->speciesIndex(speciek);
        for (size_t i = 0; i<visc_coeff_expected.size(); i++) { 
            EXPECT_FLOAT_EQ(tran.viscosityPolynomials()[k][i], visc_coeff_expected[i]);
        }
    }

    void check_cond_poly(const std::string speciek, vector_fp cond_coeff_expected, int mode) {
        MixTransport tran;
        tran.init(phase, mode);
        size_t k = phase->speciesIndex(speciek);
        for (size_t i = 0; i<cond_coeff_expected.size(); i++) {
            EXPECT_FLOAT_EQ(tran.conductivityPolynomials()[k][i], cond_coeff_expected[i]);
        }
    }

    void check_bindiff_poly(const std::string speciek, const std::string speciej, vector_fp bindiff_coeff_expected, int mode) {
        size_t k = phase->speciesIndex(speciek);
        size_t j = phase->speciesIndex(speciej);
        MixTransport tran;
        tran.init(phase, mode);
        int ic = 0;
        if (k < j) {
            for (size_t kk = 0; kk < k; kk++) {
                for (size_t jj = kk; jj < phase->nSpecies(); jj++) { ic++;}
            }
            for (size_t jj = k; jj < j; jj++) { ic++;}
        } else {
            for (size_t jj = 0; jj < j; jj++) {
                for (size_t kk = jj; kk < phase->nSpecies(); kk++) { ic++;}
            }
            for (size_t kk = j; kk < k; kk++) { ic++;}
        }    
        for (size_t i = 0; i<bindiff_coeff_expected.size(); i++) {
            EXPECT_FLOAT_EQ(tran.binDiffusivityPolynomials()[ic][i], bindiff_coeff_expected[i]);
        }
    }

    ThermoPhase* phase;
};



TEST_F(TransportPolynomialsTest, viscosityPolynomials)
{
    check_viscosity_poly("H2", {-0.00044135741,0.00054143229,-0.00010356811,9.67388e-06,-3.32304e-07}, 0);
    check_viscosity_poly("O2", {-0.00628086,0.0036755258,-0.00069891557,6.0422681e-05,-1.9517881e-06}, 0);
    check_viscosity_poly("H2O", {0.01283557,-0.0070486222,0.0014476475,-0.00012428231,3.86979e-06}, 0);
    check_viscosity_poly("H2", {-15.80391,0.853653,-0.028234165,0.00126934}, 10);
    check_viscosity_poly("O2", {-19.450466,2.6764033,-0.27227259,0.012162949}, 10);
    check_viscosity_poly("H2O", {-14.902224,-0.528783,0.3032122,-0.018641952}, 10);
}

TEST_F(TransportPolynomialsTest, conductivityPolynomials)
{
    check_cond_poly("H2", vector_fp({-1.0486523,0.622758,-0.13650987,0.01318645,-0.00047094451}), 0);
    check_cond_poly("O2", vector_fp({0.1259865,-0.0751831,0.016766921,-0.0016427659,6.02148e-05}), 0);
    check_cond_poly("H2O", vector_fp({ -0.33407864,0.20877084,-0.0484958,0.0049524531,-0.00018565441}), 0);
    check_cond_poly("H2", vector_fp({-1.86328,-0.650749,0.14256814,-0.003901962}), 10);
    check_cond_poly("O2", vector_fp({-13.457014,2.8963447,-0.272143,0.011597346}), 10);
    check_cond_poly("H2O", vector_fp({8.676,-7.62436,1.33959,-0.066989854}), 10);
}

TEST_F(TransportPolynomialsTest, binDiffusivityPolynomials)
{
    check_bindiff_poly("H2", "H2O", vector_fp({-0.00796698, 0.0029745048, -0.00023517416, -3.66864e-06, 9.44696e-07}), 0);
    check_bindiff_poly("O2", "O2", vector_fp({-0.0031017126, 0.0016181327, -0.000286911, 2.347239e-05, -7.01879e-07}), 0);
    check_bindiff_poly("H2O", "O2",  vector_fp({ -18.805519, 5.5532975, -0.48501682, 0.020188987}), 10);
    check_bindiff_poly("H2", "O2", vector_fp({-9.45132, 2.5185864, -0.11599084, 0.0051866062}), 10);
}
