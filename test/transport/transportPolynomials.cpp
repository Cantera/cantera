#include "gtest/gtest.h"

#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/transport/MixTransport.h"

using namespace Cantera;

class TransportPolynomialsTest : public testing::Test
{
public:
    TransportPolynomialsTest() {
        phase.reset(newPhase("h2o2.yaml"));
        tran.init(phase.get(), 0);
        ck_tran.init(phase.get(), CK_Mode);
    }

    void check_viscosity_poly(const std::string& speciek, const vector_fp& visc_coeff_expected, int cmode) {
        size_t k = phase->speciesIndex(speciek);
        vector_fp coeffs (cmode == CK_Mode ? 4 : 5);
        if (cmode == CK_Mode) {
            ck_tran.getViscosityPolynomial(k, &coeffs[0]);
        } else {
            tran.getViscosityPolynomial(k, &coeffs[0]);
        }
        for (size_t i = 0; i < visc_coeff_expected.size(); i++) {
            EXPECT_NEAR(coeffs[i], visc_coeff_expected[i], 1e-5);
        }
    }

    void check_cond_poly(const std::string& speciek, const vector_fp& cond_coeff_expected, int cmode) {
        MixTransport tran;
        tran.init(phase.get(), cmode);
        size_t k = phase->speciesIndex(speciek);
        vector_fp coeffs (cmode == CK_Mode ? 4 : 5);
        if (cmode == CK_Mode) {
            ck_tran.getConductivityPolynomial(k, &coeffs[0]);
        } else {
            tran.getConductivityPolynomial(k, &coeffs[0]);
        }
        for (size_t i = 0; i < cond_coeff_expected.size(); i++) {
            EXPECT_NEAR(coeffs[i], cond_coeff_expected[i], 1e-5);
        }
    }

    void check_bindiff_poly(const std::string& speciek, const std::string& speciej, const vector_fp& bindiff_coeff_expected, int cmode) {
        size_t k = phase->speciesIndex(speciek);
        size_t j = phase->speciesIndex(speciej);
        vector_fp coeffs (cmode == CK_Mode ? 4 : 5);
        if (cmode == CK_Mode) {
            ck_tran.getBinDiffusivityPolynomial(k, j, &coeffs[0]);
        } else {
            tran.getBinDiffusivityPolynomial(k, j, &coeffs[0]);
        }
        for (size_t i = 0; i < bindiff_coeff_expected.size(); i++) {
            EXPECT_NEAR(coeffs[i], bindiff_coeff_expected[i], 1e-5);
        }
    }

    std::unique_ptr<ThermoPhase> phase;
    MixTransport tran;
    MixTransport ck_tran;
};



TEST_F(TransportPolynomialsTest, viscosityPolynomials)
{
    check_viscosity_poly("H2", {-0.0003286235158, 0.0004740294433, -8.852339014e-05, 8.188000383e-06, -2.775116847e-07}, 0);
    check_viscosity_poly("O2", {-0.006186428072, 0.003618824512, -0.0006861983405, 5.91601271e-05, -1.904977178e-06}, 0);
    check_viscosity_poly("H2O", {0.009686030092, -0.005171558465, 0.001029793739, -8.310407406e-05, 2.354033086e-06}, 0);
    check_viscosity_poly("H2", {-15.74843779, 0.8286882329, -0.0245136295, 0.001085683251}, CK_Mode);
    check_viscosity_poly("O2", {-19.07096138, 2.50635136, -0.2470281853, 0.01092118192}, CK_Mode);
    check_viscosity_poly("H2O", {-15.36038093, -0.324024089, 0.2728840212, -0.01715300254}, CK_Mode);
}

TEST_F(TransportPolynomialsTest, conductivityPolynomials)
{
    check_cond_poly("H2", vector_fp({-0.9677034329, 0.574433766, -0.1257371151, 0.01212356977, -0.0004317820758}), 0);
    check_cond_poly("O2", vector_fp({0.106895521, -0.06376711344, 0.01421775618, -0.001390841153, 5.091744389e-05}), 0);
    check_cond_poly("H2O", vector_fp({-0.4100422806, 0.2548542901, -0.05893633847, 0.005999369421, -0.0002248571694}), 0);
    check_cond_poly("H2", vector_fp({0.4732693668, -1.700244996, 0.2987141261, -0.01159868465}), CK_Mode);
    check_cond_poly("O2", vector_fp({-14.49659725, 3.36457572, -0.3419869392, 0.01504842571}), CK_Mode);
    check_cond_poly("H2O", vector_fp({10.32430752, -8.362450552, 1.449100507, -0.07237437923}), CK_Mode);
}

TEST_F(TransportPolynomialsTest, binDiffusivityPolynomials)
{
    check_bindiff_poly("H2", "H2O", vector_fp({-0.009701326278, 0.004014323899, -0.0004679109588, 1.938085266e-05, 9.241023548e-08}), 0);
    check_bindiff_poly("O2", "O2", vector_fp({-0.003127289937, 0.001633232189, -0.0002902324473, 2.379515419e-05, -7.135754459e-07}), 0);
    check_bindiff_poly("H2O", "O2",  vector_fp({-18.63036291, 5.475482371, -0.4735550509, 0.01962919378}), CK_Mode);
    check_bindiff_poly("H2", "O2", vector_fp({-9.272394946, 2.438367828, -0.1040764365, 0.00460028674}), CK_Mode);
}
