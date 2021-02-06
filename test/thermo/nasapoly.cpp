#include "gtest/gtest.h"
#include "cantera/thermo/NasaPoly1.h"
#include "cantera/thermo/IdealGasPhase.h"

namespace Cantera
{

// CO2 low-temperature polynomial from GRI 3.0
static double coeffs[] = {2.35677352e+00,
                          8.98459677e-03,
                          -7.12356269e-06,
                          2.45919022e-09,
                          -1.43699548e-13,
                          -4.83719697e+04,
                          9.90105222e+00,
                         };

class NasaPoly1Test : public testing::Test
{
public:
    NasaPoly1Test()
        : poly(200.0, 1000.0, 101325.0, coeffs)
        , tpow_(6) {
    }
protected:
    void set_tpow(double T) {
        tpow_[0] = T;
        tpow_[1] = T*T;
        tpow_[2] = T*T*T;
        tpow_[3] = T*T*T*T;
        tpow_[4] = 1.0/T;
        tpow_[5] = std::log(T);
    }

    void testEquivalent(NasaPoly1& p, NasaPoly1& q) {
        EXPECT_EQ(poly.minTemp(), q.minTemp());
        EXPECT_EQ(poly.maxTemp(), q.maxTemp());
        EXPECT_EQ(poly.refPressure(), q.refPressure());

        double cp_R1, h_RT1, s_R1;
        double cp_R2, h_RT2, s_R2;
        double T = 481.99;
        p.updatePropertiesTemp(T, &cp_R1, &h_RT1, &s_R1);
        q.updatePropertiesTemp(T, &cp_R2, &h_RT2, &s_R2);
        EXPECT_DOUBLE_EQ(cp_R1, cp_R2);
        EXPECT_DOUBLE_EQ(h_RT1, h_RT2);
        EXPECT_DOUBLE_EQ(s_R1, s_R2);
    }

    NasaPoly1 poly;
    vector_fp tpow_;
};

TEST_F(NasaPoly1Test, Initialization)
{
    EXPECT_EQ(poly.minTemp(), 200.0);
    EXPECT_EQ(poly.maxTemp(), 1000.0);
    EXPECT_EQ(poly.refPressure(), 101325.0);
}

TEST_F(NasaPoly1Test, updateProperties)
{
    double cp_R, h_RT, s_R;

    // Reference values calculated using CHEMKIN II
    // Expect agreement to single-precision tolerance
    set_tpow(298.15);
    poly.updateProperties(&tpow_[0], &cp_R, &h_RT, &s_R);
    EXPECT_NEAR(4.46633496, cp_R, 1e-7);
    EXPECT_NEAR(-158.739244, h_RT, 1e-5);
    EXPECT_NEAR(25.7125777, s_R, 1e-6);

    set_tpow(876.54);
    poly.updateProperties(&tpow_[0], &cp_R, &h_RT, &s_R);
    EXPECT_NEAR(6.33029000, cp_R, 1e-7);
    EXPECT_NEAR(-50.3179924, h_RT, 1e-5);
    EXPECT_NEAR(31.5401226, s_R, 1e-6);
}

TEST_F(NasaPoly1Test, updatePropertiesTemp)
{
    double cp_R1, h_RT1, s_R1;
    double cp_R2, h_RT2, s_R2;
    double T = 481.99;

    set_tpow(T);
    poly.updatePropertiesTemp(T, &cp_R1, &h_RT1, &s_R1);
    poly.updateProperties(&tpow_[0], &cp_R2, &h_RT2, &s_R2);

    EXPECT_DOUBLE_EQ(cp_R1, cp_R2);
    EXPECT_DOUBLE_EQ(h_RT1, h_RT2);
    EXPECT_DOUBLE_EQ(s_R1, s_R2);
}


TEST(Nasa9Test, Nasa9Thermo) {
    IdealGasPhase g("gasNASA9.xml", "nasa9");
    size_t nsp = g.nSpecies();
    double pres = 1.0E5;
    vector_fp Xset(nsp, 0.0);
    Xset[0] = 0.5;
    Xset[1] = 0.5;

    vector_fp cp_R(nsp, 0.0);
    vector_fp H_RT(nsp, 0.0);
    vector_fp S_R(nsp, 0.0);

    double T0 = 300.0;
    double dT = 199.0;
    double abstol = 1e-7;

    for (size_t i = 0; i < 15; i++) {
        g.setState_TPX(T0 + i*dT, pres, &Xset[0]);
        g.getEntropy_R(&S_R[0]);
        g.getCp_R(&cp_R[0]);
        g.getEnthalpy_RT(&H_RT[0]);

        EXPECT_NEAR(cp_R[0], cp_R[1], abstol);
        EXPECT_NEAR(cp_R[0], cp_R[2], abstol);
        EXPECT_NEAR(H_RT[0], H_RT[1], abstol);
        EXPECT_NEAR(H_RT[0], H_RT[2], abstol);
        EXPECT_NEAR(S_R[0], S_R[1], abstol);
        EXPECT_NEAR(S_R[0], S_R[2], abstol);
    }
}


} // namespace Cantera

int main(int argc, char** argv)
{
    printf("Running main() from nasapoly.cpp\n");
    Cantera::make_deprecation_warnings_fatal();
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
