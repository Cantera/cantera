#include "gtest/gtest.h"
#include "cantera/thermo/NasaPoly1.h"

namespace Cantera
{

// CO2 low-temperature polynomial from GRI 3.0. Note that this order is
// different from the order used by CHEMKIN, with the 1/T and log(T)
// coefficients appearing first.
static double coeffs[] = {-4.83719697e+04,
                          9.90105222e+00,
                          2.35677352e+00,
                          8.98459677e-03,
                          -7.12356269e-06,
                          2.45919022e-09,
                          -1.43699548e-13,
                         };

class NasaPoly1Test : public testing::Test
{
public:
    NasaPoly1Test()
        : poly(0, 200.0, 1000.0, 101325.0, coeffs)
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
        EXPECT_EQ(poly.speciesIndex(), q.speciesIndex());

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
    std::vector<double> tpow_;
};

TEST_F(NasaPoly1Test, Initialization)
{
    EXPECT_EQ(poly.minTemp(), 200.0);
    EXPECT_EQ(poly.maxTemp(), 1000.0);
    EXPECT_EQ(poly.refPressure(), 101325.0);
    EXPECT_EQ(poly.speciesIndex(), (size_t) 0);
}

TEST_F(NasaPoly1Test, Copy)
{
    NasaPoly1 q(poly);
    testEquivalent(poly, q);
}

TEST_F(NasaPoly1Test, Assignment)
{
    NasaPoly1 q;
    q = poly;
    testEquivalent(poly, q);
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

} // namespace Cantera

int main(int argc, char** argv)
{
    printf("Running main() from nasapoly.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
