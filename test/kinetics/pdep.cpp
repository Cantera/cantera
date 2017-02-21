#include "gtest/gtest.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"

namespace Cantera
{

class PdepTest : public testing::Test
{
public:
    PdepTest() {}

    static void SetUpTestCase() {
        XML_Node* phase_node = get_XML_File("../data/pdep-test.xml");

        thermo_ = new IdealGasPhase();
        kin_ = new GasKinetics();

        buildSolutionFromXML(*phase_node, "gas", "phase", thermo_, kin_);
    }

    static void TearDownTestCase() {
        delete thermo_;
        thermo_ = NULL;

        delete kin_;
        kin_ = NULL;
    }

    void SetUp() {
        std::string Xref = "H:1.0, R1A:1.0, R1B:1.0, R2:1.0, "
                           "R3:1.0, R4:1.0, R5:1.0, R6:1.0";

        thermo_->setState_TPX(900.0, 101325 * 8.0, Xref);
    }

protected:
    static ThermoPhase* thermo_;
    static Kinetics* kin_;

    void set_TP(double T, double P) {
        T_ = T;
        RT_ = GasConst_cal_mol_K * T;
        P_ = P;
        thermo_->setState_TP(T_, P_);
    }

    double k(double A, double n, double Ea) {
        return A * pow(T_, n) * exp(-Ea/RT_);
    }
    double T_, RT_, P_;
};

ThermoPhase* PdepTest::thermo_ = NULL;
Kinetics* PdepTest::kin_ = NULL;

TEST_F(PdepTest, reactionCounts)
{
    EXPECT_EQ((size_t) 6, kin_->nReactions());
}

TEST_F(PdepTest, PlogLowPressure)
{
    // Test that P-log reactions have the right low-pressure limit
    set_TP(500.0, 1e-7);
    vector_fp kf(6);
    kin_->getFwdRateConstants(&kf[0]);

    // Pre-exponential factor decreases by 10^3 for second-order reaction
    // when converting from cm + mol to m + kmol
    double kf0 = k(1.212400e+13, -0.5779, 10872.7);
    double kf1 = k(1.230000e+05, 1.53, 4737.0);
    double kf2 = k(2.440000e+7, 1.04, 3980.0);
    double kf3 = k(2.889338e-17*(Avogadro/1e6), 1.98, 4521.0);

    EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0);
    EXPECT_NEAR(kf1, kf[1], 1e-9 * kf1);
    EXPECT_NEAR(kf2, kf[2], 1e-9 * kf2);
    EXPECT_NEAR(kf3, kf[3], 1e-9 * kf3);
}

TEST_F(PdepTest, PlogHighPressure)
{
    // Test that P-log reactions have the right high-pressure limit
    set_TP(500.0, 1e10);
    vector_fp kf(6);
    kin_->getFwdRateConstants(&kf[0]);

    // Pre-exponential factor decreases by 10^3 for second-order reaction
    // when converting from cm + mol to m + kmol
    double kf0 = k(5.963200e+53, -11.529, 52599.6);
    double kf3 = k(2.889338e-17*(Avogadro/1e6), 1.98, 4521.0);

    EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0);
    EXPECT_NEAR(kf3, kf[3], 1e-9 * kf3);
}

TEST_F(PdepTest, PlogDuplicatePressures)
{
    // Test that multiple rate expressions are combined when necessary
    set_TP(500.0, 1e10);
    vector_fp kf(6);

    kin_->getFwdRateConstants(&kf[0]);
    double kf1 = k(1.3700e+14, -0.79, 17603.0) + k(1.2800e+03, 1.71, 9774.0);
    double kf2 = k(-7.4100e+27, -5.54, 12108.0) + k(1.9000e+12, -0.29, 8306.0);

    EXPECT_NEAR(kf1, kf[1], 1e-9 * kf1);
    EXPECT_NEAR(kf2, kf[2], 1e-9 * kf2);
}

TEST_F(PdepTest, PlogCornerCases)
{
    // Test rate evaluation at the corner cases where the pressure
    // is exactly of the specified interpolation values
    set_TP(500.0, 101325);
    vector_fp kf(6);
    kin_->getFwdRateConstants(&kf[0]);

    double kf0 = k(4.910800e+28, -4.8507, 24772.8);
    double kf1 = k(1.2600e+17, -1.83, 15003.0) + k(1.2300e+01, 2.68, 6335.0);
    double kf2 = k(3.4600e+9, 0.442, 5463.0);

    EXPECT_NEAR(kf0, kf[0], 1e-9 * kf0);
    EXPECT_NEAR(kf1, kf[1], 1e-9 * kf1);
    EXPECT_NEAR(kf2, kf[2], 1e-9 * kf2);
}

TEST_F(PdepTest, PlogIntermediatePressure1)
{
    set_TP(1100.0, 20*101325);
    vector_fp ropf(6);
    kin_->getFwdRatesOfProgress(&ropf[0]);

    // Expected rates computed using Chemkin
    // ROP increases by 10**3 when converting from mol/cm3 to kmol/m3
    EXPECT_NEAR(3.100682e+05, ropf[0], 1e2);
    EXPECT_NEAR(2.006871e+05, ropf[1], 1e2);
    EXPECT_NEAR(4.468658e+06, ropf[2], 1e2);
    EXPECT_NEAR(1.774796e+06, ropf[3], 1e2);
}

TEST_F(PdepTest, PlogIntermediatePressure2)
{
    thermo_->setState_TP(1100.0, 0.5*101325);
    vector_fp ropf(6);
    kin_->getFwdRatesOfProgress(&ropf[0]);

    EXPECT_NEAR(5.244649e+02, ropf[0], 5e-2);
    EXPECT_NEAR(2.252537e+02, ropf[1], 2e-2);
    EXPECT_NEAR(2.985338e+03, ropf[2], 3e-1);
    EXPECT_NEAR(1.109248e+03, ropf[3], 1e-1);
}

TEST_F(PdepTest, PlogIntermediatePressure3)
{
    thermo_->setState_TP(800.0, 70*101325);
    vector_fp ropf(6);
    kin_->getFwdRatesOfProgress(&ropf[0]);

    EXPECT_NEAR(2.274501e+04, ropf[0], 1e+1);
    EXPECT_NEAR(2.307191e+05, ropf[1], 1e+2);
    EXPECT_NEAR(2.224601e+07, ropf[2], 1e+3);
    EXPECT_NEAR(1.007440e+07, ropf[3], 1e+3);
}

TEST_F(PdepTest, ChebyshevIntermediate1)
{
    // Test Chebyshev rates in the normal interpolation region
    vector_fp kf(6);

    set_TP(1100.0, 20 * 101325);
    kin_->getFwdRateConstants(&kf[0]);
    // Expected rates computed using RMG-py
    EXPECT_NEAR(3.130698657e+06, kf[4], 1e-1);
    EXPECT_NEAR(1.187949573e+00, kf[5], 1e-7);
}

TEST_F(PdepTest, ChebyshevIntermediate2)
{
    // Test Chebyshev rates in the normal interpolation region
    vector_fp kf(6);

    set_TP(400.0, 0.1 * 101325);
    kin_->getFwdRateConstants(&kf[0]);
    // Expected rates computed using RMG-py
    EXPECT_NEAR(1.713599902e+05, kf[4], 1e-3);
    EXPECT_NEAR(9.581780687e-24, kf[5], 1e-31);
}

TEST_F(PdepTest, ChebyshevIntermediateROP)
{
    set_TP(1100.0, 30 * 101325);
    vector_fp ropf(6);
    // Expected rates computed using Chemkin
    kin_->getFwdRatesOfProgress(&ropf[0]);
    EXPECT_NEAR(4.552930e+03, ropf[4], 1e-1);
    EXPECT_NEAR(4.877390e-02, ropf[5], 1e-5);
}

TEST_F(PdepTest, ChebyshevEdgeCases)
{
    vector_fp kf(6);

    // Minimum P
    set_TP(500.0, 1000.0);
    kin_->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(1.225785655e+06, kf[4], 1e-2);

    // Maximum P
    set_TP(500.0, 1.0e7);
    kin_->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(1.580981157e+03, kf[4], 1e-5);

    // Minimum T
    set_TP(300.0, 101325);
    kin_->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(5.405987017e+03, kf[4], 1e-5);

    // Maximum T
    set_TP(2000.0, 101325);
    kin_->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(3.354054351e+07, kf[4], 1e-1);
}

} // namespace Cantera

int main(int argc, char** argv)
{
    printf("Running main() from pdep.cpp\n");
    Cantera::make_deprecation_warnings_fatal();
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
