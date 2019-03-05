#include "gtest/gtest.h"
#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/ThermoFactory.h"


namespace Cantera
{

class RedlichKwongMFTP_Test : public testing::Test
{
public:
    RedlichKwongMFTP_Test() {
        test_phase.reset(newPhase("../data/co2_RK_example.cti"));
    }

    //vary the composition of a co2-h2 mixture:
    void set_r(const double r) {
        vector_fp moleFracs(7);
        moleFracs[0] = r;
        moleFracs[2] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(RedlichKwongMFTP_Test, construct_from_cti)
{
    RedlichKwongMFTP* redlich_kwong_phase = dynamic_cast<RedlichKwongMFTP*>(test_phase.get());
    EXPECT_TRUE(redlich_kwong_phase != NULL);
}

TEST_F(RedlichKwongMFTP_Test, chem_potentials)
{
    test_phase->setState_TP(298.15, 101325.);
    // Chemical potential should increase with increasing co2 mole fraction:
    //      mu = mu_0 + RT ln(gamma_k*X_k).
    // where gamma_k is the activity coefficient.  Run regression test against values calculated using
    // the model.
    const double expected_result[9] = {
        -4.573578072583122e+008,
        -4.573471168532005e+008,
        -4.573375753640399e+008,
        -4.573290069609340e+008,
        -4.573212699618942e+008,
        -4.573142489246118e+008,
        -4.573078488392255e+008,
        -4.573019907983406e+008,
        -4.572966087236250e+008
    };

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(7);
    for(int i=0; i < 9; ++i)
    {
        set_r(xmin + i*dx);
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
    }
}

TEST_F(RedlichKwongMFTP_Test, activityCoeffs)
{
    test_phase->setState_TP(298., 1.);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * 298.;
    vector_fp mu0(7);
    vector_fp activityCoeffs(7);
    vector_fp chemPotentials(7);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);

    for(int i=0; i < numSteps; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getChemPotentials(&chemPotentials[0]);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        test_phase->getStandardChemPotentials(&mu0[0]);
        EXPECT_NEAR(chemPotentials[0], mu0[0] + RT*std::log(activityCoeffs[0] * r), 1.e-6);
        EXPECT_NEAR(chemPotentials[2], mu0[2] + RT*std::log(activityCoeffs[2] * (1-r)), 1.e-6);
    }
}

TEST_F(RedlichKwongMFTP_Test, standardConcentrations)
{
    EXPECT_DOUBLE_EQ(test_phase->pressure()/(test_phase->temperature()*GasConstant), test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(test_phase->pressure()/(test_phase->temperature()*GasConstant), test_phase->standardConcentration(1));
}

TEST_F(RedlichKwongMFTP_Test, activityConcentrations)
{
    // Check to make sure activityConcentration_i == standardConcentration_i * gamma_i * X_i
    vector_fp standardConcs(7);
    vector_fp activityCoeffs(7);
    vector_fp activityConcentrations(7);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);

    for(int i=0; i < 9; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        standardConcs[0] = test_phase->standardConcentration(0);
        standardConcs[2] = test_phase->standardConcentration(2);
        test_phase->getActivityConcentrations(&activityConcentrations[0]);

        EXPECT_NEAR(standardConcs[0] * r * activityCoeffs[0], activityConcentrations[0], 1.e-6);
        EXPECT_NEAR(standardConcs[2] * (1-r) * activityCoeffs[2], activityConcentrations[2], 1.e-6);
    }
}

TEST_F(RedlichKwongMFTP_Test, setTP)
{
    // Check to make sure that the phase diagram is accurately reproduced for a few select isobars

    // All sub-cooled liquid:
    const double p1[6] = {
        1.587112190732014e+002,
        1.541966713372675e+002,
        1.501635359781652e+002,
        1.465162036435630e+002,
        1.431857735462774e+002,
        1.401207850479111e+002
    };
    // Phase change between temperatures 4 & 5:
    const double p2[6] = {
        6.267097216456422e+002,
        5.993217207540168e+002,
        5.659501111117172e+002,
        5.199644273242080e+002,
        3.393007538579040e+002,
        2.756259035569044e+002
    };
    // Supercritical; no discontinuity in rho values:
    const double p3[6] = {
        6.841288400828764e+002,
        6.668789423328959e+002,
        6.485130892980700e+002,
        6.288103574172300e+002,
        6.074749284756613e+002,
        5.841013398471708e+002
    };

    for(int i=0; i<6; ++i)
    {
        const double temp = 294 + i*2;
        set_r(0.99);
        test_phase->setState_TP(temp, 5542027.5);
        EXPECT_NEAR(test_phase->density(),p1[i],1.e-8);

        test_phase->setState_TP(temp, 7389370.);
        EXPECT_NEAR(test_phase->density(),p2[i],1.e-8);

        test_phase->setState_TP(temp, 9236712.5);
        EXPECT_NEAR(test_phase->density(),p3[i],1.e-8);
    }
}

TEST_F(RedlichKwongMFTP_Test, critPropLookup)
{
    // Check to make sure that RedlichKwongMFTP is able to properly calculate a and b
    // pureFluidParameters based on tabulated critical properties
    test_phase.reset(newPhase("../data/co2_RK_lookup.cti"));

    // Check that the critical properties (temperature and pressure) are calculated correctly for
    // pure fluids, both for those with pureFluidParameters provided in the cti file (i.e., h2) and
    // those where the pureFluidParameters are calculated based on the tabulated critical properties
    // (i.e. co2):

    // CO2 - should match tabulated values in critProperties.xml
    set_r(1.0);
    EXPECT_DOUBLE_EQ(test_phase->critTemperature(), 304.2);
    EXPECT_DOUBLE_EQ(test_phase->critPressure(), 7390000);

    // H2
    set_r(0.0);
    EXPECT_NEAR(test_phase->critTemperature(), 33.001, 1.e-3);
    EXPECT_NEAR(test_phase->critPressure(), 1347700, 100);

}
};
