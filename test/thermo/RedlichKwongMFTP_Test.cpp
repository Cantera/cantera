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
        -4.5735784132470691e+008,
        -4.5734715010829216e+008,
        -4.5733760789206791e+008,
        -4.5732903883366525e+008,
        -4.5732130124096912e+008,
        -4.5731427966336435e+008,
        -4.5730787908411121e+008,
        -4.5730202059007066e+008,
        -4.5729663809807611e+008
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
    const double rho1[6] = {
        1.5870830380619182e+002,
        1.5419384162620102e+002,
        1.5016078232989273e+002,
        1.4651351852180966e+002,
        1.4318315080653846e+002,
        1.4011821957432278e+002
    };
    // Phase change between temperatures 4 & 5:
    const double rho2[6] = {
        6.2669819090204760e+002,
        5.9931065632330956e+002,
        5.6593959797702098e+002,
        5.1995461110601525e+002,
        3.3929302641053914e+002,
        2.7562068824891088e+002
    };
    // Supercritical; no discontinuity in rho values:
    const double rho3[6] = {
        6.8411632182418634e+002,
        6.6686672949843251e+002,
        6.4850120074098390e+002,
        6.2879881554424378e+002,
        6.0746376039603331e+002,
        5.8409057903881308e+002
    };

    for(int i=0; i<6; ++i)
    {
        const double temp = 294 + i*2;
        set_r(0.99);
        test_phase->setState_TP(temp, 5542027.5);
        EXPECT_NEAR(test_phase->density(),rho1[i],1.e-8);

        test_phase->setState_TP(temp, 7389370.);
        EXPECT_NEAR(test_phase->density(),rho2[i],1.e-8);

        test_phase->setState_TP(temp, 9236712.5);
        EXPECT_NEAR(test_phase->density(),rho3[i],1.e-8);
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
