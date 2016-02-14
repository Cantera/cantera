#include "gtest/gtest.h"
#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class RedlichKister_Test : public testing::Test
{
public:
    RedlichKister_Test() {
        test_phase.reset(newPhase("../data/RedlichKisterVPSSTP_valid.xml"));
    }

    void set_r(const double r) {
        vector_fp moleFracs(2);
        moleFracs[0] = r;
        moleFracs[1] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(RedlichKister_Test, construct_from_xml)
{
    RedlichKisterVPSSTP* redlich_kister_phase = dynamic_cast<RedlichKisterVPSSTP*>(test_phase.get());
    EXPECT_TRUE(redlich_kister_phase != NULL);
}

TEST_F(RedlichKister_Test, chem_potentials)
{
    test_phase->setState_TP(298.15, 101325.);

    const double expected_result[9] = {
        -1.2791500420236044e+007,
        -1.2618554504124604e+007,
        -1.2445418272766629e+007,
        -1.2282611679165890e+007,
        -1.2134110753109487e+007,
        -1.1999465396970615e+007,
        -1.1882669410525253e+007,
        -1.1792994839484975e+007,
        -1.1730895987035934e+007
    };

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        set_r(xmin + i*dx);
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, dlnActivities)
{
    test_phase->setState_TP(298.15, 101325.);

    const double expected_result[9] = {
        0.0907127,
        0.200612,
        0.229316,
        0.193278,
        0.142257,
        0.0766133,
        -0.0712113,
        -0.309379,
        -0.492206
    };

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp dlnActCoeffdx(2);
    for(int i=0; i < 9; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getdlnActCoeffdlnX_diag(&dlnActCoeffdx[0]);
        EXPECT_NEAR(expected_result[i], dlnActCoeffdx[0], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, activityCoeffs)
{
    test_phase->setState_TP(298., 1.);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * 298.;
    vector_fp mu0(2);
    vector_fp activityCoeffs(2);
    vector_fp chemPotentials(2);
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
        EXPECT_NEAR(chemPotentials[1], mu0[1] + RT*std::log(activityCoeffs[1] * (1-r)), 1.e-6);
    }
}

TEST_F(RedlichKister_Test, standardConcentrations)
{
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(1));
}

TEST_F(RedlichKister_Test, activityConcentrations)
{
    // Check to make sure activityConcentration_i == standardConcentration_i * gamma_i * X_i
    vector_fp standardConcs(2);
    vector_fp activityCoeffs(2);
    vector_fp activityConcentrations(2);
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
        standardConcs[1] = test_phase->standardConcentration(1);
        test_phase->getActivityConcentrations(&activityConcentrations[0]);

        EXPECT_NEAR(standardConcs[0] * r * activityCoeffs[0], activityConcentrations[0], 1.e-6);
        EXPECT_NEAR(standardConcs[1] * (1-r) * activityCoeffs[1], activityConcentrations[1], 1.e-6);
    }
}

};
