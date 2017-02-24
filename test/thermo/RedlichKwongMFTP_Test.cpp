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
        -4.573578067074649e+008,
        -4.573471163377696e+008,
        -4.573375748803425e+008,
        -4.573290065058332e+008,
        -4.573212695326964e+008,
        -4.573142485189869e+008,
        -4.573078484551440e+008,
        -4.573019904340246e+008,
        -4.572966083775078e+008
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
        1.587029921158317e+002,
        1.541895558698696e+002,
        1.501572815648243e+002,
        1.465106359800041e+002,
        1.431807662747959e+002,
        1.401162435728261e+002
    };
    // Phase change between temperatures 4 & 5:
    const double p2[6] = {
        6.265136821574670e+002,
        5.991027079853330e+002,
        5.656903533839055e+002,
        5.196021189855490e+002,
        3.384435863009947e+002,
        2.755331531855265e+002
    };
    // Supercritical; no discontinuity in rho values:
    const double p3[6] = {
        6.839819449357851e+002,
        6.667277456641792e+002,
        6.483568057147166e+002,
        6.286479753170340e+002,
        6.073051275696215e+002,
        5.839223896051005e+002
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
};
