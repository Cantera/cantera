#include "gtest/gtest.h"
#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class BinarySolutionTabulatedThermo_Test : public testing::Test
{
public:
    BinarySolutionTabulatedThermo_Test(){
        test_phase.reset(newPhase("../data/BinarySolutionTabulatedThermo.cti"));
    }

    void set_defect_X(const double x) {
        vector_fp moleFracs(2);
        moleFracs[0] = x;
        moleFracs[1] = 1-x;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(BinarySolutionTabulatedThermo_Test,construct_from_cti)
{
    BinarySolutionTabulatedThermo* BinarySolutionTabulatedThermo_phase = dynamic_cast<BinarySolutionTabulatedThermo*>(test_phase.get());
    EXPECT_TRUE(BinarySolutionTabulatedThermo_phase != NULL);
}

TEST_F(BinarySolutionTabulatedThermo_Test,interp_h)
{
    test_phase->setState_TP(298.15, 101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        -1024991.831815,
        -1512199.970459,
        -2143625.893392,
        -2704188.166163,
        -2840293.936547,
        -1534983.231904,
        -1193196.003622,
        -1184444.702197,
        -1045348.216962,
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->enthalpy_mole(), 1.e-6);
        // enthalpy is temperature-independent in test data file (all species
        // use constant cp model with cp = 0)
        test_phase->setState_TP(310, 101325);
        EXPECT_NEAR(expected_result[i], test_phase->enthalpy_mole(), 1.e-6);
    }
}

TEST_F(BinarySolutionTabulatedThermo_Test,interp_s)
{
    test_phase->setState_TP(298.15, 101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        3839.8896369,
        5260.8982298,
        5764.7095442,
        7786.4293148,
        10411.4737952,
        15276.7855795,
        17900.2429773,
        22085.4823903,
        25989.1433421
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->entropy_mole(), 1.e-6);
        // entropy is temperature-independent in test data file (all species use
        // constant cp model with cp = 0)
        test_phase->setState_TP(330.0, 101325);
        EXPECT_NEAR(expected_result[i], test_phase->entropy_mole(), 1.e-6);
    }
}


TEST_F(BinarySolutionTabulatedThermo_Test,chem_potentials)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        -19347891.6985338,
        -14757822.3571570,
        -12593133.5581558,
        -12626837.8005517,
        -12131010.3944173,
        -10322881.7583731,
        - 9573869.7268930,
        -10260863.6562655,
        -10579827.0933861
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(2);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
    }
}


TEST_F(BinarySolutionTabulatedThermo_Test,mole_fractions)
{
    test_phase->setState_TP(298.15,101325.);
    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp molefracs(2);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getMoleFractions(&molefracs[0]);
        EXPECT_NEAR(xmin + i*dx, molefracs[0], 1.e-6);
    }
}

TEST_F(BinarySolutionTabulatedThermo_Test,partialMolarEntropies)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        30514.7522401,
        21514.8418794,
        14848.0284372,
        15965.4824414,
        18272.5669557,
        24453.5170723,
        25299.0032059,
        28474.6986124,
        30810.0938144
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp partialMolarEntropies(2);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getPartialMolarEntropies(&partialMolarEntropies[0]);
        EXPECT_NEAR(expected_result[i], partialMolarEntropies[0], 1.e-6);
    }
}
}
