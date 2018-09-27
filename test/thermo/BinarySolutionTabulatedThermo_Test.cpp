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
        -1019148.841268,
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
    }
}

TEST_F(BinarySolutionTabulatedThermo_Test,interp_s)
{
    test_phase->setState_TP(298.15, 101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        3852.587527,
        5260.898245,
        5764.709566,
        7786.429343,
        10411.473830,
        15276.785622,
        17900.243026,
        22085.482446,
        25989.143405
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->entropy_mole(), 1.e-6);
    }
}


TEST_F(BinarySolutionTabulatedThermo_Test,chem_potentials)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        -19327320.552727,
        -14757822.382223,
        -12593133.583222,
        -12626837.825618,
        -12131010.419483,
        -10322881.783439,
        - 9573869.751959,
        -10260863.681331,
        -10579827.118452
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
        30641.731142,
        21514.841963,
        14848.028521,
        15965.482525,
        18272.567039,
        24453.517156,
        25299.003289,
        28474.698696,
        30810.093898
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
