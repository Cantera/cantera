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
        3839.8896914480647,
        5260.8983334513332,
        5764.7097019695211,
        7786.429533070881,
        10411.474081913055,
        15276.785945165157,
        17900.243436157067,
        22085.482962782506,
        25989.144060372793
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for (int i = 0; i < numSteps; ++i)
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
        -19347891.714810669,
        -14757822.388050893,
        -12593133.605195494,
        -12626837.865623865,
        -12131010.479908356,
        -10322881.86739888,
        - 9573869.8636945337,
        -10260863.826955771,
        -10579827.307551134
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(2);
    for (int i = 0; i < numSteps; ++i)
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
    for (int i = 0; i < numSteps; ++i)
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
        30514.752294683516,
        21514.841983025333,
        14848.02859501992,
        15965.482659621264,
        18272.567242414199,
        24453.517437971925,
        25299.003664716853,
        28474.69918493319,
        30810.094532734405
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
