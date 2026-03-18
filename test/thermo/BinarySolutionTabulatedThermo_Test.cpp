#include "gtest/gtest.h"
#include "cantera/thermo/BinarySolutionTabulatedThermo.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class BinarySolutionTabulatedThermo_Test : public testing::Test
{
public:
    BinarySolutionTabulatedThermo_Test(){
        test_phase = newThermo("../data/BinarySolutionTabulatedThermo.yaml");
    }

    void set_defect_X(const double x) {
        vector<double> moleFracs(2);
        moleFracs[0] = x;
        moleFracs[1] = 1-x;
        test_phase->setMoleFractions(moleFracs);
    }

    shared_ptr<ThermoPhase> test_phase;
};

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
        3839.8897013463629,
        5260.8983501505654,
        5764.7097249117705,
        7786.4295622389736,
        10411.474117668258,
        15276.785988429341,
        17900.243488318705,
        22085.483025077028,
        25989.144133134443
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
        -19347891.744322445,
        -14757822.415520553,
        -12593133.63125352,
        -12626837.890922677,
        -12131010.504991682,
        -10322881.892878814,
        -9573869.8901660107,
        -10260863.854728648,
        -10579827.336476315
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector<double> chemPotentials(2);
    for (int i = 0; i < numSteps; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getChemPotentials(chemPotentials);
        EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
    }
}


TEST_F(BinarySolutionTabulatedThermo_Test,partialMolarEntropies)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        30514.752393666495,
        21514.842075159024,
        14848.028682418964,
        15965.482744473897,
        18272.567326544096,
        24453.517523432041,
        25299.003753502617,
        28474.699278083881,
        30810.094629749936
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector<double> partialMolarEntropies(2);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getPartialMolarEntropies(partialMolarEntropies);
        EXPECT_NEAR(expected_result[i], partialMolarEntropies[0], 1.e-6);
    }
}

TEST_F(BinarySolutionTabulatedThermo_Test,molarVolumes)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        0.03531501777842358,
        0.035715748862103429,
        0.03590414327870764,
        0.035968621429308907,
        0.035977245280539603,
        0.035995403732700486,
        0.036093852117078863,
        0.036325488894662347,
        0.036697196991506385
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->molarVolume(), 1.e-6);
    }
}

TEST_F(BinarySolutionTabulatedThermo_Test,partialMolarVolumes)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        0.04120724363741,
        0.03853288221791,
        0.03693536558788,
        0.03618236414389,
        0.03599080437984,
        0.03628136773515,
        0.03690395850931,
        0.03756972764230,
        0.03802279519842
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector<double> partialMolarVolumes(2);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getPartialMolarVolumes(partialMolarVolumes);
        EXPECT_NEAR(expected_result[i], partialMolarVolumes[0], 1.e-8);
    }
}

TEST_F(BinarySolutionTabulatedThermo_Test,calcDensity)
{
    test_phase->setState_TP(298.15,101325.);
    // These expected results are purely a regression test
    const double expected_result[9] = {
        2060.3132768194214,
        2052.9843930502343,
        2057.9170884664422,
        2069.9048793494585,
        2085.0818181061941,
        2099.6951600056354,
        2109.590568305415,
        2111.6611870644724,
        2105.6376599521886
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for (int i = 0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->density(), 1.e-6);
    }
}
}
