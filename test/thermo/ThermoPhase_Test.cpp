#include "gtest/gtest.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include <vector>

namespace Cantera
{

class ThermoPhase_Fixture : public testing::Test
{
protected:
    ThermoPhase test_phase;
public:
    ThermoPhase_Fixture() {}

    ~ThermoPhase_Fixture() {}

    void initializeElements()
    {
      test_phase.addElement("A", 1.);
      test_phase.addElement("B", 2.);
      test_phase.addElement("C", 3.);
    }
};

class TestThermoMethods : public testing::Test
{
public:
    std::unique_ptr<ThermoPhase> thermo;
    TestThermoMethods() {
        thermo.reset(newPhase("h2o2.xml"));
    }
};

TEST_F(TestThermoMethods, getMoleFractionsByName)
{
    thermo->setMoleFractionsByName("O2:0.2, H2:0.3, AR:0.5");
    compositionMap X = thermo->getMoleFractionsByName();
    EXPECT_DOUBLE_EQ(X["O2"], 0.2);
    EXPECT_DOUBLE_EQ(X["H2"], 0.3);
    EXPECT_DOUBLE_EQ(X["AR"], 0.5);

    thermo->setMoleFractionsByName("OH:1e-9, O2:0.2, h2:0.3, AR:0.5");
    X = thermo->getMoleFractionsByName();
    EXPECT_EQ(X.size(), (size_t) 4);

    X = thermo->getMoleFractionsByName(1e-5);
    EXPECT_EQ(X.size(), (size_t) 3);
}

TEST_F(TestThermoMethods, getMassFractionsByName)
{
    thermo->setMassFractionsByName("O2:0.2, H2:0.3, AR:0.5");
    compositionMap Y = thermo->getMassFractionsByName();
    EXPECT_DOUBLE_EQ(Y["O2"], 0.2);
    EXPECT_DOUBLE_EQ(Y["H2"], 0.3);
    EXPECT_DOUBLE_EQ(Y["AR"], 0.5);

    thermo->setMassFractionsByName("OH:1e-9, O2:0.2, H2:0.3, AR:0.5");
    Y = thermo->getMassFractionsByName();
    EXPECT_EQ(Y.size(), (size_t) 4);

    Y = thermo->getMassFractionsByName(1e-5);
    EXPECT_EQ(Y.size(), (size_t) 3);
}

TEST_F(TestThermoMethods, setState_nan)
{
    double nan = std::numeric_limits<double>::quiet_NaN();
    thermo->setState_TP(500, 12345);
    EXPECT_THROW(thermo->setState_TP(nan, 55555), CanteraError);
    EXPECT_THROW(thermo->setState_TP(555, nan), CanteraError);
    EXPECT_THROW(thermo->setState_HP(nan, 55555), CanteraError);
    EXPECT_THROW(thermo->setState_SV(1234, nan), CanteraError);
    EXPECT_THROW(thermo->setState_TR(555, nan), CanteraError);
}

}
