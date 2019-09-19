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

TEST_F(TestThermoMethods, setState_AnyMap)
{
    AnyMap state;
    state["temperature"] = 321;
    state["Y"] = "AR: 4, O2: 1.0";
    state["P"] = "5 bar";
    thermo->setState(state);
    EXPECT_DOUBLE_EQ(thermo->temperature(), 321);
    EXPECT_DOUBLE_EQ(thermo->pressure(), 5e5);
    EXPECT_DOUBLE_EQ(thermo->massFraction("O2"), 0.2);

    AnyMap state2;
    state2["P"] = OneAtm;
    state2["enthalpy"] = 0;
    state2["X"]["O2"] = 0.9;
    state2["X"]["AR"] = 0.1;
    thermo->setState(state2);
    EXPECT_DOUBLE_EQ(thermo->pressure(), OneAtm);
    EXPECT_NEAR(thermo->temperature(), 298.15, 1e-6);
    EXPECT_DOUBLE_EQ(thermo->moleFraction("AR"), 0.1);

    AnyMap state3;
    state3["density"] = 10;
    state3["V"] = 0.1;
    state3["mole-fractions"] = "O2: 1.0";
    EXPECT_THROW(thermo->setState(state3), CanteraError);

    AnyMap state4;
    state4["mole-fractions"] = "O2: 1.0";
    thermo->setState(state4);
    EXPECT_DOUBLE_EQ(thermo->pressure(), OneAtm);
    EXPECT_NEAR(thermo->temperature(), 298.15, 1e-6);
}

}
