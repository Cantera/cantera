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

TEST(TestMixtureMethods, getSet_EquilRatio_MixtureFraction)
{
    auto sol = newSolution("gri30.yaml", "gri30", "None");
    auto pgas = sol->thermo();
    auto& gas = *pgas;

    // start with some fuel and oxidizer compositions
    compositionMap fuel;
    compositionMap ox;
    fuel["CH4"] = 0.2;
    fuel["O2"] = 0.02;
    fuel["N2"] = 0.1;
    fuel["CO"] = 0.05;
    fuel["CO2"] = 0.02;
    ox["O2"] = 0.21;
    ox["N2"] = 0.79;
    ox["CO"] = 0.04;
    ox["CH4"] = 0.01;
    ox["CO2"] = 0.03;

    gas.setState_TPX(300.0, 1e5, fuel);
    double Y_Cf = gas.elementalMassFraction(gas.elementIndex("C"));
    double Y_Of = gas.elementalMassFraction(gas.elementIndex("O"));
    gas.setState_TPX(300.0, 1e5, ox);
    double Y_Co = gas.elementalMassFraction(gas.elementIndex("C"));
    double Y_Oo = gas.elementalMassFraction(gas.elementIndex("O"));

    // set equivalence ratio to 1.3
    gas.setEquivalenceRatio(1.3, fuel, ox, ThermoBasisType::molar);

    // set mixture to  burnt state to make sure that equivalence ratio and
    // mixture fraction are independent of reaction progress
    gas.equilibrate("HP");
    double phi = gas.getEquivalenceRatio(fuel, ox, ThermoBasisType::molar);
    double phi_loc = gas.getEquivalenceRatio();
    double mf = gas.getMixtureFraction(fuel, ox, ThermoBasisType::molar, "Bilger");
    double mf_C = gas.getMixtureFraction(fuel, ox, ThermoBasisType::molar, "C");
    double mf_O = gas.getMixtureFraction(fuel, ox, ThermoBasisType::molar, "O");
    double l = gas.getStoichAirFuelRatio(fuel, ox, ThermoBasisType::molar);

    EXPECT_NEAR(phi, 1.3, 1e-4);
    EXPECT_NEAR(phi_loc, 1.1726068608, 1e-4);
    EXPECT_NEAR(mf, 0.13415725911, 1e-4);
    EXPECT_NEAR(mf_C, (gas.elementalMassFraction(gas.elementIndex("C"))-Y_Co)/(Y_Cf-Y_Co), 1e-4);
    EXPECT_NEAR(mf_O, (gas.elementalMassFraction(gas.elementIndex("O"))-Y_Oo)/(Y_Of-Y_Oo), 1e-4);
    EXPECT_NEAR(l, 6.5972850678733, 1e-4);

    // set mixture according to mixture fraction
    gas.setMixtureFraction(mf, fuel, ox, ThermoBasisType::molar);
    gas.equilibrate("HP");
    phi = gas.getEquivalenceRatio(fuel, ox, ThermoBasisType::molar);
    phi_loc = gas.getEquivalenceRatio();
    mf = gas.getMixtureFraction(fuel, ox, ThermoBasisType::molar, "Bilger");
    mf_C = gas.getMixtureFraction(fuel, ox, ThermoBasisType::molar, "C");
    mf_O = gas.getMixtureFraction(fuel, ox, ThermoBasisType::molar, "O");
    l = gas.getStoichAirFuelRatio(fuel, ox, ThermoBasisType::molar);
    double p = gas.pressure(); // make sure the pressure has not been altered

    EXPECT_NEAR(phi, 1.3, 1e-4);
    EXPECT_NEAR(phi_loc, 1.1726068608, 1e-4);
    EXPECT_NEAR(mf, 0.13415725911, 1e-4);
    EXPECT_NEAR(mf_C, (gas.elementalMassFraction(gas.elementIndex("C"))-Y_Co)/(Y_Cf-Y_Co), 1e-4);
    EXPECT_NEAR(mf_O, (gas.elementalMassFraction(gas.elementIndex("O"))-Y_Oo)/(Y_Of-Y_Oo), 1e-4);
    EXPECT_NEAR(l, 6.5972850678733, 1e-4);
    EXPECT_NEAR(p, 1e5, 1e-4);

    // do the same for mass fractions as input

    // convert fuel and oxidizer compositions to (non-normalized) mass fractions
    gas.setState_TPX(300, 1e5, fuel);
    fuel.clear();
    for (size_t i=0; i!=gas.nSpecies(); ++i)
        fuel[gas.speciesName(i)] = gas.massFraction(i)*3;

    gas.setState_TPX(300, 1e5, ox);
    ox.clear();
    for (size_t i=0; i!=gas.nSpecies(); ++i)
        ox[gas.speciesName(i)] = gas.massFraction(i)*7;

    gas.setState_TPY(300.0, 1e5, fuel);
    Y_Cf = gas.elementalMassFraction(gas.elementIndex("C"));
    Y_Of = gas.elementalMassFraction(gas.elementIndex("O"));
    gas.setState_TPY(300.0, 1e5, ox);
    Y_Co = gas.elementalMassFraction(gas.elementIndex("C"));
    Y_Oo = gas.elementalMassFraction(gas.elementIndex("O"));

    gas.setEquivalenceRatio(1.3, fuel, ox, ThermoBasisType::mass);

    gas.equilibrate("HP");

    phi = gas.getEquivalenceRatio(fuel, ox, ThermoBasisType::mass);
    phi_loc = gas.getEquivalenceRatio();
    mf = gas.getMixtureFraction(fuel, ox, ThermoBasisType::mass, "Bilger");
    mf_C = gas.getMixtureFraction(fuel, ox, ThermoBasisType::mass, "C");
    mf_O = gas.getMixtureFraction(fuel, ox, ThermoBasisType::mass, "O");
    l = gas.getStoichAirFuelRatio(fuel, ox, ThermoBasisType::mass);

    EXPECT_NEAR(phi, 1.3, 1e-4);
    EXPECT_NEAR(phi_loc, 1.1726068608, 1e-4);
    EXPECT_NEAR(mf, 0.13415725911, 1e-4);
    EXPECT_NEAR(mf_C, (gas.elementalMassFraction(gas.elementIndex("C"))-Y_Co)/(Y_Cf-Y_Co), 1e-4);
    EXPECT_NEAR(mf_O, (gas.elementalMassFraction(gas.elementIndex("O"))-Y_Oo)/(Y_Of-Y_Oo), 1e-4);
    EXPECT_NEAR(l, 6.5972850678733, 1e-4);

    gas.setMixtureFraction(mf, fuel, ox, ThermoBasisType::mass);

    gas.equilibrate("HP");

    phi = gas.getEquivalenceRatio(fuel, ox, ThermoBasisType::mass);
    phi_loc = gas.getEquivalenceRatio();
    mf = gas.getMixtureFraction(fuel, ox, ThermoBasisType::mass, "Bilger");
    mf_C = gas.getMixtureFraction(fuel, ox, ThermoBasisType::mass, "C");
    mf_O = gas.getMixtureFraction(fuel, ox, ThermoBasisType::mass, "O");
    l = gas.getStoichAirFuelRatio(fuel, ox, ThermoBasisType::mass);

    p = gas.pressure(); // make sure the pressure has not been altered

    EXPECT_NEAR(phi, 1.3, 1e-4);
    EXPECT_NEAR(phi_loc, 1.1726068608, 1e-4);
    EXPECT_NEAR(mf, 0.13415725911, 1e-4);
    EXPECT_NEAR(mf_C, (gas.elementalMassFraction(gas.elementIndex("C"))-Y_Co)/(Y_Cf-Y_Co), 1e-4);
    EXPECT_NEAR(mf_O, (gas.elementalMassFraction(gas.elementIndex("O"))-Y_Oo)/(Y_Of-Y_Oo), 1e-4);
    EXPECT_NEAR(l, 6.5972850678733, 1e-4);
    EXPECT_NEAR(p, 1e5, 1e-4);

    // test some special cases

    vector_fp v_ox(gas.nSpecies());
    vector_fp v_fuel(gas.nSpecies());
    v_ox[gas.speciesIndex("O2")] = 21.0;
    v_ox[gas.speciesIndex("N2")] = 79.0;
    v_fuel[gas.speciesIndex("CH4")] = 10.0;

    // special case 1: pure oxidizer
    gas.setState_TPX(300.0, 1e5, v_ox.data());
    EXPECT_NEAR(gas.getEquivalenceRatio(v_fuel.data(), v_ox.data(), ThermoBasisType::mass), 0.0, 1e-4);
    EXPECT_NEAR(gas.getEquivalenceRatio(v_fuel.data(), v_ox.data(), ThermoBasisType::molar), 0.0, 1e-4);
    EXPECT_NEAR(gas.getEquivalenceRatio(), 0.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::molar, "Bilger"), 0.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::mass, "Bilger"), 0.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::molar, "C"), 0.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::mass, "C"), 0.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::molar, "O"), 0.0, 1e-4);

    double Yo2 = gas.massFraction(gas.speciesIndex("O2"));
    double Ych4 = gas.massFraction(gas.speciesIndex("CH4"));
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setMixtureFraction(0.0, v_fuel.data(), v_ox.data(), ThermoBasisType::molar);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setEquivalenceRatio(0.0, v_fuel.data(), v_ox.data(), ThermoBasisType::molar);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);

    gas.setState_TPY(300.0, 1e5, v_ox.data());
    Yo2 = gas.massFraction(gas.speciesIndex("O2"));
    Ych4 = gas.massFraction(gas.speciesIndex("CH4"));
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setEquivalenceRatio(0.0, v_fuel.data(), v_ox.data(), ThermoBasisType::mass);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setMixtureFraction(0.0, v_fuel.data(), v_ox.data(), ThermoBasisType::mass);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);

    // special case 2: pure fuel
    gas.setState_TPX(300.0, 1e5, v_fuel.data());
    ASSERT_EQ(gas.getEquivalenceRatio(v_fuel.data(), v_ox.data(), ThermoBasisType::mass) > 1e10, true);
    ASSERT_EQ(gas.getEquivalenceRatio(v_fuel.data(), v_ox.data(), ThermoBasisType::molar) > 1e10, true);
    ASSERT_EQ(gas.getEquivalenceRatio() > 1e10, true);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::molar, "Bilger"), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::mass, "Bilger"), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::molar, "C"), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::mass, "C"), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::molar, "O"), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(v_fuel.data(), v_ox.data(), ThermoBasisType::mass, "O"), 1.0, 1e-4);

    Yo2 = gas.massFraction(gas.speciesIndex("O2"));
    Ych4 = gas.massFraction(gas.speciesIndex("CH4"));
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setMixtureFraction(1.0, v_fuel.data(), v_ox.data(), ThermoBasisType::mass);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setMixtureFraction(1.0, v_fuel.data(), v_ox.data(), ThermoBasisType::molar);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setEquivalenceRatio(1e10, v_fuel.data(), v_ox.data(), ThermoBasisType::mass);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
    gas.setState_TPX(300.0, 1e5, "N2:1");
    gas.setEquivalenceRatio(1e10, v_fuel.data(), v_ox.data(), ThermoBasisType::molar);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("O2")), Yo2, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);

    // special case 3: stoichiometric mixture, input is string
    std::string sfuel = "CH4";
    std::string sox = "O2:21,N2:79";
    gas.setEquivalenceRatio(1.0, sfuel, sox, ThermoBasisType::molar);
    gas.setMixtureFraction(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::molar, "Bilger"), sfuel, sox, ThermoBasisType::molar);
    EXPECT_NEAR(gas.getEquivalenceRatio(sfuel, sox, ThermoBasisType::molar), 1.0, 1e-4);
    EXPECT_NEAR(gas.getEquivalenceRatio(), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::molar, "Bilger"), 0.05516607283, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::molar, "C"), 0.05516607283, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::molar, "O"), 0.05516607283, 1e-4);
    EXPECT_NEAR(gas.getStoichAirFuelRatio(sfuel, sox, ThermoBasisType::molar), 9.52380952380, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), 0.05516607283, 1e-4);

    gas.setEquivalenceRatio(1.0, sfuel, sox, ThermoBasisType::mass);
    gas.setMixtureFraction(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::mass, "Bilger"), sfuel, sox, ThermoBasisType::mass);
    EXPECT_NEAR(gas.getEquivalenceRatio(sfuel, sox, ThermoBasisType::mass), 1.0, 1e-4);
    EXPECT_NEAR(gas.getEquivalenceRatio(), 1.0, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::mass, "Bilger"), 0.0500096579, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::mass, "C"), 0.0500096579, 1e-4);
    EXPECT_NEAR(gas.getMixtureFraction(sfuel, sox, ThermoBasisType::mass, "O"), 0.0500096579, 1e-4);
    EXPECT_NEAR(gas.getStoichAirFuelRatio(sfuel, sox, ThermoBasisType::mass), 10.593805138247204, 1e-4);
    EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), 0.0500096579, 1e-4);
}

}
