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

class EquilRatio_MixFrac_Test : public testing::Test
{
public:
    void initSolution() {
        m_sol = newSolution("gri30.yaml", "gri30", "None");
    }

    void set_arbitrary_mixture(ThermoBasis basis) {
        auto& gas = *m_sol->thermo();
        m_fuel.clear();
        m_fuel["CH4"] = 0.2;
        m_fuel["O2"] = 0.02;
        m_fuel["N2"] = 0.1;
        m_fuel["CO"] = 0.05;
        m_fuel["CO2"] = 0.02;
        m_ox.clear();
        m_ox["O2"] = 0.21;
        m_ox["N2"] = 0.79;
        m_ox["CO"] = 0.04;
        m_ox["CH4"] = 0.01;
        m_ox["CO2"] = 0.03;

        if (basis == ThermoBasis::mass) {
            // convert fuel and oxidizer compositions to (non-normalized) mass fractions
            gas.setState_TPX(300, 1e5, m_fuel);
            m_fuel.clear();
            for (size_t i=0; i!=gas.nSpecies(); ++i) {
                m_fuel[gas.speciesName(i)] = gas.massFraction(i)*3;
            }

            gas.setState_TPX(300, 1e5, m_ox);
            m_ox.clear();
            for (size_t i=0; i!=gas.nSpecies(); ++i) {
                m_ox[gas.speciesName(i)] = gas.massFraction(i)*7;
            }
        }
    }

    void test_arbitrary_equilRatio_MixFrac(ThermoBasis basis) {
        auto& gas = *m_sol->thermo();
        if (basis == ThermoBasis::mass) {
            gas.setState_TPY(300.0, 1e5, m_fuel);
        } else {
            gas.setState_TPX(300.0, 1e5, m_fuel);
        }
        double Y_Cf = gas.elementalMassFraction(gas.elementIndex("C"));
        double Y_Of = gas.elementalMassFraction(gas.elementIndex("O"));
        if (basis == ThermoBasis::mass) {
            gas.setState_TPY(300.0, 1e5, m_ox);
        } else {
            gas.setState_TPX(300.0, 1e5, m_ox);
        }
        double Y_Co = gas.elementalMassFraction(gas.elementIndex("C"));
        double Y_Oo = gas.elementalMassFraction(gas.elementIndex("O"));

        gas.setEquivalenceRatio(1.3, m_fuel, m_ox, basis);
        double T = gas.temperature();

        // set mixture to burnt state to make sure that equivalence ratio and
        // mixture fraction are independent of reaction progress
        gas.equilibrate("HP");
        test_mixture_results(T, basis, 1.3, 1.1726068608195617, 0.13415725911057605,
                (gas.elementalMassFraction(gas.elementIndex("C"))-Y_Co)/(Y_Cf-Y_Co),
                (gas.elementalMassFraction(gas.elementIndex("O"))-Y_Oo)/(Y_Of-Y_Oo),
                8.3901204498353561, m_fuel, m_ox);

        gas.setState_TP(300.0,1e5);
        gas.setMixtureFraction(gas.mixtureFraction(m_fuel, m_ox, basis, "Bilger"), m_fuel, m_ox, basis);
        T = gas.temperature();
        gas.equilibrate("HP");
        test_mixture_results(T, basis, 1.3, 1.1726068608195617, 0.13415725911057605,
                (gas.elementalMassFraction(gas.elementIndex("C"))-Y_Co)/(Y_Cf-Y_Co),
                (gas.elementalMassFraction(gas.elementIndex("O"))-Y_Oo)/(Y_Of-Y_Oo),
                8.3901204498353561, m_fuel, m_ox);
    }

    template<typename T>
    void test_mixture_results(double Temp, ThermoBasis basis, double phi, double loc_phi, double mf_Bilger,
                              double mf_C, double mf_O, double AFR_st, const T& fuel, const T& ox) {
        auto& gas = *m_sol->thermo();
        EXPECT_NEAR(gas.equivalenceRatio(fuel, ox, basis), phi, 1e-4);
        EXPECT_NEAR(gas.equivalenceRatio(), loc_phi, 1e-4);
        EXPECT_NEAR(gas.mixtureFraction(fuel, ox, basis, "Bilger"), mf_Bilger, 1e-4);
        EXPECT_NEAR(gas.mixtureFraction(fuel, ox, basis, "C"), mf_C, 1e-4);
        EXPECT_NEAR(gas.mixtureFraction(fuel, ox, basis, "O"), mf_O, 1e-4);
        EXPECT_NEAR(gas.stoichAirFuelRatio(fuel, ox, basis), AFR_st, 1e-4);
        EXPECT_NEAR(gas.pressure(), 1e5, 1e-4);
        EXPECT_NEAR(Temp, 300.0, 1e-4);
    }

    void test_pure_mixture(ThermoBasis basis, bool oxidizer, double phi, double mf) {
        auto& gas = *m_sol->thermo();
        vector_fp v_ox(gas.nSpecies());
        vector_fp v_fuel(gas.nSpecies());
        v_ox[gas.speciesIndex("O2")] = 21.0;
        v_ox[gas.speciesIndex("N2")] = 79.0;
        v_fuel[gas.speciesIndex("CH4")] = 10.0;

        if (oxidizer) {
            gas.setState_TPX(300.0, 1e5, v_ox.data());
            EXPECT_NEAR(gas.equivalenceRatio(v_fuel.data(), v_ox.data(), basis), phi, 1e-4);
            EXPECT_NEAR(gas.equivalenceRatio(), 0.0, 1e-4);
        } else {
            gas.setState_TPX(300.0, 1e5, v_fuel.data());
            ASSERT_EQ(gas.equivalenceRatio(v_fuel.data(), v_ox.data(), basis) > phi, true);
            ASSERT_EQ(gas.equivalenceRatio() > phi, true);
        }
        EXPECT_NEAR(gas.mixtureFraction(v_fuel.data(), v_ox.data(), basis, "Bilger"), mf, 1e-4);
        EXPECT_NEAR(gas.mixtureFraction(v_fuel.data(), v_ox.data(), basis, "C"), mf, 1e-4);

        double Ych4 = gas.massFraction(gas.speciesIndex("CH4"));
        gas.setState_TPX(300.0, 1e5, "N2:1");
        gas.setMixtureFraction(mf, v_fuel.data(), v_ox.data(), basis);
        EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
        gas.setState_TPX(300.0, 1e5, "N2:1");
        gas.setEquivalenceRatio(phi, v_fuel.data(), v_ox.data(), basis);
        EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), Ych4, 1e-4);
    }

    void test_stoich_mixture(ThermoBasis basis, double mf, double AFR_st) {
        auto& gas = *m_sol->thermo();
        std::string sfuel = "CH4";
        std::string sox = "O2:21,N2:79";
        gas.setState_TP(300.0, 1e5);
        gas.setEquivalenceRatio(1.0, sfuel, sox, basis);
        gas.setMixtureFraction(gas.mixtureFraction(sfuel, sox, basis, "Bilger"), sfuel, sox, basis);
        test_mixture_results(gas.temperature(), basis, 1.0, 1.0, mf, mf, mf, AFR_st, sfuel, sox);
        EXPECT_NEAR(gas.massFraction(gas.speciesIndex("CH4")), mf, 1e-4);
    }

    shared_ptr<Solution> m_sol;
    compositionMap m_fuel;
    compositionMap m_ox;
};

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_Arbitrary_Mixture_Molar)
{
    initSolution();
    set_arbitrary_mixture(ThermoBasis::molar);
    test_arbitrary_equilRatio_MixFrac(ThermoBasis::molar);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_Arbitrary_Mixture_Mass)
{
    initSolution();
    set_arbitrary_mixture(ThermoBasis::mass);
    test_arbitrary_equilRatio_MixFrac(ThermoBasis::mass);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_PureOx_Molar)
{
    initSolution();
    test_pure_mixture(ThermoBasis::molar, true, 0.0, 0.0);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_PureOx_Mass)
{
    initSolution();
    test_pure_mixture(ThermoBasis::mass, true, 0.0, 0.0);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_PureFuel_Molar)
{
    initSolution();
    test_pure_mixture(ThermoBasis::molar, false, 1e10, 1.0);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_PureFuel_Mass)
{
    initSolution();
    test_pure_mixture(ThermoBasis::mass, false, 1e10, 1.0);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_StoichMix_Molar)
{
    initSolution();
    test_stoich_mixture(ThermoBasis::molar, 0.055166413925195397, 17.126971264726048);
}

TEST_F(EquilRatio_MixFrac_Test, EquilRatio_MixFrac_StoichMix_Mass)
{
    initSolution();
    test_stoich_mixture(ThermoBasis::mass, 0.050011556441079318, 18.995378491732041);
}

}
