#include "gtest/gtest.h"
#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"

#include <cmath>


namespace Cantera
{

class PlasmaPhase_Test : public testing::Test
{
public:
    PlasmaPhase_Test() {
        auto thermo = newThermo("example_data/oxygen-plasma-itikawa.yaml");
        test_phase = std::dynamic_pointer_cast<PlasmaPhase>(thermo);
    }

    shared_ptr<PlasmaPhase> test_phase;
};

TEST_F(PlasmaPhase_Test, setState_TP)
{
    test_phase->setState_TP(1000, 1e5);
    EXPECT_DOUBLE_EQ(test_phase->temperature(), 1000);
    // Default electron temperature is set to the value corresponding to a
    // mean energy of 5 eV (Te=2/3 * kb * T).
    EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 38681.7270718336);
    EXPECT_DOUBLE_EQ(test_phase->pressure(), 1e5);

    // Setting gas temperature first will not work.
    // test_phase->setState_TP(2000, 2e5);
    // test_phase->setElectronTemperature(2001);
    // EXPECT_DOUBLE_EQ(test_phase->temperature(), 2000);
    // EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 2001);
    // EXPECT_DOUBLE_EQ(test_phase->pressure(), 2e5);

    // Test setting electron temperature first.
    test_phase->setElectronTemperature(3000);
    test_phase->setState_TP(2002, 3e5);
    EXPECT_DOUBLE_EQ(test_phase->temperature(), 2002);
    EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 3000);
    EXPECT_DOUBLE_EQ(test_phase->pressure(), 3e5);
}

TEST_F(PlasmaPhase_Test, setState_TD)
{
    test_phase->setState_TD(1000, 1e-2);
    EXPECT_DOUBLE_EQ(test_phase->temperature(), 1000);
    // Default electron temperature is set to the value corresponding to a
    // mean energy of 5 eV (Te=2/3 * kb * T).
    EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 38681.7270718336);
    EXPECT_DOUBLE_EQ(test_phase->density(), 1e-2);

    // Setting gas temperature first will not work.
    // test_phase->setState_TD(2000, 2e-2);
    // test_phase->setElectronTemperature(2001);
    // EXPECT_DOUBLE_EQ(test_phase->temperature(), 2000);
    // EXPECT_DOUBLE_EQ(test_phase->density(), 2e-2);
    // EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 2001);

    // Test setting electron temperature first.
    test_phase->setElectronTemperature(3000);
    test_phase->setState_TD(2002, 3e-2);
    EXPECT_DOUBLE_EQ(test_phase->temperature(), 2002);
    EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 3000);
    EXPECT_DOUBLE_EQ(test_phase->density(), 3e-2);
}

TEST_F(PlasmaPhase_Test, setState)
{
    AnyMap state;
    state["gas-temperature"] = 1234;
    state["P"] = "5 bar";
    state["X"] = "O2: 0.8, O2+: 0.1, E: 0.1";
    test_phase->setElectronTemperature(1234);
    test_phase->setState(state);

    EXPECT_DOUBLE_EQ(test_phase->temperature(), 1234);
    EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 1234);
    EXPECT_DOUBLE_EQ(test_phase->pressure(), 5e5);
    EXPECT_DOUBLE_EQ(test_phase->moleFraction("O2"), 0.8);
    EXPECT_DOUBLE_EQ(test_phase->moleFraction("O2+"), 0.1);
    EXPECT_DOUBLE_EQ(test_phase->moleFraction("E"), 0.1);

    AnyMap state2;
    state2["T"] = 1234;
    state2["Te"] = 4321;
    state2["D"] = 1e-2;
    state2["Y"] = "O2: 0.8, O2+: 0.1, E: 0.1";
    test_phase->setState(state2);

    EXPECT_DOUBLE_EQ(test_phase->temperature(), 1234);
    EXPECT_DOUBLE_EQ(test_phase->electronTemperature(), 4321);
    EXPECT_DOUBLE_EQ(test_phase->density(), 1e-2);
    EXPECT_DOUBLE_EQ(test_phase->massFraction("O2"), 0.8);
    EXPECT_DOUBLE_EQ(test_phase->massFraction("O2+"), 0.1);
    EXPECT_DOUBLE_EQ(test_phase->massFraction("E"), 0.1);
}

};
