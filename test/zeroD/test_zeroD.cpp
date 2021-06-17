#include <cstdio>
#include <string>
#include "gtest/gtest.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zerodim.h"

using namespace Cantera;

//This test ensures that prior reactor initialization of a reactor does not affect later integration within a network. This example was adapted from test_reactor.py::test_equilibrium_HP.
TEST(ZeroDim, test_individual_reactor_initialization)
{
    //Initial conditions
    double T0 = 1100.0;
    double P0 = 1013250;
    std::string X0 = "H2:1.0, O2:0.5, AR:8.0";
    //IdealGasConstPressureReactor solution, phase, and kinetics objects
    std::shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    std::shared_ptr<ThermoPhase> gas1 = sol1->thermo();
    std::shared_ptr<Kinetics> kin1 = sol1->kinetics();
    gas1->setState_TPX(T0, P0, X0);
    //Set up reactor object
    IdealGasConstPressureReactor reactor1;
    reactor1.setKineticsMgr(*kin1);
    reactor1.setThermoMgr(*gas1);
    //initialize reactor at arbitrary value
    reactor1.initialize(0.1);
    //Setup reactor network and integrating
    ReactorNet network;
    network.addReactor(reactor1);
    network.advance(1.0);
    //Secondary gas for comparison
    std::shared_ptr<Solution> sol2 = newSolution("h2o2.yaml");
    std::shared_ptr<ThermoPhase> gas2 = sol2->thermo();
    std::shared_ptr<Kinetics> kin2 = sol2->kinetics();
    gas2->setState_TPX(T0, P0, X0);
    gas2->equilibrate("HP");
    //Compare the two gases.
    double tol = 1e-8; //relative tolerance from python assertNear test
    EXPECT_NEAR(gas1->temperature(), gas2->temperature(), tol);
    EXPECT_NEAR(gas1->density(), gas2->density(), tol);
    EXPECT_NEAR(gas1->pressure(), P0, tol);
    for(size_t i = 0; i < gas1->stateSize()-2; i++)
    {
        EXPECT_NEAR(gas1->massFraction(i), gas2->massFraction(i), tol);
    }
}

int main(int argc, char** argv)
{
    printf("Running main() from test_zeroD.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
