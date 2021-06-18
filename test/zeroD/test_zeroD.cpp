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
    //Reactor solution, phase, and kinetics objects
    std::shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    std::shared_ptr<ThermoPhase> gas1 = sol1->thermo();
    std::shared_ptr<Kinetics> kin1 = sol1->kinetics();
    gas1->setState_TPX(T0, P0, X0);
    //Set up reactor object
    Reactor reactor1;
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
    gas2->equilibrate("UV");
    //Secondary reactor for comparison
    Reactor reactor2;
    reactor2.setKineticsMgr(*kin2);
    reactor2.setThermoMgr(*gas2);
    reactor2.initialize(0);
    //Get state of reactors
    double state1 [reactor1.neq()];
    double state2 [reactor2.neq()];
    reactor1.getState(state1);
    reactor2.getState(state2);
    //Compare the reactors.
    EXPECT_EQ(reactor1.neq(), reactor2.neq());
    double tol = 1e-14;
    EXPECT_NEAR(state1[0], state2[0], tol);
    EXPECT_NEAR(state1[1], state2[1], tol);
    for(size_t i = 3; i < reactor2.neq(); i++)
    {
        EXPECT_NEAR(state1[i], state2[i], tol);
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
