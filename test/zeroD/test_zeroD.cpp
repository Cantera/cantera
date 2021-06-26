#include <cstdio>
#include <string>
#include "gtest/gtest.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zerodim.h"

using namespace Cantera;

// This test ensures that prior reactor initialization of a reactor does
// not affect later integration within a network. This example was
// adapted from test_reactor.py::test_equilibrium_HP.
TEST(ZeroDim, test_individual_reactor_initialization)
{
    // Initial conditions
    double T0 = 1100.0;
    double P0 = 10 * OneAtm;
    double tol = 1e-7;
    std::string X0 = "H2:1.0, O2:0.5, AR:8.0";
    // Reactor solution, phase, and kinetics objects
    std::shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(T0, P0, X0);
    // Set up reactor object
    Reactor reactor1;
    reactor1.insert(sol1);
    // Initialize reactor prior to integration to ensure no impact
    reactor1.initialize();
    // Setup reactor network and integrate
    ReactorNet network;
    network.addReactor(reactor1);
    network.advance(1.0);
    // Secondary gas for comparison
    std::shared_ptr<Solution> sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(T0, P0, X0);
    sol2->thermo()->equilibrate("UV");
    // Secondary reactor for comparison
    Reactor reactor2;
    reactor2.insert(sol2);
    reactor2.initialize(0.0);
    // Get state of reactors
    std::vector<double> state1 (reactor1.neq());
    std::vector<double> state2 (reactor2.neq());
    reactor1.getState(state1.data());
    reactor2.getState(state2.data());
    // Compare the reactors.
    EXPECT_EQ(reactor1.neq(), reactor2.neq());
    for (size_t i = 0; i < reactor1.neq(); i++)
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
