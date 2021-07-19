#include "gtest/gtest.h"
#include "cantera/numerics/CVodesIntegrator.h"
#include "cantera/zerodim.h"

using namespace Cantera;

TEST(GeneralNumericsTests, test_integrators_funceval)
{
    // Required Variables
    Integrator integ;
    auto cvodesInteg = newIntegrator("CVODE");
    // Setting up solution to insert in reactor
    auto sol = newSolution("h2o2.yaml");
    // Set up reactor object
    Reactor reactor;
    reactor.insert(sol);
    ReactorNet network;
    network.addReactor(reactor);
    network.initialize();
    // Testing parameters
    EXPECT_THROW(integ.getIntegratorTimeStep(), CanteraError);
    cvodesInteg->initialize(0.0, network);
    EXPECT_EQ(cvodesInteg->getIntegratorTimeStep(), 0);
}

int main(int argc, char** argv)
{
    printf("Running main() from test_preconditioners.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
