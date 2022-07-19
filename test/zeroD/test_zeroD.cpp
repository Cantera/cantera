#include "gtest/gtest.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/zerodim.h"
#include "cantera/base/Interface.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/numerics/PreconditionerFactory.h"
#include "cantera/numerics/AdaptivePreconditioner.h"

using namespace Cantera;

// This test ensures that prior reactor initialization of a reactor does
// not affect later integration within a network. This example was
// adapted from test_reactor.py::test_equilibrium_HP.
TEST(ZeroDim, test_individual_reactor_initialization)
{
    // initial conditions
    double T0 = 1100.0;
    double P0 = 10 * OneAtm;
    double tol = 1e-7;
    std::string X0 = "H2:1.0, O2:0.5, AR:8.0";
    // reactor solution, phase, and kinetics objects
    std::shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(T0, P0, X0);
    // set up reactor object
    Reactor reactor1;
    reactor1.insert(sol1);
    // initialize reactor prior to integration to ensure no impact
    reactor1.initialize();
    // setup reactor network and integrate
    ReactorNet network;
    network.addReactor(reactor1);
    network.initialize();
    network.advance(1.0);
    // secondary gas for comparison
    std::shared_ptr<Solution> sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(T0, P0, X0);
    sol2->thermo()->equilibrate("UV");
    // secondary reactor for comparison
    Reactor reactor2;
    reactor2.insert(sol2);
    reactor2.initialize(0.0);
    // get state of reactors
    vector_fp state1(reactor1.neq());
    vector_fp state2(reactor2.neq());
    reactor1.getState(state1.data());
    reactor2.getState(state2.data());
    // compare the reactors.
    EXPECT_EQ(reactor1.neq(), reactor2.neq());
    for (size_t i = 0; i < reactor1.neq(); i++) {
        EXPECT_NEAR(state1[i], state2[i], tol);
    }
}

TEST(MoleReactorTestSet, test_mole_reactor_get_state)
{
    // setting up solution object and thermo/kinetics pointers
    double tol = 1e-8;
    auto sol = newSolution("h2o2.yaml");
    sol->thermo()->setState_TPY(1000.0, OneAtm, "H2:0.5, O2:0.5");
    IdealGasConstPressureMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(0.5);
    reactor.setEnergy(false);
    reactor.initialize();
    vector_fp state(reactor.neq());
    // test get state
    const ThermoPhase& thermo = reactor.contents();
    const vector_fp& imw = thermo.inverseMolecularWeights();
    // prescribed state
    double mass = reactor.volume() * thermo.density();
    size_t H2I = reactor.componentIndex("H2")-1;
    size_t O2I = reactor.componentIndex("O2")-1;
    double O2_Moles = imw[O2I] * 0.5 * mass;
    double  H2_Moles = imw[H2I] * 0.5 * mass;
    // test getState
    reactor.getState(state.data());
    EXPECT_NEAR(state[reactor.componentIndex("H2")], H2_Moles, tol);
    EXPECT_NEAR(state[reactor.componentIndex("O2")], O2_Moles, tol);
    EXPECT_NEAR(reactor.volume(), 0.5, tol);
    EXPECT_NEAR(reactor.pressure(), OneAtm, tol);
}

TEST(AdaptivePreconditionerTests, test_adaptive_precon_utils)
{
    // setting the tolerance
    double tol = 1e-8;
    size_t testSize = 4;
    AdaptivePreconditioner precon;
    precon.initialize(testSize);
    // test get and set utils
    double droptol = 1e-4;
    precon.setIlutDropTol(droptol);
    EXPECT_NEAR(precon.ilutDropTol(), droptol, tol);

    double fillfactor = testSize/2;
    precon.setIlutFillFactor(fillfactor);
    EXPECT_NEAR(precon.ilutFillFactor(), fillfactor, tol);

    double gamma = 1;
    precon.setGamma(gamma);
    EXPECT_NEAR(precon.gamma(), gamma, tol);
    // test setup and getting the matrix
    precon.setup();
    Eigen::SparseMatrix<double> identity(testSize, testSize);
    identity.setIdentity();
    EXPECT_TRUE(precon.matrix().isApprox(identity));
    // test solve
    vector_fp output(testSize, 0);
    vector_fp rhs_vector(testSize, 10);
    precon.solve(testSize, rhs_vector.data(), output.data());
    for (size_t i = 0; i < testSize; i++) {
        EXPECT_NEAR(rhs_vector[i], output[i], tol);
    }
    // test prune preconditioner and threshold
    double thresh = 0.5;
    precon.setThreshold(thresh);
    EXPECT_NEAR(precon.threshold(), thresh, tol);
    for (size_t i = 0; i < testSize; i++) {
        for (size_t j = 0; j < testSize; j++) {
            precon.setValue(i, j, thresh * 0.9);
        }
    }
    Eigen::MatrixXd testMat(testSize, testSize);
    testMat.setIdentity();
    testMat.fill(thresh * 0.9);
    EXPECT_TRUE(precon.jacobian().isApprox(testMat));
    precon.setup();
    EXPECT_TRUE(precon.matrix().isApprox(identity * (thresh * 1.1)));
    // reset and setup then test again
    precon.reset();
    precon.setup();
    EXPECT_TRUE(precon.matrix().isApprox(identity));
}

TEST(AdaptivePreconditionerTests, test_precon_solver_stats)
{
    // setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("h2o2.yaml");
    sol->thermo()->setState_TPY(1000.0, OneAtm, "H2:0.5, O2:0.5");
    IdealGasMoleReactor reactor;
    reactor.insert(sol);
    reactor.setInitialVolume(0.5);
    // setup reactor network and integrate
    ReactorNet network;
    network.addReactor(reactor);
    // setup preconditioner
    std::shared_ptr<PreconditionerBase> precon_ptr = newPreconditioner("Adaptive");
    network.setPreconditioner(precon_ptr);
    EXPECT_THROW(network.step(), CanteraError);
    // take a step
    network.setLinearSolverType("GMRES");
    // get linear solver stats
    network.step();
    AnyMap linearStats = network.linearSolverStats();
    EXPECT_GE(linearStats["jac_evals"].asInt(), 0);
    EXPECT_GE(linearStats["lin_rhs_evals"].asInt(), 0);
    EXPECT_GE(linearStats["lin_iters"].asInt(), 0);
    EXPECT_GE(linearStats["lin_conv_fails"].asInt(), 0);
    EXPECT_GE(linearStats["prec_evals"].asInt(), 0);
    EXPECT_GE(linearStats["prec_solves"].asInt(), 0);
    EXPECT_GE(linearStats["jt_vec_setup_evals"].asInt(), 0);
    EXPECT_GE(linearStats["jt_vec_prod_evals"].asInt(), 0);
    // nonlinear solver stats
    AnyMap nonlinearStats = network.nonlinearSolverStats();
    EXPECT_GE(nonlinearStats["nonlinear_iters"].asInt(), 0);
    EXPECT_GE(nonlinearStats["nonlinear_conv_fails"].asInt(), 0);
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
