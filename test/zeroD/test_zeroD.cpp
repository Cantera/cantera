#include "gtest/gtest.h"
#include "cantera/zerodim.h"
#include "cantera/numerics/eigen_sparse.h"
#include "cantera/numerics/SystemJacobianFactory.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/thermo/PlasmaPhase.h"
#include <numeric>

using namespace Cantera;

// This test is an (almost) exact equivalent of a clib test
// (clib::test_ctreactor.cpp::ctreactor::reactor_simple)
TEST(zerodim, simple)
{
    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    auto sol = newSolution("gri30.yaml", "gri30", "none");
    sol->thermo()->setState_TPX(T, P, X);
    auto reactor = newReactor4("IdealGasReactor", sol, true, "simple");
    ASSERT_EQ(reactor->name(), "simple");
    reactor->initialize();
    ReactorNet network(reactor);
    network.initialize();

    double t = 0.0;
    while (t < 0.1) {
        ASSERT_GE(reactor->temperature(), T);
        t = network.time() + 5e-3;
        network.advance(t);
    }
}

// Test guards preventing segfaults for uninitialized zerodim objects
TEST(zerodim, test_guards)
{
    // Wall with no adjacent reactors
    Wall wall;
    EXPECT_THROW(wall.heatRate(), CanteraError);
    EXPECT_THROW(wall.expansionRate(), CanteraError);

    // FlowDevice with no adjacent reactors
    EXPECT_THROW(FlowDevice().massFlowRate(), CanteraError);
    EXPECT_THROW(MassFlowController().updateMassFlowRate(0.), CanteraError);
    EXPECT_THROW(PressureController().updateMassFlowRate(0.), CanteraError);
    EXPECT_THROW(Valve().updateMassFlowRate(0.), CanteraError);
}

TEST(zerodim, reactor)
{
    auto gas = newSolution("h2o2.yaml", "ohmech", "none");
    auto reactor = newReactor4("Reactor", gas, true, "my-reactor");

    ASSERT_THROW(reactor->componentName(200), IndexError);
    ASSERT_THROW(reactor->componentIndex("spam"), CanteraError);
}

TEST(zerodim, reservoir)
{
    auto gas = newSolution("gri30.yaml", "gri30", "none");
    auto res = newReservoir(gas, true, "my-reservoir");
    ASSERT_EQ(res->type(), "Reservoir");
    ASSERT_EQ(res->name(), "my-reservoir");
}

TEST(zerodim, flowdevice)
{
    auto gas = newSolution("gri30.yaml", "gri30", "none");

    auto node0 = newReactor4("IdealGasReactor", gas, true, "upstream");
    auto node1 = newReactor4("IdealGasReactor", gas, true, "downstream");

    auto valve = newFlowDevice("Valve", node0, node1, "valve");
    ASSERT_EQ(valve->name(), "valve");
    ASSERT_EQ(valve->in().name(), "upstream");
    ASSERT_EQ(valve->out().name(), "downstream");

    ASSERT_EQ(node0->nInlets(), 0U);
    ASSERT_EQ(node0->nOutlets(), 1U);
    ASSERT_EQ(node1->nInlets(), 1U);
    ASSERT_EQ(node1->nOutlets(), 0U);
}

TEST(zerodim, wall)
{
    auto gas = newSolution("gri30.yaml", "gri30", "none");

    auto node0 = newReactor4("IdealGasReactor", gas, true, "left");
    auto node1 = newReactor4("IdealGasReactor", gas, true, "right");

    auto wall = newWall("Wall", node0, node1, "wall");
    ASSERT_EQ(wall->name(), "wall");
    ASSERT_EQ(wall->left().name(), "left");
    ASSERT_EQ(wall->right().name(), "right");

    ASSERT_EQ(node0->nWalls(), 1U);
    ASSERT_EQ(node1->nWalls(), 1U);
}

TEST(zerodim, mole_reactor)
{
    // simplified version of continuous_reactor.py
    auto gas = newSolution("h2o2.yaml", "ohmech", "none");

    auto tank = make_shared<Reservoir>(gas, true, "fuel-air-tank");
    auto exhaust = make_shared<Reservoir>(gas, true, "exhaust");

    auto stirred = make_shared<IdealGasMoleReactor>(gas, true, "stirred-reactor");
    stirred->setEnergyEnabled(false);
    stirred->setInitialVolume(30.5 * 1e-6);

    auto mfc = make_shared<MassFlowController>(tank, stirred, "mass-flow-controller");
    double residenceTime = 2.;
    double mass = stirred->mass();
    mfc->setMassFlowRate(mass/residenceTime);

    auto preg = make_shared<PressureController>(stirred, exhaust, "pressure-regulator");
    preg->setPrimary(mfc);
    preg->setDeviceCoefficient(1e-3);

    auto net = ReactorNet(stirred);
    net.initialize();
}

TEST(zerodim, mole_reactor_2)
{
    // simplified version of continuous_reactor.py
    auto gas = newSolution("h2o2.yaml", "ohmech", "none");

    auto tank = newReservoir(gas, true, "fuel-air-tank");
    auto exhaust = newReservoir(gas, true, "exhaust");

    auto stirred = newReactor4(
        "IdealGasMoleReactor", gas, true, "stirred-reactor");
    stirred->setEnergyEnabled(false);
    stirred->setInitialVolume(30.5 * 1e-6);

    auto mfc = newFlowDevice(
        "MassFlowController", tank, stirred, "mass-flow-controller");
    double residenceTime = 2.;
    double mass = stirred->mass();
    mfc->setMassFlowRate(mass/residenceTime);

    auto preg = newFlowDevice(
        "PressureController", stirred, exhaust, "pressure-regulator");
    preg->setPrimary(mfc);
    preg->setDeviceCoefficient(1e-3);

    vector<shared_ptr<ReactorBase>> reactors{stirred};
    auto net = newReactorNet(reactors);
    net->initialize();
}

// This test ensures that prior reactor initialization of a reactor does
// not affect later integration within a network. This example was
// adapted from test_reactor.py::test_equilibrium_HP.
TEST(zerodim, test_individual_reactor_initialization)
{
    // initial conditions
    double T0 = 1100.0;
    double P0 = 10 * OneAtm;
    double tol = 1e-7;
    string X0 = "H2:1.0, O2:0.5, AR:8.0";
    // reactor solution, phase, and kinetics objects
    shared_ptr<Solution> sol1 = newSolution("h2o2.yaml");
    sol1->thermo()->setState_TPX(T0, P0, X0);
    // set up reactor object
    auto reactor1 = newReactor4("Reactor", sol1, true);
    // initialize reactor prior to integration to ensure no impact
    reactor1->initialize();
    // setup reactor network and integrate
    ReactorNet network(reactor1);
    network.initialize();
    network.advance(1.0);
    // secondary gas for comparison
    shared_ptr<Solution> sol2 = newSolution("h2o2.yaml");
    sol2->thermo()->setState_TPX(T0, P0, X0);
    sol2->thermo()->equilibrate("UV");
    // secondary reactor for comparison
    auto reactor2 = newReactor4("Reactor", sol2, true);
    reactor2->initialize(0.0);
    // get state of reactors
    vector<double> state1(reactor1->neq());
    vector<double> state2(reactor2->neq());
    reactor1->getState(state1);
    reactor2->getState(state2);
    // compare the reactors.
    EXPECT_EQ(reactor1->neq(), reactor2->neq());
    for (size_t i = 0; i < reactor1->neq(); i++) {
        EXPECT_NEAR(state1[i], state2[i], tol);
    }
}

TEST(zerodim, plasma_reactor_energy)
{
    auto sol = newSolution("air-plasma.yaml", "air-plasma-Phelps", "none");
    auto thermo = sol->thermo();
    auto* plasma = dynamic_cast<PlasmaPhase*>(thermo.get());
    ASSERT_NE(plasma, nullptr);
    const double T0 = 300.0;
    const double P0 = OneAtm;
    thermo->setState_TPX(T0, P0, "N2:0.8, O2:0.2, Electron:1E-11");
    // Set up some non-trivial & arbitrary plasma parameters so that
    // jouleHeatingPower() is exercised
    plasma->setElectricField(1e5);
    plasma->updateElectronEnergyDistribution();
    auto reactor = newReactor4("IdealGasReactor", sol, false, "plasma-reactor");
    reactor->setEnergyEnabled(true);
    ReactorNet net(reactor);
    net.initialize();
    const double t_end = 1e-3;
    // The plasma reactor must integrate the corrected (Joule-only) intrinsic
    // heating without error. Without the EEDF correction the electron mobility
    // (and thus the conductivity driving jouleHeatingPower()) is small, so the
    // temperature rise is negligible here; the quantitative temperature-rise
    // regression is validated in the EEDF-correction branch.
    ASSERT_NO_THROW(net.advance(t_end));
    const double T_final = reactor->temperature();
    EXPECT_NEAR(T_final, T0, 1.0);
}

TEST(MoleReactorTestSet, test_mole_reactor_get_state)
{
    // setting up solution object and thermo/kinetics pointers
    double tol = 1e-8;
    auto sol = newSolution("h2o2.yaml");
    sol->thermo()->setState_TPY(1000.0, OneAtm, "H2:0.5, O2:0.5");
    IdealGasConstPressureMoleReactor reactor(sol, true);
    reactor.setInitialVolume(0.5);
    reactor.setEnergyEnabled(false);
    reactor.initialize();
    vector<double> state(reactor.neq());
    // test get state
    const ThermoPhase& thermo = *reactor.phase()->thermo();
    auto imw = thermo.inverseMolecularWeights();
    // prescribed state
    double mass = reactor.volume() * thermo.density();
    size_t H2I = reactor.componentIndex("H2")-1;
    size_t O2I = reactor.componentIndex("O2")-1;
    double O2_Moles = imw[O2I] * 0.5 * mass;
    double  H2_Moles = imw[H2I] * 0.5 * mass;
    // test getState
    reactor.getState(state);
    EXPECT_NEAR(state[reactor.componentIndex("H2")], H2_Moles, tol);
    EXPECT_NEAR(state[reactor.componentIndex("O2")], O2_Moles, tol);
    EXPECT_NEAR(reactor.volume(), 0.5, tol);
    EXPECT_NEAR(reactor.pressure(), OneAtm, tol);
}

class TestReactorNet : public ReactorNet
{
public:
    using ReactorNet::ReactorNet;

    const vector<double>& atolVector() const {
        return m_atol;
    }
};

TEST(MoleReactorTestSet, test_bulk_mole_reactor_default_atol)
{
    auto makeNetwork = [](double volume) {
        auto sol = newSolution("h2o2.yaml", "ohmech", "none");
        sol->thermo()->setState_TPX(1000.0, 2 * OneAtm, "H2:0.5, O2:0.5");
        auto reactor = make_shared<IdealGasMoleReactor>(sol, true);
        reactor->setInitialVolume(volume);
        return make_pair(reactor, make_shared<TestReactorNet>(reactor));
    };

    auto [reactor1, net1] = makeNetwork(1.0e-09);
    auto [reactor2, net2] = makeNetwork(1.0e+04);
    net1->initialize();
    net2->initialize();

    vector<double> state1(reactor1->neq());
    vector<double> state2(reactor2->neq());
    reactor1->getState(state1);
    reactor2->getState(state2);
    size_t sidx = reactor1->speciesOffset();
    double moles1 = std::accumulate(state1.begin() + sidx, state1.end(), 0.0);
    double moles2 = std::accumulate(state2.begin() + sidx, state2.end(), 0.0);

    // Non-species components should have the default atol, while species components
    // should have atol scaled by the total number of moles in the reactor.
    EXPECT_DOUBLE_EQ(net1->atolVector()[0], net1->atol());
    EXPECT_DOUBLE_EQ(net1->atolVector()[1], net1->atol());
    for (size_t k = 0; k < reactor1->phase()->thermo()->nSpecies(); k++) {
        EXPECT_DOUBLE_EQ(net1->atolVector()[sidx + k], net1->atol() * moles1);
        EXPECT_DOUBLE_EQ(net2->atolVector()[sidx + k], net2->atol() * moles2);
        EXPECT_DOUBLE_EQ(net2->atolVector()[sidx + k] / net1->atolVector()[sidx + k],
                         moles2 / moles1);
    }
}

TEST(MoleReactorTestSet, test_mole_reactor_surface_default_atol)
{
    auto makeNetwork = [](double area) {
        auto surf = newSolution("methane_pox_on_pt.yaml", "Pt_surf");
        auto gas = surf->adjacent("gas");
        gas->thermo()->setState_TPX(900.0, OneAtm, "CH4:1.0, O2:2.0, AR:7.52");
        auto reactor = make_shared<IdealGasMoleReactor>(gas, true);
        vector<shared_ptr<ReactorBase>> adjacent{reactor};
        auto rsurf = newReactorSurface("MoleReactorSurface", surf, adjacent, true);
        rsurf->setArea(area);
        return make_tuple(reactor, rsurf, make_shared<TestReactorNet>(reactor));
    };

    auto [reactor1, rsurf1, net1] = makeNetwork(1.0e-07);
    auto [reactor2, rsurf2, net2] = makeNetwork(3.0e+05);
    net1->initialize();
    net2->initialize();

    vector<double> state1(rsurf1->neq());
    vector<double> state2(rsurf2->neq());
    rsurf1->getState(state1);
    rsurf2->getState(state2);
    double moles1 = std::accumulate(state1.begin(), state1.end(), 0.0);
    double moles2 = std::accumulate(state2.begin(), state2.end(), 0.0);

    for (size_t k = 0; k < rsurf1->neq(); k++) {
        size_t i1 = rsurf1->offset() + k;
        size_t i2 = rsurf2->offset() + k;
        EXPECT_DOUBLE_EQ(net1->atolVector()[i1], net1->atol() * moles1);
        EXPECT_DOUBLE_EQ(net2->atolVector()[i2], net2->atol() * moles2);
        EXPECT_DOUBLE_EQ(net2->atolVector()[i2] / net1->atolVector()[i1],
                         moles2 / moles1);
    }
}

TEST(MoleReactorTestSet, test_user_atol_is_not_scaled)
{
    auto surf = newSolution("methane_pox_on_pt.yaml", "Pt_surf");
    auto gas = surf->adjacent("gas");
    gas->thermo()->setState_TPX(900.0, OneAtm, "CH4:1.0, O2:2.0, AR:7.52");
    auto reactor = make_shared<IdealGasMoleReactor>(gas, true);
    reactor->setInitialVolume(4000.0);
    vector<shared_ptr<ReactorBase>> adjacent{reactor};
    auto rsurf = newReactorSurface("MoleReactorSurface", surf, adjacent, true);
    rsurf->setArea(0.0001);
    TestReactorNet net(reactor);
    double atol = 2.0e-20;
    net.setAbsoluteTolerance(atol);
    net.initialize();

    for (const auto& value : net.atolVector()) {
        EXPECT_DOUBLE_EQ(value, atol);
    }

    net.clearAbsoluteTolerance();
    vector<double> reactorState(reactor->neq());
    reactor->getState(reactorState);
    size_t sidx = reactor->speciesOffset();
    double bulkMoles = std::accumulate(reactorState.begin() + sidx,
                                       reactorState.end(), 0.0);
    EXPECT_DOUBLE_EQ(net.atolVector()[sidx], net.atol() * bulkMoles);
}

TEST(MoleReactorTestSet, test_reactor_atol_vector)
{
    auto surf = newSolution("methane_pox_on_pt.yaml", "Pt_surf");
    auto gas = surf->adjacent("gas");
    gas->thermo()->setState_TPX(900.0, OneAtm, "CH4:1.0, O2:2.0, AR:7.52");
    auto reactor = make_shared<IdealGasMoleReactor>(gas, true);
    vector<shared_ptr<ReactorBase>> adjacent{reactor};
    auto rsurf = newReactorSurface("MoleReactorSurface", surf, adjacent, true);

    vector<double> reactorAtol(reactor->neq());
    vector<double> surfaceAtol(rsurf->neq());
    for (size_t i = 0; i < reactorAtol.size(); i++) {
        reactorAtol[i] = (i + 1) * 1.0e-20;
    }
    for (size_t i = 0; i < surfaceAtol.size(); i++) {
        surfaceAtol[i] = (i + 1) * 1.0e-25;
    }
    reactor->setAbsoluteTolerances(reactorAtol);
    rsurf->setAbsoluteTolerances(surfaceAtol);

    TestReactorNet net(reactor);
    net.setAbsoluteTolerance(2.0e-30);
    net.initialize();

    for (size_t i = 0; i < reactorAtol.size(); i++) {
        EXPECT_DOUBLE_EQ(net.atolVector()[reactor->offset() + i], reactorAtol[i]);
    }
    for (size_t i = 0; i < surfaceAtol.size(); i++) {
        EXPECT_DOUBLE_EQ(net.atolVector()[rsurf->offset() + i], surfaceAtol[i]);
    }
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

    int fillfactor = static_cast<int>(testSize) / 2;
    precon.setIlutFillFactor(fillfactor);
    EXPECT_NEAR(precon.ilutFillFactor(), fillfactor, tol);

    double gamma = 1;
    precon.setGamma(gamma);
    EXPECT_NEAR(precon.gamma(), gamma, tol);
    // test setup and getting the matrix
    precon.updatePreconditioner();
    Eigen::SparseMatrix<double> identity(testSize, testSize);
    identity.setIdentity();
    EXPECT_TRUE(precon.matrix().isApprox(identity));
    // test solve
    vector<double> output(testSize, 0);
    vector<double> rhs_vector(testSize, 10);
    precon.solve(rhs_vector, output);
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
    precon.updatePreconditioner();
    EXPECT_TRUE(precon.matrix().isApprox(identity * (thresh * 1.1)));
    // reset and setup then test again
    precon.reset();
    precon.updatePreconditioner();
    EXPECT_TRUE(precon.matrix().isApprox(identity));
}

TEST(AdaptivePreconditionerTests, test_precon_solver_stats)
{
    // setting up solution object and thermo/kinetics pointers
    auto sol = newSolution("h2o2.yaml");
    sol->thermo()->setState_TPY(1000.0, OneAtm, "H2:0.5, O2:0.5");
    auto reactor = newReactor4("IdealGasMoleReactor", sol, true);
    reactor->setInitialVolume(0.5);
    // setup reactor network and integrate
    ReactorNet network(reactor);
    // setup preconditioner
    shared_ptr<SystemJacobian> precon_ptr = newSystemJacobian("Adaptive");
    network.setPreconditioner(precon_ptr);
    EXPECT_THROW(network.step(), CanteraError);
    // take a step
    network.setLinearSolverType("GMRES");
    // get solver stats
    network.step();
    AnyMap stats = network.solverStats();
    EXPECT_GE(stats["jac_evals"].asInt(), 0);
    EXPECT_GE(stats["lin_rhs_evals"].asInt(), 0);
    EXPECT_GE(stats["lin_iters"].asInt(), 0);
    EXPECT_GE(stats["lin_conv_fails"].asInt(), 0);
    EXPECT_GE(stats["prec_evals"].asInt(), 0);
    EXPECT_GE(stats["prec_solves"].asInt(), 0);
    EXPECT_GE(stats["jt_vec_setup_evals"].asInt(), 0);
    EXPECT_GE(stats["jt_vec_prod_evals"].asInt(), 0);
    EXPECT_GE(stats["nonlinear_iters"].asInt(), 0);
    EXPECT_GE(stats["nonlinear_conv_fails"].asInt(), 0);
}

TEST(AdaptivePreconditionerTests, multi_reactor_valve_pressure_coupling)
{
    auto gas = newSolution("h2o2.yaml");

    gas->thermo()->setState_TPX(1000.0, 2.0 * OneAtm, "H2:2.0, O2:1.0, AR:8.0");
    auto upstream = newReactor4("IdealGasMoleReactor", gas, true, "upstream");
    upstream->setEnergyEnabled(false);
    upstream->setInitialVolume(0.5);

    gas->thermo()->setState_TPX(900.0, OneAtm, "H2:1.0, O2:1.0, AR:8.0");
    auto downstream = newReactor4(
        "IdealGasConstPressureMoleReactor", gas, true, "downstream");
    downstream->setEnergyEnabled(false);
    downstream->setInitialVolume(0.5);

    double kValve = 1e-5;
    auto valve = make_shared<Valve>(upstream, downstream);
    valve->setDeviceCoefficient(kValve);

    vector<shared_ptr<ReactorBase>> reactors{upstream, downstream};
    ReactorNet network(reactors);
    Eigen::SparseMatrix<double> jac = network.jacobian();

    size_t h2 = downstream->phase()->thermo()->speciesIndex("H2");
    size_t row = upstream->neq() + downstream->componentIndex("H2");
    double imwH2 = downstream->phase()->thermo()->inverseMolecularWeights()[h2];
    double Yh2 = upstream->phase()->thermo()->massFraction(h2);
    double coeff = imwH2 * Yh2 * kValve;
    double dPdT = upstream->pressure() / upstream->temperature();
    double dPdV = -upstream->pressure() / upstream->volume();

    EXPECT_NEAR(jac.coeff(row, upstream->componentIndex("temperature")),
                coeff * dPdT, 1e-10 * std::abs(coeff * dPdT));
    EXPECT_NEAR(jac.coeff(row, upstream->componentIndex("volume")),
                coeff * dPdV, 1e-10 * std::abs(coeff * dPdV));
}

TEST(AdaptivePreconditionerTests, connector_composition_coupling_flag)
{
    auto gas = newSolution("h2o2.yaml");

    gas->thermo()->setState_TPX(1000.0, OneAtm, "H2:2.0, O2:1.0, AR:8.0");
    auto upstream = newReactor4("IdealGasMoleReactor", gas, true, "upstream");
    upstream->setEnergyEnabled(false);
    upstream->setInitialVolume(0.5);

    gas->thermo()->setState_TPX(900.0, OneAtm, "H2:1.0, O2:1.0, AR:8.0");
    auto downstream = newReactor4("IdealGasMoleReactor", gas, true, "downstream");
    downstream->setEnergyEnabled(false);
    downstream->setInitialVolume(0.5);

    double mdot = 0.02;
    auto mfc = make_shared<MassFlowController>(upstream, downstream);
    mfc->setMassFlowRate(mdot);

    vector<shared_ptr<ReactorBase>> reactors{upstream, downstream};
    ReactorNet network(reactors);

    size_t h2 = downstream->phase()->thermo()->speciesIndex("H2");
    size_t row = upstream->neq() + downstream->componentIndex("H2");
    size_t col = upstream->componentIndex("H2");

    auto thermo = upstream->phase()->thermo();
    auto mw = thermo->molecularWeights();
    double Yh2 = thermo->massFraction(h2);
    double dYdn = mw[h2] * (1.0 - Yh2) / upstream->mass();
    double expected = downstream->phase()->thermo()->inverseMolecularWeights()[h2]
                      * mdot * dYdn;

    // By default, composition coupling is included
    Eigen::SparseMatrix<double> fullJac = network.jacobian();
    EXPECT_NEAR(fullJac.coeff(row, col), expected, 1e-10 * std::abs(expected));

    // With the skip flag, composition coupling is excluded
    AnyMap settings;
    settings["skip-connector-composition-dependence"] = true;
    network.setDerivativeSettings(settings);
    Eigen::SparseMatrix<double> sparseJac = network.jacobian();
    EXPECT_EQ(sparseJac.coeff(row, col), 0.0);
}

TEST(AdaptivePreconditionerTests, multi_reactor_wall_coupling)
{
    auto gas = newSolution("h2o2.yaml");

    gas->thermo()->setState_TPX(1000.0, 2.0 * OneAtm, "H2:2.0, O2:1.0, AR:8.0");
    auto left = newReactor4("IdealGasMoleReactor", gas, true, "left");
    left->setInitialVolume(0.5);

    gas->thermo()->setState_TPX(900.0, OneAtm, "H2:1.0, O2:1.0, AR:8.0");
    auto right = newReactor4("IdealGasMoleReactor", gas, true, "right");
    right->setInitialVolume(0.5);

    double area = 2.0;
    double heatTransferCoeff = 20.0;
    double expansionCoeff = 1e-9;
    auto wall = make_shared<Wall>(left, right);
    wall->setArea(area);
    wall->setHeatTransferCoeff(heatTransferCoeff);
    wall->setEmissivity(0.0);
    wall->setExpansionRateCoeff(expansionCoeff);

    vector<shared_ptr<ReactorBase>> reactors{left, right};
    ReactorNet network(reactors);
    Eigen::SparseMatrix<double> jac = network.jacobian();

    double totalCv = left->mass() * left->phase()->thermo()->cv_mass();
    size_t rightTemperature = left->neq() + right->componentIndex("temperature");
    double expansionWorkTerm = left->pressure() * area * expansionCoeff
        * right->pressure() / right->temperature() / totalCv;
    EXPECT_NEAR(jac.coeff(left->componentIndex("temperature"), rightTemperature),
                area * heatTransferCoeff / totalCv + expansionWorkTerm,
                1e-10 * (area * heatTransferCoeff / totalCv + expansionWorkTerm));

    double dPdT = left->pressure() / left->temperature();
    double dPdV = -left->pressure() / left->volume();
    EXPECT_NEAR(jac.coeff(left->componentIndex("volume"),
                          left->componentIndex("temperature")),
                area * expansionCoeff * dPdT,
                1e-8 * std::abs(area * expansionCoeff * dPdT));
    EXPECT_NEAR(jac.coeff(left->componentIndex("volume"),
                          left->componentIndex("volume")),
                area * expansionCoeff * dPdV,
                1e-10 * std::abs(area * expansionCoeff * dPdV));
}

TEST(AdaptivePreconditionerTests, const_pressure_outlet_enthalpy_jacobian)
{
    // Constant pressure reactor as upstream with outlet; energy enabled. Verifies that:
    // (1) the reactor outlet contributes cross-reactor species terms, and
    // (2) addEnthalpyJacobian is called on the upstream reactor when the
    //     downstream reactor processes its inlet with energy on.
    auto gas = newSolution("h2o2.yaml");

    gas->thermo()->setState_TPX(1000.0, 2.0 * OneAtm, "H2:2.0, O2:1.0, AR:8.0");
    auto upstream = newReactor4(
        "IdealGasConstPressureMoleReactor", gas, true, "upstream");

    gas->thermo()->setState_TPX(900.0, OneAtm, "H2:1.0, O2:1.0, N2:8.0");
    auto downstream = newReactor4("IdealGasMoleReactor", gas, true, "downstream");
    downstream->setInitialVolume(0.5);

    double kValve = 1e-5;
    auto valve = make_shared<Valve>(upstream, downstream);
    valve->setDeviceCoefficient(kValve);

    vector<shared_ptr<ReactorBase>> reactors{upstream, downstream};
    ReactorNet network(reactors);
    Eigen::SparseMatrix<double> jac = network.jacobian();

    // Constant pressure reactor inlet-enthalpy Jacobian: cross-reactor term
    // Expected: mdot * cp_mass_upstream / TotalCv_downstream
    double deltaP = upstream->pressure() - downstream->pressure();
    double mdot = kValve * deltaP;
    double cp_up = upstream->phase()->thermo()->cp_mass();
    double totalCvDown = downstream->mass() * downstream->phase()->thermo()->cv_mass();
    double expected_dE_dT_up = mdot * cp_up / totalCvDown;
    size_t downstream_energy = upstream->neq() + downstream->componentIndex("temperature");
    size_t upstream_T = upstream->componentIndex("temperature");
    EXPECT_NEAR(jac.coeff(downstream_energy, upstream_T),
                expected_dE_dT_up, 1e-8 * std::abs(expected_dE_dT_up));

    // Reactor outlet species: coupling from downstream pressure to
    // species equation. The downstream reactor's addPressureJacobian
    // contributes to d(dn_H2/dt)/dT_downstream.
    size_t h2 = upstream->phase()->thermo()->speciesIndex("H2");
    double imwH2 = upstream->phase()->thermo()->inverseMolecularWeights()[h2];
    double Yh2 = upstream->phase()->thermo()->massFraction(h2);
    double dPdT_down = downstream->pressure() / downstream->temperature();
    size_t IGCPMR_h2 = upstream->componentIndex("H2");
    size_t downstream_T = upstream->neq() + downstream->componentIndex("temperature");
    EXPECT_NEAR(jac.coeff(IGCPMR_h2, downstream_T),
                imwH2 * Yh2 * kValve * dPdT_down,
                1e-8 * std::abs(imwH2 * Yh2 * kValve * dPdT_down));
}

TEST(AdaptivePreconditionerTests, flow_device_pressure_function_jacobian)
{
    // Valve with a non-trivial pressure function (quadratic). Verifies that
    // pressureFunction_ddP uses the analytic derivative of the pressure function
    // rather than the linear default.
    auto gas = newSolution("h2o2.yaml");

    gas->thermo()->setState_TPX(1000.0, 2.0 * OneAtm, "H2:2.0, O2:1.0, AR:8.0");
    auto upstream = newReactor4("IdealGasMoleReactor", gas, true, "upstream");
    upstream->setEnergyEnabled(false);
    upstream->setInitialVolume(0.5);

    gas->thermo()->setState_TPX(900.0, OneAtm, "H2:1.0, O2:1.0, AR:8.0");
    auto downstream = newReactor4("IdealGasMoleReactor", gas, true, "downstream");
    downstream->setEnergyEnabled(false);
    downstream->setInitialVolume(0.5);

    double kValve = 1e-5;
    auto valve = make_shared<Valve>(upstream, downstream);
    valve->setDeviceCoefficient(kValve);
    // Quadratic pressure function: mdot = kValve * deltaP^2
    // Derivative: dmdot/d(deltaP) = kValve * 2 * deltaP
    auto pfunc = make_shared<Pow1>(2.0);
    valve->setPressureFunction(pfunc);

    vector<shared_ptr<ReactorBase>> reactors{upstream, downstream};
    ReactorNet network(reactors);
    Eigen::SparseMatrix<double> jac = network.jacobian();

    // Cross-reactor term: d(dn_H2/dt_down)/dT_up.
    size_t h2 = upstream->phase()->thermo()->speciesIndex("H2");
    double imwH2 = downstream->phase()->thermo()->inverseMolecularWeights()[h2];
    double Yh2 = upstream->phase()->thermo()->massFraction(h2);
    double deltaP = upstream->pressure() - downstream->pressure();
    double dmdot_ddP = kValve * 2.0 * deltaP;  // pfunc derivative = 2*deltaP
    double dPdT = upstream->pressure() / upstream->temperature();

    size_t row = upstream->neq() + downstream->componentIndex("H2");
    size_t col = upstream->componentIndex("temperature");
    double expected = imwH2 * Yh2 * dmdot_ddP * dPdT;
    EXPECT_NEAR(jac.coeff(row, col), expected, 1e-8 * std::abs(expected));
}

TEST(AdaptivePreconditionerTests, wall_emissivity_jacobian)
{
    // Wall with non-zero emissivity. Verifies that addHeatRateJacobian includes
    // the Stefan-Boltzmann radiation term in addition to convective heat transfer.
    auto gas = newSolution("h2o2.yaml");

    gas->thermo()->setState_TPX(1000.0, 2.0 * OneAtm, "H2:2.0, O2:1.0, AR:8.0");
    auto left = newReactor4("IdealGasMoleReactor", gas, true, "left");
    left->setInitialVolume(0.5);

    gas->thermo()->setState_TPX(900.0, OneAtm, "H2:1.0, O2:1.0, AR:8.0");
    auto right = newReactor4("IdealGasMoleReactor", gas, true, "right");
    right->setInitialVolume(0.5);

    double area = 2.0;
    double emissivity = 0.5;
    auto wall = make_shared<Wall>(left, right);
    wall->setArea(area);
    wall->setHeatTransferCoeff(0.0);
    wall->setEmissivity(emissivity);
    wall->setExpansionRateCoeff(0.0);

    vector<shared_ptr<ReactorBase>> reactors{left, right};
    ReactorNet network(reactors);
    Eigen::SparseMatrix<double> jac = network.jacobian();

    // d(Q)/dT_right = -4 * emiss * area * StefanBoltz * T_right^3
    // For the left reactor: coeff = -1/TotalCv (f = -1 for the wall's left side).
    // Result: jac[left_T, right_T] = (-1/TotalCv_left) * rightCoeff
    //       = 4 * emiss * area * StefanBoltz * T_right^3 / TotalCv_left
    double totalCvLeft = left->mass() * left->phase()->thermo()->cv_mass();
    double T_right = right->phase()->thermo()->temperature();
    double expected = 4.0 * emissivity * area * StefanBoltz
                      * std::pow(T_right, 3) / totalCvLeft;
    size_t left_T = left->componentIndex("temperature");
    size_t right_T = left->neq() + right->componentIndex("temperature");
    EXPECT_NEAR(jac.coeff(left_T, right_T), expected, 1e-8 * std::abs(expected));
}

int main(int argc, char** argv)
{
    printf("Running main() from test_zeroD.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    Cantera::CanteraError::setStackTraceDepth(20);
    Cantera::addDataDirectory("test/data");
    Cantera::addDataDirectory("data");
    printStackTraceOnSegfault();
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
