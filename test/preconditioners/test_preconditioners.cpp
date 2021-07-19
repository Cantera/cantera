#include "gtest/gtest.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo.h"
#include "cantera/kinetics.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/zerodim.h"
#include <stdlib.h>
#include <time.h>
#include <queue>

using namespace Cantera;

#if CT_SUNDIALS_VERSION >= 40
    TEST(AdaptivePreconditioning, test_run_sim)
    {
        // Setting up solution object and thermo/kinetics pointers
        auto sol = newSolution("methane_twostep.yaml");
        sol->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor;
        reactor.insert(sol);
        reactor.setInitialVolume(1.0);
        // Creating inlet reservoir object and adding gas
        Reservoir inlet;
        inlet.insert(sol);
        //Creating exhaust reservoir object and adding gas
        Reservoir exhaust;
        exhaust.insert(sol);
        // Creating mass flow controllers
        MassFlowController inletMassFlowController;
        PressureController outletMassFlowController;
        //Connecting reactors
        inletMassFlowController.install(inlet,reactor);
        outletMassFlowController.install(reactor,exhaust);
        outletMassFlowController.setMaster(&inletMassFlowController);
        outletMassFlowController.setPressureCoeff(0.01);
        // Set constant massflow rate
        inletMassFlowController.setMassFlowRate(1.0);
        // Creating reactor network
        ReactorNet network;
        // Create and add preconditioner
        AdaptivePreconditioner precon;
        network.setIntegratorType(&precon,GMRES);
        // network->setVerbose(); //Setting verbose to be true
        network.addReactor(reactor); //Adding reactor to network
        // Setting up simulation
        network.setInitialTime(0.0);
        network.setMaxTimeStep(0.1);
        network.setMaxSteps(10000);
        network.setTolerances(1e-6,1e-6);
        network.setSensitivityTolerances(1e-6,1e-6);
        network.step();
    }

    TEST(AdaptivePreconditioning, test_two_step_mechanism)
    {
        // Constants
        double volume = 1.0;
        double startTime = 0.0;
        size_t reactorStart = 0;
        // Setting up solution object and thermo/kinetics pointers
        auto sol = newSolution("methane_twostep.yaml");
        auto thermo = sol->thermo();
        auto kinetics = sol->kinetics();
        thermo->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0, CO:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor;
        reactor.insert(sol);
        reactor.setInitialVolume(volume);
        reactor.initialize();
        // State produced within CVODES for this example
        std::vector<double> y(reactor.neq(), 0.0);
        std::vector<double> ydot(reactor.neq(), 0.0);
        reactor.getState(y.data());
        // Internal preconditioner
        AdaptivePreconditioner internalPrecon;
        std::vector<size_t> preconDims{reactor.neq(), reactor.neq()};
        internalPrecon.initialize(&preconDims, 1e-15);
        internalPrecon.setReactorStart(reactorStart);
        internalPrecon.reactorLevelSetup(&reactor, reactorStart, startTime, y.data(), ydot.data(), nullptr);
        // Creating external preconditioner for comparison
        AdaptivePreconditioner externalPrecon;
        externalPrecon.initialize(internalPrecon.getDimensions(), 1e-15);
        externalPrecon.setReactorStart(reactorStart);
        // Getting data for manual calculations
        std::vector<double> concs(thermo->nSpecies(), 0.0);
        std::vector<double> kf(kinetics->nReactions(), 0.0);
        std::vector<double> kr(kinetics->nReactions(), 0.0);
        thermo->getConcentrations(concs.data());
        kinetics->getFwdRateConstants(kf.data());
        kinetics->getRevRateConstants(kr.data());
        // Assign concentrations to species
        double O2 = concs[0];
        double CH4 = concs[2];
        double CO = concs[4];
        // Setting elements
        // Mass
        externalPrecon.setElement(0, 0, 0);
        // O2
        externalPrecon.setElement(2, 2, (- 2.25 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.25 * kf[1] * CO * std::pow(O2, -0.5)) / volume); // dO2/dO2
        externalPrecon.setElement(2, 3, 0); // dO2/dH20
        externalPrecon.setElement(2, 4, - 1.5 * kf[0] * std::pow(O2, 1.5) / volume); // dO2/dCH4
        externalPrecon.setElement(2, 5, 0.5 * kr[1]); // dO2/dCO2
        externalPrecon.setElement(2, 6, - 0.5 * kf[1] * std::pow(O2, 0.5)); // dO2/dCO
        // H2O
        externalPrecon.setElement(3, 2, 3 * kf[0] * CH4 * std::pow(O2, 0.5) / volume); // dH2O/dO2
        externalPrecon.setElement(3, 3, 0); // dH2O/dH20
        externalPrecon.setElement(3, 4, 2 * kf[0] * std::pow(O2, 1.5) / volume); // dH2O/dCH4
        externalPrecon.setElement(3, 5, 0); // dH2O/dCO2
        externalPrecon.setElement(3, 6, 0); // dH2O/dCO
        // CH4
        externalPrecon.setElement(4, 2, - 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) / volume); // dCH4/dO2
        externalPrecon.setElement(4, 3, 0); // dCH4/dH20
        externalPrecon.setElement(4, 4, - kf[0] * std::pow(O2, 1.5) / volume); // dCH4/dCH4
        externalPrecon.setElement(4, 5, 0); // dCH4/dCO2
        externalPrecon.setElement(4, 6, 0); // dCH4/dCO
        // CO2
        externalPrecon.setElement(5, 2, (0.5 * kf[1] * CO * std::pow(O2, -0.5)) / volume); // dCO2/dO2
        externalPrecon.setElement(5, 3, 0); // dCO2/dH20
        externalPrecon.setElement(5, 4, 0); // dCO2/dCH4
        externalPrecon.setElement(5, 5, - kr[1]); // dCO2/dCO2
        externalPrecon.setElement(5, 6, kf[1] * std::pow(O2, 0.5)); // dCO2/CO
        //CO
        externalPrecon.setElement(6, 2, 1.5 * kf[0] * CH4 * std::pow(O2, 0.5) - 0.5 * kf[1] * CO * std::pow(O2, -0.5) / volume); // dCO/dO2
        externalPrecon.setElement(6, 3, 0); // dCO/dH20
        externalPrecon.setElement(6, 4, kf[0] * std::pow(O2, 1.5) / volume); // dCO/dCH4
        externalPrecon.setElement(6, 5, kr[1]); // dCO/dCO2
        externalPrecon.setElement(6, 6, - kf[1] * std::pow(O2,0.5)); // dCO/CO
        // Temperature Derivatives
        externalPrecon.TemperatureDerivatives(&reactor, startTime, y.data(), ydot.data(), nullptr);
        // Check that the two are equal
        EXPECT_EQ(externalPrecon == internalPrecon, true);
        // Reset Internal and test acceptPreconditioner
        internalPrecon.reset();
        reactor.acceptPreconditioner(&internalPrecon, reactorStart, startTime, y.data(), ydot.data(), nullptr);
        EXPECT_EQ(externalPrecon == internalPrecon, true);
    }

    TEST(AdaptivePreconditioning, test_one_step_mechanism_network)
    {
        // Constants
        double volume = 1.0;
        double startTime = 0.0;
        double sharedThreshold = 1e-16;
        // Setting up solution object and thermo/kinetics pointers one
        auto sol1 = newSolution("methane_onestep.yaml");
        auto thermo1 = sol1->thermo();
        auto kinetics1 = sol1->kinetics();
        thermo1->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor1;
        reactor1.insert(sol1);
        reactor1.setInitialVolume(volume);
        // Setting up solution object and thermo/kinetics pointers two
        auto sol2 = newSolution("methane_onestep.yaml");
        sol2->thermo()->setState_TPX(300.0, 101325, "CH4:1.0, O2:1.0, H2O:0, CO2:0");
        // Set up reactor object
        IdealGasConstPressureReactor reactor2;
        reactor2.insert(sol2);
        reactor2.setInitialVolume(volume);
        // Network
        ReactorNet network;
        network.addReactor(reactor1);
        network.addReactor(reactor2);
        // Create and add preconditioner
        AdaptivePreconditioner internalPrecon;
        internalPrecon.setThreshold(sharedThreshold);
        network.setIntegratorType(&internalPrecon,GMRES);
        network.initialize();
        // Use this to reset absolute tolerance so test passes
        double atol = -1e10;
        internalPrecon.setAbsoluteTolerance(atol);
        // State produced within CVODES for this example
        std::vector<double> ydot(network.neq(), 0.0);
        std::vector<double> y(network.neq(), 0.0);
        network.getState(y.data());
        network.preconditionerSetup(startTime, y.data(), ydot.data(), nullptr);
        // Creating external preconditioner for comparison
        AdaptivePreconditioner externalPrecon;
        externalPrecon.initialize(internalPrecon.getDimensions(), atol);
        externalPrecon.setThreshold(sharedThreshold);
        externalPrecon.setReactorStart(0);
        // Getting data for manual calculations
        std::vector<double> concs(thermo1->nSpecies(), 0.0);
        std::vector<double> kf(kinetics1->nReactions(), 0.0);
        thermo1->getConcentrations(concs.data());
        kinetics1->getFwdRateConstants(kf.data());
        // Assign concentrations to species
        double O2 = concs[0];
        double CH4 = concs[2];
        // Setting elements
        // Mass
        externalPrecon.setElement(0, 0, 0);
        // O2
        externalPrecon.setElement(2, 2, -4 * kf[0] * O2 * CH4 / volume); // dO2/dO2
        externalPrecon.setElement(2, 3, 0); // dO2/dH20
        externalPrecon.setElement(2, 4, -2 * kf[0] * O2 * O2 / volume); // dO2/dCH4
        externalPrecon.setElement(2, 5, 0); // dO2/dCO2
        // H2O
        externalPrecon.setElement(3, 2, 4 * kf[0]*O2*CH4/volume); // dH2O/dO2
        externalPrecon.setElement(3, 3, 0); // dH2O/dH20
        externalPrecon.setElement(3, 4, 2 * kf[0] * O2 * O2 / volume); // dH2O/dCH4
        externalPrecon.setElement(3, 5, 0); // dH2O/dCO2
        // CH4
        externalPrecon.setElement(4, 2, -2 * kf[0] * O2 * CH4 / volume); // dCH4/dO2
        externalPrecon.setElement(4, 3, 0); // dCH4/dH20
        externalPrecon.setElement(4, 4, -kf[0] * O2 * O2 / volume); // dCH4/dCH4
        externalPrecon.setElement(4, 5, 0); // dCH4/dCO2
        // CO2
        externalPrecon.setElement(5, 2, 2 * kf[0] * O2 * CH4 / volume); // dCO2/dO2
        externalPrecon.setElement(5, 3, 0); // dCO2/dH20
        externalPrecon.setElement(5, 4, kf[0] * O2 * O2 / volume); // dCO2/dCO2
        externalPrecon.setElement(5, 5, 0); // dCO2/dCO2
        externalPrecon.TemperatureDerivatives(&reactor1, startTime, y.data(), ydot.data(), nullptr);
        size_t neq = reactor1.neq();
        for (size_t i = 0; i < neq; i++)
        {
            for (size_t j = 0; j < neq; j++)
            {
                externalPrecon.setElement(i+neq,j+neq,externalPrecon.getElement(i,j));
            }
        }
        // Make into preconditioner as P = (I - gamma * J_bar)
        externalPrecon.transformJacobianToPreconditioner();
        // Check that the two are equal
        EXPECT_EQ(externalPrecon==internalPrecon, true);
    }

    TEST(AdaptivePreconditioning, test_preconditioned_hydrogen_auto_ignition)
    {
        // create an ideal gas mixture that corresponds to GRI-Mech 3.0
        auto sol = newSolution("gri30.yaml", "gri30", "None");
        auto gas = sol->thermo();
        // set the state
        gas->setState_TPX(1001.0, OneAtm, "H2:2.0, O2:1.0, N2:4.0");
        // create a reactor
        IdealGasConstPressureReactor r;
        // 'insert' the gas into the reactor and environment.
        r.insert(sol);
        // create preconditioner
        AdaptivePreconditioner precon;
        // create reactor network and set to use preconditioner
        ReactorNet sim;
        sim.addReactor(r);
        sim.setIntegratorType(&precon,GMRES);
        // main loop
        double dt = 1.e-5; // interval at which output is written
        int nsteps = 100; // number of intervals
        for (int i = 1; i <= nsteps; i++) {
            double tm = i*dt;
            sim.advance(tm);
        }
    }

#endif

// This tests initialization and dimension setting
TEST(AdaptivePreconditioning, test_initialize_set_dimensions)
{
    // Seed random number generator
    std::srand(std::time(NULL));
    // Generate random test dimensions between 5 and 15
    size_t base = 5;
    size_t limit = 15;
    size_t randomNum = std::rand() % (limit-base+1) + base;
    std::vector<size_t> dims{randomNum, randomNum};
    // Create preconditioner object
    AdaptivePreconditioner precon;
    // Initialize matrix
    precon.initialize(&dims, 1e-15);
    // Test get dimensions
    std::vector<size_t> *preconDims = precon.getDimensions();
    // Test that the dimensions are set properly via initialize
    EXPECT_EQ (dims.at(0),preconDims->at(0));
    EXPECT_EQ (dims.at(1),preconDims->at(1));
    // Checking inside of sparse matrix
    Eigen::SparseMatrix<double>* SpMat = precon.getMatrix();
    EXPECT_EQ (dims.at(0),(size_t) SpMat->innerSize());
    EXPECT_EQ (dims.at(1),(size_t) SpMat->outerSize());
    // Test Set Dimensions
    size_t newRandom = std::rand() % (limit-base+1) + base;
    std::vector<size_t> newDims{newRandom, newRandom};
    // Initialize with different dimensions
    precon.setDimensions(&newDims);
    // Test that the dimensions are set properly via setDimensions
    EXPECT_EQ (newDims.at(0),preconDims->at(0));
    EXPECT_EQ (newDims.at(1),preconDims->at(1));
    // Test failure of initialize
    dims[1] -= 1;
    EXPECT_THROW(precon.initialize(&dims, 1e-15), CanteraError);
}

TEST(AdaptivePreconditioning, test_get_set_copy_assignment_compare)
{
    // Seed random number generator
    std::srand(std::time(NULL));
    // Generate random test dimensions between 10 and 50
    size_t base = 10;
    size_t limit = 50;
    size_t randomNum = std::rand() % (limit-base+1) + base;
    std::vector<size_t> dims{randomNum, randomNum};
    // Create preconditioner object
    AdaptivePreconditioner precon;
    // Initialize matrix
    double atol = 1e-15;
    precon.initialize(&dims, atol);
    // Set threshold
    double thresh = (double) base+2;
    precon.setThreshold(thresh);
    // Random values to put in matrix
    std::queue<size_t> values;
    for(size_t i=0; i<base; i++)
    {
        values.push(std::rand() % 100);
        values.push(std::rand() % dims[1]);
        values.push(std::rand() % dims[0]);
    }
    // Check set and get elements
    for(size_t i=0; i<base; i++)
    {
        double currElement = (double) values.front();
        values.pop();
        size_t col = values.front();
        values.pop();
        size_t row = values.front();
        values.pop();
        precon.setElement(row,col,currElement);
        double returnedElement = precon.getElement(row,col);
        if(std::abs(currElement) >= thresh || row==col)
        {
            EXPECT_EQ (returnedElement,currElement);
        }
        else
        {
            EXPECT_EQ (returnedElement,0.0);
        }
    }
    // Create preconditioner object for copy and compare
    AdaptivePreconditioner preconCopy = precon;
    EXPECT_EQ(preconCopy==precon, true);
    // Check tolerance matches
    EXPECT_EQ(precon.getAbsoluteTolerance(), atol);
    // Reset preconditioner then compare again
    precon.reset();
    EXPECT_EQ(preconCopy==precon, false);
    // Call assignment then compare again
    precon = preconCopy;
    EXPECT_EQ(preconCopy==precon,true);
}

TEST(AdaptivePreconditioning, test_preconditioner_base)
{
    // Required variables
    PreconditionerBase precon;
    IdealGasConstPressureReactor pressureReactor;
    Reactor generalReactor;
    // Test error throwing of preconditioner base
    EXPECT_THROW(precon.solve(nullptr, nullptr, nullptr, nullptr, 0), CanteraError);
    EXPECT_THROW(precon.setup(nullptr, nullptr, 0.0, nullptr, nullptr, nullptr), CanteraError);
    EXPECT_THROW(generalReactor.acceptPreconditioner(&precon, 0, 0.0, nullptr, nullptr, nullptr), CanteraError);
    EXPECT_THROW(pressureReactor.acceptPreconditioner(&precon, 0, 0.0, nullptr, nullptr, nullptr), CanteraError);
    EXPECT_EQ(PRECONDITIONER_NOT_SET, precon.getPreconditionerType());
    EXPECT_THROW(precon.initialize(nullptr, 1e-15), CanteraError);
    EXPECT_THROW(precon.reset(), CanteraError);
    EXPECT_THROW(precon.setElement(0, 0, 0.0), CanteraError);
    EXPECT_THROW(precon.getElement(0, 0), CanteraError);
    EXPECT_THROW(precon.setTimeStep(0.0), CanteraError);
    EXPECT_THROW(precon.getTimeStep(), CanteraError);
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
