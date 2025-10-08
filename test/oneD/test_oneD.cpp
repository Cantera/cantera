#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/onedim.h"
#include "cantera/oneD/DomainFactory.h"
#include "cantera/oneD/IonFlow.h"

using namespace Cantera;

// This test is an exact equivalent of a clib test
// (clib::test_ctonedim.cpp::ctonedim::freeflame_from_parts)
TEST(onedim, freeflame)
{
    auto sol = newSolution("h2o2.yaml", "ohmech", "mixture-averaged");
    auto gas = sol->thermo();
    size_t nsp = gas->nSpecies();

    // reactants
    double uin = .3;
    double T = 300;
    double P = 101325;
    string X = "H2:0.65, O2:0.5, AR:2";
    gas->setState_TPX(T, P, X);
    double rho_in = gas->density();
    vector<double> yin(nsp);
    gas->getMassFractions(&yin[0]);

    // product estimate
    gas->equilibrate("HP");
    vector<double> yout(nsp);
    gas->getMassFractions(&yout[0]);
    double rho_out = gas->density();
    double Tad = gas->temperature();

    // flow
    auto flow = newFlow1D("free-flow", sol, "flow");

    // grid
    int nz = 21;
    double lz = 0.02;
    flow->setupUniformGrid(nz, lz);

    // inlet
    auto inlet = newBoundary1D("inlet", sol);
    inlet->setMoleFractions(X);
    inlet->setMdot(uin * rho_in);
    inlet->setTemperature(T);

    // outlet
    auto outlet = newBoundary1D("outlet", sol);
    double uout = inlet->mdot() / rho_out;

    // set up simulation
    vector<shared_ptr<Domain1D>> domains { inlet, flow, outlet };
    auto flame = newSim1D(domains);
    int dom = static_cast<int>(flame->domainIndex("flow"));
    ASSERT_EQ(dom, 1);

    // set up initial guess
    vector<double> locs{0.0, 0.3, 0.7, 1.0};
    vector<double> value{uin, uin, uout, uout};
    flow->setProfile("velocity", locs, value);
    value = {T, T, Tad, Tad};
    flow->setProfile("T", locs, value);
    for (size_t i = 0; i < nsp; i++) {
        value = {yin[i], yin[i], yout[i], yout[i]};
        flow->setProfile(gas->speciesName(i), locs, value);
    }

    // simulation settings
    double ratio = 15.0;
    double slope = 0.3;
    double curve = 0.5;
    flame->setRefineCriteria(dom, ratio, slope, curve);
    flame->setFixedTemperature(0.85 * T + .15 * Tad);

    // solve
    flow->solveEnergyEqn();
    bool refine_grid = false;
    int loglevel = 0;
    flame->solve(loglevel, refine_grid);
    flame->save("gtest-freeflame.yaml", "cpp", "Solution from C++ interface", true);
    if (usesHDF5()) {
        flame->save("gtest-freeflame.h5", "cpp", "Solution from C++ interface", true);
    }

    ASSERT_EQ(flow->nPoints(), static_cast<size_t>(nz + 1));
    double Tprev = flow->value("T", 0);
    for (size_t n = 0; n < flow->nPoints(); n++) {
        T = flow->value("T", n);
        ASSERT_GE(T, Tprev);
        Tprev = T;
    }
}

TEST(onedim, flame_types)
{
    auto sol = newSolution("h2o2.yaml", "ohmech", "mixture-averaged");

    auto free = newDomain<Flow1D>("free-flow", sol, "flow");
    ASSERT_EQ(free->type(), "free-flow");
    auto symm = newDomain<Flow1D>("axisymmetric-flow", sol, "flow");
    ASSERT_EQ(symm->type(), "axisymmetric-flow");
    auto burner = newDomain<Flow1D>("unstrained-flow", sol, "flow");
    ASSERT_EQ(burner->type(), "unstrained-flow");
}

TEST(onedim, ion_flame_types)
{
    auto sol = newSolution("ch4_ion.yaml");
    ASSERT_EQ(sol->transport()->transportModel(), "ionized-gas");

    auto free = newDomain<IonFlow>("free-flow", sol, "flow");
    ASSERT_EQ(free->type(), "free-ion-flow");
    auto symm = newDomain<IonFlow>("axisymmetric-flow", sol, "flow");
    ASSERT_EQ(symm->type(), "axisymmetric-ion-flow");
    auto burner = newDomain<IonFlow>("unstrained-flow", sol, "flow");
    ASSERT_EQ(burner->type(), "unstrained-ion-flow");
}

int main(int argc, char** argv)
{
    printf("Running main() from test_oneD.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    make_deprecation_warnings_fatal();
    printStackTraceOnSegfault();
    CanteraError::setStackTraceDepth(20);
    int result = RUN_ALL_TESTS();
    appdelete();
    return result;
}
