#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "cantera/core.h"
#include "cantera/zerodim.h"

#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/ctreactor.h"
#include "cantera_clib/ctreactornet.h"

using namespace Cantera;

namespace {

//! Text of the most recent %Cantera error
string lastError()
{
    int32_t buflen = ct_getCanteraError(0, 0);
    vector<char> buf(buflen);
    ct_getCanteraError(buflen, buf.data());
    return string(buf.data());
}

}

TEST(ctreactor, reactor_soln)
{
    int32_t sol = sol_newSolution("gri30.yaml", "gri30", "none");
    int32_t reactor = reactor_new("IdealGasReactor", sol, 1, "test");
    ASSERT_GE(reactor, 0);

    int32_t ret = reactor_setName(reactor, "spam");
    ASSERT_EQ(ret, 0);
    int32_t buflen = reactor_name(reactor, 0, 0);
    vector<char> buf(buflen);
    reactor_name(reactor, buflen, buf.data());
    string rName(buf.data());
    ASSERT_EQ(rName, "spam");
}

vector<double> T_ctreactor = {
    1050.000, 1050.064, 1050.197, 1050.369, 1050.593, 1050.881, 1051.253, 1051.736,
    1052.370, 1053.216, 1054.372, 1056.007, 1058.448, 1062.431, 1070.141, 1094.331,
    2894.921, 2894.921, 2894.921, 2894.921, 2894.921};

TEST(ctreactor, reactor_simple)
{
    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    int32_t sol = sol_newSolution("gri30.yaml", "gri30", "none");
    int32_t thermo = sol_thermo(sol);
    thermo_setMoleFractionsByName(thermo, X.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);

    int32_t reactor = reactor_new("IdealGasReactor", sol, 0, "test");
    vector<int32_t> reactors{reactor};
    int32_t net = reactornet_new(1, reactors.data());
    ASSERT_GE(net, 0);

    double t = 0.0;
    int32_t count = 0;
    while (t < 0.1) {
        double Tref = T_ctreactor[count];
        ASSERT_NEAR(reactor_temperature(reactor), Tref, 1e-2);
        t = reactornet_time(net) + 5e-3;
        int32_t ret = reactornet_advance(net, t);
        ASSERT_EQ(ret, 0);
        count++;
    }
    // Reactor contents should be the same Solution & ThermoPhase objects
    int32_t phase = reactor_phase(reactor);
    int32_t phase_thermo = sol_thermo(phase);
    EXPECT_DOUBLE_EQ(thermo_temperature(phase_thermo), thermo_temperature(thermo));
    thermo_setTemperature(phase_thermo, 345.0);
    EXPECT_DOUBLE_EQ(thermo_temperature(thermo), 345.0);
}

TEST(ctreactor, reactor_clone)
{
    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    int32_t sol = sol_newSolution("gri30.yaml", "gri30", "none");
    int32_t thermo = sol_thermo(sol);
    thermo_setMoleFractionsByName(thermo, X.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);

    int32_t reactor = reactor_new("IdealGasReactor", sol, 1, "test");
    vector<int32_t> reactors{reactor};
    int32_t net = reactornet_new(1, reactors.data());
    ASSERT_GE(net, 0);

    double t = 0.0;
    int32_t count = 0;
    while (t < 0.1) {
        double Tref = T_ctreactor[count];
        ASSERT_NEAR(reactor_temperature(reactor), Tref, 1e-2);
        t = reactornet_time(net) + 5e-3;
        int32_t ret = reactornet_advance(net, t);
        ASSERT_EQ(ret, 0);
        count++;
    }
    // Reactor contents should be independent objects with a distinct state, and the
    // original thermo object should be unmodified
    int32_t phase = reactor_phase(reactor);
    ASSERT_GE(phase, 0);
    int32_t phase_thermo = sol_thermo(phase);
    ASSERT_GE(phase_thermo, 0);
    EXPECT_DOUBLE_EQ(thermo_temperature(thermo), T);
    EXPECT_GT(thermo_temperature(phase_thermo), T);
}

TEST(ctreactor, surface)
{
    vector<double> Tref = {
        1050.000, 1050.517, 1051.107, 1051.744, 1052.443, 1053.220, 1054.102, 1055.124,
        1056.344, 1057.849, 1059.791, 1062.453, 1066.443, 1073.424, 1091.018, 2918.287,
        2918.287, 2918.287, 2918.287, 2918.287};

    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    int32_t surf = sol_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);
    int32_t gas = sol_adjacent(surf, 0);
    int32_t gas_thermo = sol_thermo(gas);
    thermo_setMoleFractionsByName(gas_thermo, X.c_str());
    thermo_setTemperature(gas_thermo, T);
    thermo_setPressure(gas_thermo, P);

    int32_t reactor = reactor_new("IdealGasReactor", gas, 1, "test");
    vector<int32_t> adjacent = {reactor};
    int32_t rsurf = reactor_newSurface(surf, 1, adjacent.data(), 1, "surf");
    ASSERT_GE(rsurf, 0);
    int32_t net = reactornet_new(1, adjacent.data());
    ASSERT_GE(net, 0);

    double t = 0.0;
    int count = 0;
    while (t < 0.1) {
        ASSERT_NEAR(reactor_temperature(reactor), Tref[count], 1e-2);
        t = reactornet_time(net) + 5e-3;
        int32_t status = reactornet_advance(net, t);
        ASSERT_EQ(status, 0);
        count++;
    }
}

TEST(ctreactor, sensitivity)
{
    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    int32_t sol = sol_newSolution("gri30.yaml", "gri30", "none");
    int32_t thermo = sol_thermo(sol);
    thermo_setMoleFractionsByName(thermo, X.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);

    int32_t reactor = reactor_new("IdealGasReactor", sol, 1, "test");
    vector<int32_t> reactors{reactor};
    int32_t net = reactornet_new(1, reactors.data());
    ASSERT_GE(net, 0);
    ASSERT_EQ(reactor_addSensitivityReaction(reactor, 2), 0);
    ASSERT_EQ(reactornet_setSensitivityTolerances(net, 1e-6, 1e-6), 0);
    ASSERT_EQ(reactornet_advance(net, 1e-3), 0);

    // Locate the index of the temperature component in the global state vector
    int32_t neq = reactornet_neq(net);
    ASSERT_GT(neq, 0);
    int32_t k = -1;
    for (int32_t i = 0; i < neq; i++) {
        int32_t len = reactornet_componentName(net, i, 0, 0);
        vector<char> name(len);
        reactornet_componentName(net, i, len, name.data());
        // ReactorNet component names are prefixed with the reactor name
        if (string(name.data()) == "test: temperature") {
            k = i;
            break;
        }
    }
    ASSERT_GE(k, 0);

    // The index-based and name-based overloads agree for the same component
    double s_index = reactornet_sensitivity(net, k, 0);
    double s_name = reactornet_sensitivityByName(net, "temperature", 0, 0);
    ASSERT_NE(s_index, DERR);
    ASSERT_NE(s_name, DERR);
    EXPECT_DOUBLE_EQ(s_index, s_name);
    EXPECT_TRUE(std::isfinite(s_index));

    // Out-of-range requests report an error
    ASSERT_EQ(reactornet_sensitivity(net, k, 99), DERR);
    ASSERT_EQ(reactornet_sensitivityByName(net, "spam", 0, 0), DERR);
    EXPECT_THAT(lastError(), testing::HasSubstr("spam"));
}

TEST(ctreactor, components)
{
    int32_t sol = sol_newSolution("h2o2.yaml", "", "none");
    int32_t reactor = reactor_new("IdealGasReactor", sol, 1, "first");
    vector<int32_t> reactors{reactor};
    int32_t net = reactornet_new(1, reactors.data());
    ASSERT_GE(net, 0);

    // The state vector of an IdealGasReactor is [mass, volume, temperature,
    // mass fractions...]
    int32_t neq = reactor_neq(reactor);
    ASSERT_EQ(neq, thermo_nSpecies(sol_thermo(sol)) + 3);
    ASSERT_EQ(reactornet_neq(net), neq);

    // Names and indices round-trip for every component of the reactor
    for (int32_t k = 0; k < neq; k++) {
        int32_t len = reactor_componentName(reactor, k, 0, 0);
        ASSERT_GT(len, 0);
        vector<char> name(len);
        reactor_componentName(reactor, k, len, name.data());
        EXPECT_EQ(reactor_componentIndex(reactor, name.data()), k);
    }
    EXPECT_EQ(reactor_componentIndex(reactor, "temperature"), 2);

    // An unknown name and an out-of-range index each report an error
    ASSERT_EQ(reactor_componentIndex(reactor, "spam"), ERR);
    EXPECT_THAT(lastError(), testing::HasSubstr("spam"));

    ASSERT_EQ(reactor_componentName(reactor, neq, 0, 0), -1);
    EXPECT_THAT(lastError(), testing::HasSubstr("IndexError"));
}

TEST(ctreactor, global_component_index)
{
    int32_t sol = sol_newSolution("h2o2.yaml", "", "none");
    int32_t first = reactor_new("IdealGasReactor", sol, 1, "first");
    int32_t second = reactor_new("IdealGasReactor", sol, 1, "second");
    vector<int32_t> reactors{first, second};
    int32_t net = reactornet_new(2, reactors.data());
    ASSERT_GE(net, 0);

    // Components of the second reactor are offset by the length of the first
    int32_t offset = reactor_neq(first);
    ASSERT_GT(offset, 0);
    int32_t local = reactor_componentIndex(second, "OH");
    ASSERT_GE(local, 0);
    int32_t global = reactornet_globalComponentIndex(net, "OH", 1);
    EXPECT_EQ(global, offset + local);

    // globalComponentIndex takes an unprefixed name, while componentName returns
    // the name prefixed with that of the reactor it belongs to
    int32_t len = reactornet_componentName(net, global, 0, 0);
    ASSERT_GT(len, 0);
    vector<char> name(len);
    reactornet_componentName(net, global, len, name.data());
    EXPECT_EQ(string(name.data()), "second: OH");

    // The first reactor is the default
    EXPECT_EQ(reactornet_globalComponentIndex(net, "temperature", 0),
              reactor_componentIndex(first, "temperature"));

    // An unknown name and an out-of-range index each report an error
    ASSERT_EQ(reactornet_globalComponentIndex(net, "spam", 0), ERR);
    EXPECT_THAT(lastError(), testing::HasSubstr("spam"));

    ASSERT_EQ(reactornet_componentName(net, reactornet_neq(net), 0, 0), -1);
    EXPECT_THAT(lastError(), testing::HasSubstr("IndexError"));
}
