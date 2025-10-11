#include <gtest/gtest.h>

#include "cantera/core.h"
#include "cantera/zerodim.h"

#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/ctreactor.h"
#include "cantera_clib/ctreactornet.h"

using namespace Cantera;

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
