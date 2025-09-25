#include <gtest/gtest.h>

#include "cantera/core.h"
#include "cantera/zerodim.h"

#include "cantera/clib/ct.h"
#include "cantera/clib/ctreactor.h"

using namespace Cantera;

TEST(ctreactor, reactor_soln)
{
    int sol = soln_newSolution("gri30.yaml", "gri30", "none");
    int reactor = reactor_new("IdealGasReactor", sol, "test");
    ASSERT_EQ(reactor, 0);

    int ret = reactor_setName(reactor, "spam");
    ASSERT_EQ(ret, 0);
    int buflen = reactor_name(reactor, 0, 0);
    vector<char> buf(buflen);
    reactor_name(reactor, buflen, buf.data());
    string rName(buf.data());
    ASSERT_EQ(rName, "spam");
}


TEST(ctreactor, reactor_simple)
{
    vector<double> Tref = {
        1050.000, 1050.064, 1050.197, 1050.369, 1050.593, 1050.881, 1051.253, 1051.736,
        1052.370, 1053.216, 1054.372, 1056.007, 1058.448, 1062.431, 1070.141, 1094.331,
        2894.921, 2894.921, 2894.921, 2894.921, 2894.921};

    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    int sol = soln_newSolution("gri30.yaml", "gri30", "none");
    int thermo = soln_thermo(sol);
    thermo_setMoleFractionsByName(thermo, X.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);

    int reactor = reactor_new("IdealGasReactor", sol, "test");
    int net = reactornet_new();
    suppress_deprecation_warnings();
    int ret = reactornet_addreactor(net, reactor);
    make_deprecation_warnings_fatal();
    ASSERT_EQ(ret, 0);

    double t = 0.0;
    int count = 0;
    while (t < 0.1) {
        ASSERT_NEAR(reactor_temperature(reactor), Tref[count], 1e-2);
        t = reactornet_time(net) + 5e-3;
        ret = reactornet_advance(net, t);
        ASSERT_EQ(ret, 0);
        count++;
    }
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

    int surf = soln_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);
    int gas = soln_adjacent(surf, 0);
    int gas_thermo = soln_thermo(gas);
    thermo_setMoleFractionsByName(gas_thermo, X.c_str());
    thermo_setTemperature(gas_thermo, T);
    thermo_setPressure(gas_thermo, P);

    int reactor = reactor_new("IdealGasReactor", gas, "test");
    suppress_deprecation_warnings();
    int rsurf = reactor_new("ReactorSurface", surf, "surf");
    ASSERT_GE(rsurf, 0);
    int status = reactorsurface_install(rsurf, reactor);
    ASSERT_EQ(status, 0);
    int net = reactornet_new();

    status = reactornet_addreactor(net, reactor);
    make_deprecation_warnings_fatal();
    ASSERT_EQ(status, 0);

    double t = 0.0;
    int count = 0;
    while (t < 0.1) {
        ASSERT_NEAR(reactor_temperature(reactor), Tref[count], 1e-2);
        t = reactornet_time(net) + 5e-3;
        status = reactornet_advance(net, t);
        ASSERT_EQ(status, 0);
        count++;
    }
}
