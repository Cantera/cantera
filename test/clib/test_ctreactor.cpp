#include <gtest/gtest.h>

#include "cantera/core.h"
#include "cantera/zerodim.h"

#include "cantera/clib/ct.h"
#include "cantera/clib/ctreactor.h"

using namespace Cantera;

TEST(ctreactor, reactor_basic)
{
    double T = 1050;
    double P = 5 * 101325;
    string Y = "H2:1.0, O2:3.0";

    int sol = soln_newSolution("h2o2.yaml", "ohmech", "none");
    int thermo = soln_thermo(sol);
    thermo_setMassFractionsByName(thermo, Y.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);
    int h2 = thermo_speciesIndex(thermo, "H2");

    int reactor = reactor_new("IdealGasReactor", sol, "my-reactor");
    ASSERT_EQ(reactor, 0);

    int buflen = reactor_name(reactor, 0, 0);
    vector<char> buf(buflen);
    reactor_name(reactor, buflen, buf.data());
    string str(buf.data());
    ASSERT_EQ(str, "my-reactor");

    int ret = reactor_setName(reactor, "spam");
    ASSERT_EQ(ret, 0);
    buflen = reactor_name(reactor, 0, 0);
    buf.resize(buflen);
    reactor_name(reactor, buflen, buf.data());
    str = buf.data();
    ASSERT_EQ(str, "spam");

    buflen = reactor_type(reactor, 0, 0);
    buf.resize(buflen);
    reactor_type(reactor, buflen, buf.data());
    str = buf.data();
    ASSERT_EQ(str, "IdealGasReactor");

    int reservoir = reactor_new("Reservoir", sol, "reservoir");
    int wall = connector_new("Wall", reservoir, reactor, "my-wall");

    buflen = connector_name(wall, 0, 0);
    buf.resize(buflen);
    connector_name(wall, buflen, buf.data());
    str = buf.data();
    ASSERT_EQ(str, "my-wall");

    ret = connector_setName(wall, "eggs");
    ASSERT_EQ(ret, 0);
    buflen = connector_name(wall, 0, 0);
    buf.resize(buflen);
    connector_name(wall, buflen, buf.data());
    str = buf.data();
    ASSERT_EQ(str, "eggs");

    buflen = connector_type(wall, 0, 0);
    buf.resize(buflen);
    connector_type(wall, buflen, buf.data());
    str = buf.data();
    ASSERT_EQ(str, "Wall");

    reactor_setInitialVolume(reactor, 1.234);
    ASSERT_DOUBLE_EQ(reactor_volume(reactor), 1.234);

    ASSERT_DOUBLE_EQ(reactor_pressure(reactor), P);
    ASSERT_DOUBLE_EQ(reactor_temperature(reactor), T);
    ASSERT_DOUBLE_EQ(reactor_massFraction(reactor, h2), .25);
    ASSERT_DOUBLE_EQ(reactor_density(reactor), thermo_density(thermo));
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

    int sol = soln_newSolution("gri30.yaml", "gri30", "none");
    int thermo = soln_thermo(sol);
    thermo_setMoleFractionsByName(thermo, X.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);

    int reactor = reactor_new("IdealGasReactor", sol, "test");
    int net = reactornet_new();
    int ret = reactornet_addreactor(net, reactor);
    ASSERT_EQ(ret, 0);

    double t = 0.0;
    int count = 0;
    while (t < 0.1) {
        double Tref = T_ctreactor[count];
        ASSERT_NEAR(reactor_temperature(reactor), Tref, 1e-2);
        t = reactornet_time(net) + 5e-3;
        ret = reactornet_advance(net, t);
        ASSERT_EQ(ret, 0);
        count++;
    }
}
