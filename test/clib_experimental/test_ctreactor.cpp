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
    int sol = sol_newSolution("gri30.yaml", "gri30", "none");
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

vector<double> T_ctreactor = {
    1050.000, 1050.064, 1050.197, 1050.369, 1050.593, 1050.881, 1051.253, 1051.736,
    1052.370, 1053.216, 1054.372, 1056.007, 1058.448, 1062.431, 1070.141, 1094.331,
    2894.921, 2894.921, 2894.921, 2894.921, 2894.921};

TEST(ctreactor, reactor_simple)
{
    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";

    int sol = sol_newSolution("gri30.yaml", "gri30", "none");
    int thermo = sol_thermo(sol);
    thermo_setMoleFractionsByName(thermo, X.c_str());
    thermo_setTemperature(thermo, T);
    thermo_setPressure(thermo, P);

    int reactor = reactor_new("IdealGasReactor", sol, "test");
    vector<int> reactors{reactor};
    int net = reactornet_new(1, reactors.data());
    ASSERT_EQ(net, 0);

    double t = 0.0;
    int count = 0;
    while (t < 0.1) {
        double Tref = T_ctreactor[count];
        ASSERT_NEAR(reactor_temperature(reactor), Tref, 1e-2);
        t = reactornet_time(net) + 5e-3;
        int ret = reactornet_advance(net, t);
        ASSERT_EQ(ret, 0);
        count++;
    }
}
