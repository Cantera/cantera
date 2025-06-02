#include <gtest/gtest.h>

#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"

#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/ctdomain.h"
#include "cantera_clib/ctonedim.h"

using namespace Cantera;

TEST(ctonedim3, freeflow)
{
    ct3_resetStorage();

    int sol = sol3_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = sol3_thermo(sol);
    ASSERT_GE(ph, 0);

    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";
    thermo3_setMoleFractionsByName(ph, X.c_str());
    thermo3_setTemperature(ph, T);
    thermo3_setPressure(ph, P);

    int flow = domain3_newFlow1D("free-flow", sol, "");
    ASSERT_GE(flow, 0);
    domain3_setID(flow, "flow");
    ASSERT_NEAR(flow3_pressure(flow), P, 1e-5);

    int buflen = domain3_type(flow, 0, 0);
    vector<char> buf(buflen);
    domain3_type(flow, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "free-flow");
}

TEST(ctonedim3, inlet)
{
    int sol = sol3_newSolution("h2o2.yaml", "ohmech", "default");
    int inlet = domain3_newBoundary1D("inlet", sol, "");
    ASSERT_GE(inlet, 0);

    int buflen = domain3_type(inlet, 0, 0);
    vector<char> buf(buflen);
    domain3_type(inlet, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "inlet");
}

TEST(ctonedim3, outlet)
{
    int sol = sol3_newSolution("h2o2.yaml", "ohmech", "default");
    int outlet = domain3_newBoundary1D("outlet", sol, "");
    ASSERT_GE(outlet, 0);

    int buflen = domain3_type(outlet, 0, 0);
    vector<char> buf(buflen);
    domain3_type(outlet, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "outlet");
}

TEST(ctonedim3, reacting_surface)
{
    int interface = sol3_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);
    int surf = domain3_newBoundary1D("reacting-surface", interface, "");
    ASSERT_GE(surf, 0);

    int buflen = domain3_type(surf, 0, 0);
    vector<char> buf(buflen);
    domain3_type(surf, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "reacting-surface");
}

TEST(ctonedim3, catcomb)
{
    int sol = sol3_newSolution("ptcombust.yaml", "gas", "default");
    int interface = sol3_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);

    // inlet and flow domains
    int inlet = domain3_newBoundary1D("inlet", sol, "inlet");
    int flow = domain3_newFlow1D("axisymmetric-flow", sol, "flow");
    ASSERT_EQ(flow, inlet+1);

    // reacting surface
    int reac_surf = domain3_newBoundary1D("reacting-surface", interface, "surface");
    ASSERT_EQ(reac_surf, flow+1);

    // set up stack
    vector<int> doms{inlet, flow, reac_surf};
    int flame = sim1D3_newSim1D(3, doms.data());
    ASSERT_GE(flame, 0);
    int dom = sim1D3_domainIndex(flame, "flow");
    ASSERT_GE(dom, 0);
}

TEST(ctonedim3, freeflame)
{
    ct3_resetStorage();
    auto gas = newThermo("h2o2.yaml", "ohmech");

    int sol = sol3_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = sol3_thermo(sol);
    size_t nsp = thermo3_nSpecies(ph);

    // reactants
    double uin = .3;
    double T = 300;
    double P = 101325;
    string X = "H2:0.65, O2:0.5, AR:2";
    thermo3_setMoleFractionsByName(ph, X.c_str());
    thermo3_setTemperature(ph, T);
    thermo3_setPressure(ph, P);
    double rho_in = thermo3_density(ph);
    vector<double> yin(nsp);
    thermo3_getMassFractions(ph, nsp, yin.data());

    // product estimate
    int ret = thermo3_equilibrate(ph, "HP", "auto", 1e-9, 50000, 1000, 0);
    ASSERT_GE(ret, 0);
    double rho_out = thermo3_density(ph);
    double Tad = thermo3_temperature(ph);
    vector<double> yout(nsp);
    thermo3_getMassFractions(ph, nsp, yout.data());

    // flow
    int flow = domain3_newFlow1D("free-flow", sol, "flow");

    // grid
    int nz = 21;
    double lz = 0.02;
    domain3_setupUniformGrid(flow, nz, lz, 0.);
    ASSERT_EQ(domain3_nPoints(flow), nz);

    // inlet
    int reac = domain3_newBoundary1D("inlet", sol, "inlet");
    bdry3_setMoleFractions(reac, X.c_str());
    bdry3_setMdot(reac, uin * rho_in);
    bdry3_setTemperature(reac, T);

    // outlet
    int prod = domain3_newBoundary1D("outlet", sol, "outlet");
    double uout = bdry3_mdot(reac) / rho_out;

    // set up stack
    vector<int> doms{reac, flow, prod};
    int flame = sim1D3_newSim1D(3, doms.data());
    int dom = sim1D3_domainIndex(flame, "flow");
    ASSERT_EQ(dom, 1);

    // set up initial guess
    vector<double> locs{0.0, 0.3, 0.7, 1.0};
    vector<double> value{uin, uin, uout, uout};
    sim1D3_setInitialGuess(flame, "velocity", 4, locs.data(), 4, value.data());
    value = {T, T, Tad, Tad};
    sim1D3_setInitialGuess(flame, "T", 4, locs.data(), 4, value.data());
    for (size_t i = 0; i < nsp; i++) {
        value = {yin[i], yin[i], yout[i], yout[i]};
        int buflen = thermo3_getSpeciesName(ph, i, 0, 0);
        vector<char> buf(buflen);
        thermo3_getSpeciesName(ph, i, buflen, buf.data());
        string name(buf.data());
        ASSERT_EQ(name, gas->speciesName(i));
        sim1D3_setInitialGuess(flame, buf.data(), 4, locs.data(), 4, value.data());
    }

    // simulation settings
    double ratio = 15.0;
    double slope = 0.3;
    double curve = 0.5;
    sim1D3_setRefineCriteria(flame, dom, ratio, slope, curve, 0.);
    sim1D3_setFixedTemperature(flame, 0.85 * T + .15 * Tad);

    // solve and save
    flow3_solveEnergyEqn(flow, 1);
    bool refine_grid = false;
    int loglevel = 0;
    sim1D3_solve(flame, loglevel, refine_grid);
    ret = sim1D3_save(flame, "gtest-freeflame.yaml", "clib",
                      "Solution from CLib interface");
    ASSERT_GE(ret, 0);
    if (usesHDF5()) {
        ret = sim1D3_save(flame, "gtest-freeflame.h5", "clib",
                          "Solution from CLib interface");
        ASSERT_GE(ret, 0);
    }

    ASSERT_EQ(domain3_nPoints(flow), nz + 1);
    ASSERT_EQ(domain3_nComponents(flow), 16);
    int comp = domain3_componentIndex(flow, "T");
    ASSERT_EQ(comp, 2);
    double Tprev = sim1D3_value(flame, dom, comp, 0);
    for (int n = 0; n < domain3_nPoints(flow); n++) {
        T = sim1D3_value(flame, dom, comp, n);
        ASSERT_GE(T + 1e-3, Tprev);
        Tprev = T;
    }
}
