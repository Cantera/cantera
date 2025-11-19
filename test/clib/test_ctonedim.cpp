#include <gtest/gtest.h>

#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"

#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/ctdomain.h"
#include "cantera_clib/ctonedim.h"

using namespace Cantera;

TEST(ctonedim, freeflow)
{
    ct_resetStorage();

    int32_t sol = sol_newSolution("h2o2.yaml", "ohmech", "default");
    int32_t ph = sol_thermo(sol);
    ASSERT_GE(ph, 0);

    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";
    thermo_setMoleFractionsByName(ph, X.c_str());
    thermo_setTemperature(ph, T);
    thermo_setPressure(ph, P);

    int32_t flow = domain_newFlow1D("free-flow", sol, "");
    ASSERT_GE(flow, 0);
    domain_setID(flow, "flow");
    ASSERT_NEAR(flow_pressure(flow), P, 1e-5);

    int32_t buflen = domain_domainType(flow, 0, 0);
    vector<char> buf(buflen);
    domain_domainType(flow, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "free-flow");
}

TEST(ctonedim, inlet)
{
    int32_t sol = sol_newSolution("h2o2.yaml", "ohmech", "default");
    int32_t inlet = domain_newBoundary1D("inlet", sol, "");
    ASSERT_GE(inlet, 0);

    int32_t buflen = domain_domainType(inlet, 0, 0);
    vector<char> buf(buflen);
    domain_domainType(inlet, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "inlet");
}

TEST(ctonedim, outlet)
{
    int32_t sol = sol_newSolution("h2o2.yaml", "ohmech", "default");
    int32_t outlet = domain_newBoundary1D("outlet", sol, "");
    ASSERT_GE(outlet, 0);

    int32_t buflen = domain_domainType(outlet, 0, 0);
    vector<char> buf(buflen);
    domain_domainType(outlet, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "outlet");
}

TEST(ctonedim, reacting_surface)
{
    int32_t interface = sol_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);
    int32_t surf = domain_newBoundary1D("reacting-surface", interface, "");
    ASSERT_GE(surf, 0);

    int32_t buflen = domain_domainType(surf, 0, 0);
    vector<char> buf(buflen);
    domain_domainType(surf, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "reacting-surface");
}

TEST(ctonedim, catcomb)
{
    int32_t sol = sol_newSolution("ptcombust.yaml", "gas", "default");
    int32_t interface = sol_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);

    // inlet and flow domains
    int32_t inlet = domain_newBoundary1D("inlet", sol, "inlet");
    int32_t flow = domain_newFlow1D("axisymmetric-flow", sol, "flow");
    ASSERT_EQ(flow, inlet+1);

    // reacting surface
    int32_t reac_surf = domain_newBoundary1D("reacting-surface", interface, "surface");
    ASSERT_EQ(reac_surf, flow+1);

    // set up stack
    vector<int32_t> doms{inlet, flow, reac_surf};
    int32_t flame = sim1D_newSim1D(3, doms.data());
    ASSERT_GE(flame, 0);
    int32_t dom = sim1D_domainIndex(flame, "flow");
    ASSERT_GE(dom, 0);
}

TEST(ctonedim, freeflame)
{
    ct_resetStorage();
    auto gas = newThermo("h2o2.yaml", "ohmech");

    int32_t sol = sol_newSolution("h2o2.yaml", "ohmech", "default");
    int32_t ph = sol_thermo(sol);
    size_t nsp = thermo_nSpecies(ph);

    // reactants
    double uin = .3;
    double T = 300;
    double P = 101325;
    string X = "H2:0.65, O2:0.5, AR:2";
    thermo_setMoleFractionsByName(ph, X.c_str());
    thermo_setTemperature(ph, T);
    thermo_setPressure(ph, P);
    double rho_in = thermo_density(ph);
    vector<double> yin(nsp);
    thermo_getMassFractions(ph, nsp, yin.data());

    // product estimate
    int32_t ret = thermo_equilibrate(ph, "HP", "auto", 1e-9, 50000, 1000, 0);
    ASSERT_GE(ret, 0);
    double rho_out = thermo_density(ph);
    double Tad = thermo_temperature(ph);
    vector<double> yout(nsp);
    thermo_getMassFractions(ph, nsp, yout.data());

    // flow
    int32_t flow = domain_newFlow1D("free-flow", sol, "flow");

    // grid
    int32_t nz = 21;
    double lz = 0.02;
    domain_setupUniformGrid(flow, nz, lz, 0.);
    ASSERT_EQ(domain_nPoints(flow), nz);

    // inlet
    int32_t reac = domain_newBoundary1D("inlet", sol, "inlet");
    bdry_setMoleFractionsByName(reac, X.c_str());
    bdry_setMdot(reac, uin * rho_in);
    bdry_setTemperature(reac, T);

    // outlet
    int32_t prod = domain_newBoundary1D("outlet", sol, "outlet");
    double uout = bdry_mdot(reac) / rho_out;

    // set up stack
    vector<int32_t> doms{reac, flow, prod};
    int32_t flame = sim1D_newSim1D(3, doms.data());
    int32_t dom = sim1D_domainIndex(flame, "flow");
    ASSERT_EQ(dom, 1);

    // set up initial guess
    vector<double> locs{0.0, 0.3, 0.7, 1.0};
    vector<double> value{uin, uin, uout, uout};
    domain_setProfile(flow, "velocity", 4, locs.data(), 4, value.data());
    value = {T, T, Tad, Tad};
    domain_setProfile(flow, "T", 4, locs.data(), 4, value.data());
    for (size_t i = 0; i < nsp; i++) {
        value = {yin[i], yin[i], yout[i], yout[i]};
        int32_t buflen = thermo_speciesName(ph, i, 0, 0);
        vector<char> buf(buflen);
        thermo_speciesName(ph, i, buflen, buf.data());
        string name(buf.data());
        ASSERT_EQ(name, gas->speciesName(i));
        domain_setProfile(flow, buf.data(), 4, locs.data(), 4, value.data());
    }

    // simulation settings
    double ratio = 15.0;
    double slope = 0.3;
    double curve = 0.5;
    sim1D_setRefineCriteria(flame, dom, ratio, slope, curve, 0.);
    sim1D_setFixedTemperature(flame, 0.85 * T + .15 * Tad);

    // solve and save
    flow_solveEnergyEqn(flow, 1);
    bool refine_grid = false;
    int32_t loglevel = 0;
    sim1D_solve(flame, loglevel, refine_grid);
    ret = sim1D_save(flame, "gtest-freeflame.yaml", "clib",
                     "Solution from CLib interface", true);
    ASSERT_GE(ret, 0);
    if (usesHDF5()) {
        ret = sim1D_save(flame, "gtest-freeflame.h5", "clib",
                         "Solution from CLib interface", true);
        ASSERT_GE(ret, 0);
    }

    ASSERT_EQ(domain_nPoints(flow), nz + 1);
    ASSERT_EQ(domain_nComponents(flow), 16);
    int32_t comp = domain_componentIndex(flow, "T");
    ASSERT_EQ(comp, 2);
    nz = domain_nPoints(flow);
    vector<double> Tvec(nz);
    domain_values(dom, "T", nz, Tvec.data());
    double Tprev = Tvec[0];
    for (int32_t n = 0; n < domain_nPoints(flow); n++) {
        T = Tvec[n];
        ASSERT_GE(T + 1e-3, Tprev);
        Tprev = T;
    }

    ret = sim1D_restore(flame, "gtest-freeflame.yaml", "clib");
    ASSERT_EQ(ret, 0);
}
