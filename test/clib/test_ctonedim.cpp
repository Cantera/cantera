#include <gtest/gtest.h>

#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"

#include "cantera/clib/ct.h"
#include "cantera/clib/ctonedim.h"

using namespace Cantera;

TEST(ctonedim, freeflow)
{
    ct_resetStorage();

    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = soln_thermo(sol);
    ASSERT_GE(ph, 0);

    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";
    thermo_setMoleFractionsByName(ph, X.c_str());
    thermo_setTemperature(ph, T);
    thermo_setPressure(ph, P);

    int flow = domain_new("free-flow", sol, "");
    ASSERT_GE(flow, 0);
    domain_setID(flow, "flow");
    ASSERT_NEAR(flow1D_pressure(flow), P, 1e-5);

    int buflen = domain_type(flow, 0, 0);
    vector<char> buf(buflen);
    domain_type(flow, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "free-flow");
}

TEST(ctonedim, inlet)
{
    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int inlet = domain_new("inlet", sol, "");
    ASSERT_GE(inlet, 0);

    int buflen = domain_type(inlet, 0, 0);
    vector<char> buf(buflen);
    domain_type(inlet, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "inlet");
}

TEST(ctonedim, outlet)
{
    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int outlet = domain_new("outlet", sol, "");
    ASSERT_GE(outlet, 0);

    int buflen = domain_type(outlet, 0, 0);
    vector<char> buf(buflen);
    domain_type(outlet, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "outlet");
}

TEST(ctonedim, reacting_surface)
{
    int interface = soln_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);
    int surf = domain_new("reacting-surface", interface, "");
    ASSERT_GE(surf, 0);

    int buflen = domain_type(surf, 0, 0);
    vector<char> buf(buflen);
    domain_type(surf, buflen, buf.data());
    string domType(buf.data());
    ASSERT_EQ(domType, "reacting-surface");
}

TEST(ctonedim, catcomb)
{
    int sol = soln_newSolution("ptcombust.yaml", "gas", "default");
    int interface = soln_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);

    // inlet and flow domains
    int inlet = domain_new("inlet", sol, "inlet");
    int flow = domain_new("axisymmetric-flow", sol, "flow");
    ASSERT_EQ(flow, inlet+1);

    // reacting surface
    int reac_surf = domain_new("reacting-surface", interface, "surface");
    ASSERT_EQ(reac_surf, flow+1);

    // set up stack
    vector<int> doms{inlet, flow, reac_surf};
    int flame = sim1D_new(3, doms.data());
    ASSERT_GE(flame, 0);
    int dom = sim1D_domainIndex(flame, "flow");
    ASSERT_GE(dom, 0);
}

TEST(ctonedim, freeflame)
{
    ct_resetStorage();
    auto gas = newThermo("h2o2.yaml", "ohmech");

    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = soln_thermo(sol);
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
    int ret = thermo_equilibrate(ph, "HP", "auto", 1e-9, 50000, 1000, 0);
    ASSERT_GE(ret, 0);
    double rho_out = thermo_density(ph);
    double Tad = thermo_temperature(ph);
    vector<double> yout(nsp);
    thermo_getMassFractions(ph, nsp, yout.data());

    // flow
    int flow = domain_new("free-flow", sol, "flow");

    // grid
    int nz = 21;
    double lz = 0.02;
    vector<double> z(nz);
    double dz = lz;
    dz /= (double)(nz - 1);
    for (int iz = 0; iz < nz; iz++) {
        z[iz] = iz * dz;
    }
    domain_setupGrid(flow, nz, z.data());

    // inlet
    int reac = domain_new("inlet", sol, "inlet");
    bdry_setMoleFractions(reac, X.c_str());
    bdry_setMdot(reac, uin * rho_in);
    bdry_setTemperature(reac, T);

    // outlet
    int prod = domain_new("outlet", sol, "outlet");
    double uout = bdry_mdot(reac) / rho_out;

    // set up stack
    vector<int> doms{reac, flow, prod};
    int flame = sim1D_new(3, doms.data());
    int dom = sim1D_domainIndex(flame, "flow");
    ASSERT_EQ(dom, 1);

    // set up initial guess
    vector<double> locs{0.0, 0.3, 0.7, 1.0};
    vector<double> value{uin, uin, uout, uout};
    int comp = static_cast<int>(domain_componentIndex(flow, "velocity"));
    sim1D_setProfile(flame, dom, comp, 4, locs.data(), 4, value.data());
    value = {T, T, Tad, Tad};
    comp = static_cast<int>(domain_componentIndex(flow, "T"));
    sim1D_setProfile(flame, dom, comp, 4, locs.data(), 4, value.data());
    for (size_t i = 0; i < nsp; i++) {
        value = {yin[i], yin[i], yout[i], yout[i]};
        int buflen = thermo_getSpeciesName(ph, i, 0, 0) + 1; // include \0
        vector<char> buf(buflen);
        thermo_getSpeciesName(ph, i, buflen, buf.data());
        string name(buf.data());
        ASSERT_EQ(name, gas->speciesName(i));
        comp = static_cast<int>(domain_componentIndex(flow, buf.data()));
        sim1D_setProfile(flame, dom, comp, 4, locs.data(), 4, value.data());
    }

    // simulation settings
    double ratio = 15.0;
    double slope = 0.3;
    double curve = 0.5;
    sim1D_setRefineCriteria(flame, dom, ratio, slope, curve, 0.);
    sim1D_setFixedTemperature(flame, 0.85 * T + .15 * Tad);

    // solve and save
    flow1D_solveEnergyEqn(flow, 1);
    bool refine_grid = false;
    int loglevel = 0;
    sim1D_solve(flame, loglevel, refine_grid);
    ret = sim1D_save(flame, "gtest-freeflame.yaml", "clib",
                     "Solution from CLib interface");
    ASSERT_GE(ret, 0);
    if (usesHDF5()) {
        ret = sim1D_save(flame, "gtest-freeflame.h5", "clib",
                         "Solution from CLib interface");
        ASSERT_GE(ret, 0);
    }

    ASSERT_EQ(domain_nPoints(flow), static_cast<size_t>(nz + 1));
    comp = static_cast<int>(domain_componentIndex(dom, "T"));
    double Tprev = sim1D_value(flame, dom, comp, 0);
    for (size_t n = 0; n < domain_nPoints(flow); n++) {
        T = sim1D_value(flame, dom, comp, static_cast<int>(n));
        ASSERT_GE(T, Tprev);
        Tprev = T;
    }
}
