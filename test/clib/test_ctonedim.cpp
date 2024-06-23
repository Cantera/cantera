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

    int flow = domain_new("free-flow", sol, "flow");
    ASSERT_GE(flow, 0);
    domain_setID(flow, "flow");
    ASSERT_NEAR(flow1D_pressure(flow), P, 1e-5);

    int buflen = domain_type3(flow, 0, 0);
    char* buf = new char[buflen];
    domain_type3(flow, buflen, buf);
    string domName = buf;
    ASSERT_EQ(domName, "free-flow");
    delete[] buf;
}

TEST(ctonedim, freeflow_from_parts)
{
    ct_resetStorage();

    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = soln_thermo(sol);
    ASSERT_GE(ph, 0);
    int kin = soln_kinetics(sol);
    ASSERT_GE(kin, 0);
    int tr = soln_transport(sol);
    ASSERT_GE(tr, 0);

    double T = 1050;
    double P = 5 * 101325;
    string X = "CH4:1.0, O2:2.0, N2:7.52";
    thermo_setMoleFractionsByName(ph, X.c_str());
    thermo_setTemperature(ph, T);
    thermo_setPressure(ph, P);

    int itype = 2; // free flow
    int flow = flow1D_new(ph, kin, tr, itype);
    ASSERT_GE(flow, 0);
    domain_setID(flow, "flow");
    ASSERT_NEAR(flow1D_pressure(flow), P, 1e-5);

    int buflen = domain_type3(flow, 0, 0);
    char* buf = new char[buflen];
    domain_type3(flow, buflen, buf);
    string domName = buf;
    ASSERT_EQ(domName, "free-flow");
    delete[] buf;
}

TEST(ctonedim, inlet)
{
    int inlet = inlet_new();
    ASSERT_GE(inlet, 0);

    int buflen = domain_type3(inlet, 0, 0);
    char* buf = new char[buflen];
    domain_type3(inlet, buflen, buf);
    string domName = buf;
    ASSERT_EQ(domName, "inlet");
    delete[] buf;
}

TEST(ctonedim, outlet)
{
    int index = outlet_new();
    ASSERT_GE(index, 0);
}

TEST(ctonedim, reacting_surface)
{
    int surf = soln_newInterface("ptcombust.yaml", "Pt_surf", 0, 0);
    int index = domain_new("reacting-surface", surf, "surface");
    ASSERT_GE(index, 0);
}

TEST(ctonedim, reacting_surface_from_parts)
{
    int index = reactingsurf_new();
    ASSERT_GE(index, 0);

    int gas = thermo_newFromFile("ptcombust.yaml", "gas");
    int surf = thermo_newFromFile("ptcombust.yaml", "Pt_surf");
    int kin = kin_newFromFile("ptcombust.yaml", "", surf, gas, -1, -1, -1);
    ASSERT_GE(kin, 0);

    int ret = reactingsurf_setkineticsmgr(index, kin);
    ASSERT_EQ(ret, 0);
}

TEST(ctonedim, catcomb_stack)
{
    int sol = soln_newSolution("ptcombust.yaml", "gas", "default");
    int gas = soln_thermo(sol);
    int gas_kin = soln_kinetics(sol);
    int trans = soln_transport(sol);

    int surf = thermo_newFromFile("ptcombust.yaml", "Pt_surf");
    int surf_kin = kin_newFromFile("ptcombust.yaml", "", surf, gas, -1, -1, -1);

    // inlet
    int inlet = inlet_new();

    // flow
    int itype = 1; // free flow
    int flow = flow1D_new(gas, gas_kin, trans, itype);
    domain_setID(flow, "flow");

    // reacting surface
    int reac_surf = reactingsurf_new();
    int ret = reactingsurf_setkineticsmgr(reac_surf, surf_kin);
    ASSERT_EQ(ret, 0);

    // set up stack
    vector<int> doms{inlet, flow, reac_surf};
    int flame = sim1D_new(3, doms.data());
    ASSERT_GE(flame, 0);
    int dom = sim1D_domainIndex(flame, "flow");
    ASSERT_GE(dom, 0);
}

TEST(ctonedim, freeflame_from_parts)
{
    ct_resetStorage();
    auto gas = newThermo("h2o2.yaml", "ohmech");

    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = soln_thermo(sol);
    int kin = soln_kinetics(sol);
    int tr = soln_transport(sol);
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
    int ret = thermo_equilibrate(ph, "HP", 0, 1e-9, 50000, 1000, 0);
    ASSERT_GE(ret, 0);
    double rho_out = thermo_density(ph);
    double Tad = thermo_temperature(ph);
    vector<double> yout(nsp);
    thermo_getMassFractions(ph, nsp, yout.data());

    // flow
    int itype = 2; // free flow
    int flow = flow1D_new(ph, kin, tr, itype);
    domain_setID(flow, "flow");

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
    int reac = inlet_new();
    domain_setID(reac, "inlet");
    bdry_setMoleFractions(reac, X.c_str());
    bdry_setMdot(reac, uin * rho_in);
    bdry_setTemperature(reac, T);

    // outlet
    int prod = outlet_new();
    domain_setID(prod, "outlet");
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
        char* buf = new char[buflen];
        thermo_getSpeciesName(ph, i, buflen, buf);
        string name = buf;
        ASSERT_EQ(name, gas->speciesName(i));
        comp = static_cast<int>(domain_componentIndex(flow, buf));
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

TEST(ctonedim, stflow_tests)
{
    //! @todo: To be removed after Cantera 3.1
    ct_resetStorage();
    auto gas = newThermo("h2o2.yaml", "ohmech");

    int sol = soln_newSolution("h2o2.yaml", "ohmech", "default");
    int ph = soln_thermo(sol);
    int kin = soln_kinetics(sol);
    int tr = soln_transport(sol);

    // spot check some errors
    int itype = 2; // free flow
    int ret = stflow_new(ph, kin, tr, itype);
    ASSERT_EQ(ret, -1);  // -1 is an error code

    int flow = flow1D_new(ph, kin, tr, itype);
    ASSERT_EQ(stflow_setTransport(flow, tr), -1);
    ASSERT_EQ(stflow_pressure(flow), DERR);  // DERR is an error code
    ASSERT_EQ(stflow_setPressure(flow, OneAtm), -1);
}
