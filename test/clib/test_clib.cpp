#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/ct.h"

using namespace Cantera;
using ::testing::HasSubstr;

string reportError()
{
    vector<char> output_buf;
    int buflen = ct_getCanteraError(0, output_buf.data()) + 1;
    output_buf.resize(buflen);
    ct_getCanteraError(buflen, output_buf.data());
    return string(output_buf.data());
}

TEST(ct, cabinet_exceptions)
{
    soln_newSolution("h2o2.yaml", "ohmech", "default");
    soln_name(999, 0, 0);

    string err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 999 out of range."));

    soln_thermo(998);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 998 out of range."));

    int ret = soln_del(999);
    ASSERT_EQ(ret, -1);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("delete a non-existing object."));

    int ref = soln_newSolution("h2o2.yaml", "ohmech", "default");
    soln_del(ref);
    int thermo = soln_thermo(ref);
    EXPECT_EQ(thermo, -2);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    ct_resetStorage();
    ret = soln_del(0);
    ASSERT_EQ(ret, -1);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("delete a non-existing object."));
}

TEST(ct, new_solution)
{
    string name = "ohmech";
    int ref = soln_newSolution("h2o2.yaml", name.c_str(), "default");

    int buflen = soln_name(ref, 0, 0) + 1; // include \0
    ASSERT_EQ(buflen, int(name.size() + 1));

    int thermo = soln_thermo(ref);
    ASSERT_EQ(thermo_parent(thermo), ref);

    vector<char> buf(buflen);
    soln_name(ref, buflen, buf.data());
    string solName(buf.data());
    ASSERT_EQ(solName, name);
}

TEST(ct, soln_objects)
{
    ct_resetStorage();

    int ref = soln_newSolution("gri30.yaml", "gri30", "none");
    ASSERT_EQ(ref, 0);
    ASSERT_EQ(thermo_size(), 1); // one ThermoPhase object

    int ref2 = soln_newSolution("h2o2.yaml", "ohmech", "default");
    ASSERT_EQ(ref2, 1);
    ASSERT_EQ(thermo_size(), 2); // two ThermoPhase objects

    int thermo = soln_thermo(ref);
    ASSERT_EQ(thermo_parent(thermo), ref);

    int thermo2 = soln_thermo(ref2);
    ASSERT_EQ(thermo2, 1); // references stored object with index '1'
    ASSERT_EQ(thermo_nSpecies(thermo2), 10u);
    ASSERT_EQ(thermo_parent(thermo2), ref2);

    int kin = soln_kinetics(ref);

    int kin2 = soln_kinetics(ref2);
    ASSERT_EQ(kin2, 1);
    ASSERT_EQ(kin_nReactions(kin2), 29u);
    ASSERT_EQ(kin_parent(kin2), ref2);
    ASSERT_EQ(kin_parent(kin), ref);

    int trans = soln_transport(ref);
    ASSERT_EQ(trans_parent(trans), ref);

    int trans2 = soln_transport(ref2);
    ASSERT_EQ(trans2, 1);
    int buflen = trans_transportModel(trans2, 0, 0);
    vector<char> buf(buflen);
    trans_transportModel(trans2, buflen, buf.data());
    string trName(buf.data());
    ASSERT_EQ(trName, "mixture-averaged");
    ASSERT_EQ(trans_parent(trans2), ref2);

    soln_del(ref2);
    size_t nsp = thermo_nSpecies(thermo2);
    ASSERT_EQ(nsp, npos);
    string err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    size_t nr = thermo_nSpecies(thermo2);
    ASSERT_EQ(nr, npos);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    trans2 = soln_setTransportModel(ref, "mixture-averaged");
    ASSERT_EQ(trans2, 2);
    buflen = trans_transportModel(trans2, 0, 0);
    buf.resize(buflen);
    trans_transportModel(trans2, buflen, buf.data());
    trName = buf.data();
    ASSERT_EQ(trName, "mixture-averaged");
}

TEST(ct, new_interface)
{
    ct_resetStorage();

    int sol = soln_newSolution("ptcombust.yaml", "gas", "none");
    ASSERT_EQ(sol, 0);

    vector<int> adj{sol};
    int surf = soln_newInterface("ptcombust.yaml", "Pt_surf", 1, adj.data());
    ASSERT_EQ(surf, 1);

    int ph_surf = soln_thermo(surf);
    int buflen = soln_name(ph_surf, 0, 0) + 1; // include \0
    vector<char> buf(buflen);
    soln_name(ph_surf, buflen, buf.data());
    string solName(buf.data());
    ASSERT_EQ(solName, "Pt_surf");

    int kin_surf = soln_kinetics(surf);
    buflen = kin_getType(kin_surf, 0, 0) + 1; // include \0
    buf.resize(buflen);
    kin_getType(ph_surf, buflen, buf.data());
    string kinType(buf.data());
    ASSERT_EQ(kinType, "surface");
}

TEST(ct, new_interface_auto)
{
    ct_resetStorage();

    vector<int> adj;
    int surf = soln_newInterface("ptcombust.yaml", "Pt_surf", 0, adj.data());
    ASSERT_EQ(surf, 0);

    ASSERT_EQ(soln_nAdjacent(surf), 1u);
    int gas = soln_adjacent(surf, 0);
    ASSERT_EQ(gas, 1);

    int buflen = soln_name(gas, 0, 0) + 1; // include \0
    vector<char> buf(buflen);
    soln_name(gas, buflen, buf.data());
    string solName(buf.data());
    ASSERT_EQ(solName, "gas");

    buflen = soln_adjacentName(surf, 0, 0, 0) + 1;
    buf.resize(buflen);
    soln_adjacentName(surf, 0, buflen, buf.data());
    solName = buf.data();
    ASSERT_EQ(solName, "gas");
}

TEST(ct, thermo)
{
    int ret;
    int sol = soln_newSolution("gri30.yaml", "gri30", "none");
    int thermo = soln_thermo(sol);
    ASSERT_GE(thermo, 0);
    size_t nsp = thermo_nSpecies(thermo);
    ASSERT_EQ(nsp, 53u);

    ret = thermo_setTemperature(thermo, 500);
    ASSERT_EQ(ret, 0);
    ret = thermo_setPressure(thermo, 5 * 101325);
    ASSERT_EQ(ret, 0);
    ret = thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    ASSERT_EQ(ret, 0);

    ret = thermo_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    ASSERT_EQ(ret, 0);
    double T = thermo_temperature(thermo);
    ASSERT_GT(T, 2200);
    ASSERT_LT(T, 2300);

    size_t ns = thermo_nSpecies(thermo);
    vector<double> work(ns);
    vector<double> X(ns);
    thermo_getMoleFractions(thermo, ns, X.data());
    double prod;

    thermo_getPartialMolarEnthalpies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_enthalpy_mole(thermo), 1e-6);

    thermo_getPartialMolarEntropies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_entropy_mole(thermo), 1e-6);

    thermo_getPartialMolarIntEnergies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_intEnergy_mole(thermo), 1e-6);

    thermo_getPartialMolarCp(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_cp_mole(thermo), 1e-6);

    thermo_getPartialMolarVolumes(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, 1./thermo_molarDensity(thermo), 1e-6);
}

TEST(ct, kinetics)
{
    int sol0 = soln_newSolution("gri30.yaml", "gri30", "none");
    int thermo = soln_thermo(sol0);
    int kin = soln_kinetics(sol0);
    ASSERT_GE(kin, 0);

    size_t nr = kin_nReactions(kin);
    ASSERT_EQ(nr, 325u);

    thermo_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    double T = thermo_temperature(thermo);
    thermo_setTemperature(thermo, T - 200);

    auto sol = newSolution("gri30.yaml", "gri30", "none");
    auto phase = sol->thermo();
    auto kinetics = sol->kinetics();

    phase->equilibrate("HP");
    ASSERT_NEAR(T, phase->temperature(), 1e-2);
    phase->setTemperature(T - 200);

    vector<double> c_ropf(nr);
    kin_getFwdRatesOfProgress(kin, 325, c_ropf.data());
    vector<double> cpp_ropf(nr);
    kinetics->getFwdRatesOfProgress(cpp_ropf.data());

    for (size_t n = 0; n < nr; n++) {
        ASSERT_NEAR(cpp_ropf[n], c_ropf[n], 1e-6);
    }
}

TEST(ct, transport)
{
    int sol0 = soln_newSolution("gri30.yaml", "gri30", "default");
    int thermo = soln_thermo(sol0);
    int tran = soln_transport(sol0);

    size_t nsp = thermo_nSpecies(thermo);
    vector<double> c_dkm(nsp);
    int ret = trans_getMixDiffCoeffs(tran, 53, c_dkm.data());
    ASSERT_EQ(ret, 0);

    vector<double> cpp_dkm(nsp);
    auto sol = newSolution("gri30.yaml", "gri30");
    auto transport = sol->transport();
    transport->getMixDiffCoeffs(cpp_dkm.data());

    for (size_t n = 0; n < nsp; n++) {
        ASSERT_NEAR(cpp_dkm[n], c_dkm[n], 1e-10);
    }
}


int main(int argc, char** argv)
{
    printf("Running main() from test_clib.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    make_deprecation_warnings_fatal();
    printStackTraceOnSegfault();
    Cantera::CanteraError::setStackTraceDepth(20);
    vector<string> fileNames = {"gtest-freeflame.yaml", "gtest-freeflame.h5"};
    for (const auto& fileName : fileNames) {
        if (std::ifstream(fileName).good()) {
            std::remove(fileName.c_str());
        }
    }
    int result = RUN_ALL_TESTS();
    appdelete();
    return result;
}
