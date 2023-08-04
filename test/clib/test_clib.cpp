#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/ct.h"

using namespace Cantera;
using ::testing::HasSubstr;

string reportError()
{
    int buflen = 0;
    char* output_buf = 0;
    buflen = ct_getCanteraError(buflen, output_buf) + 1;
    output_buf = new char[buflen];
    ct_getCanteraError(buflen, output_buf);
    string err = output_buf;
    delete[] output_buf;
    return err;
}

TEST(ct, cabinet_exceptions)
{
    soln_newSolution("h2o2.yaml", "ohmech", "");
    soln_name(999, 0, 0);

    string err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 999 out of range."));

    int ret = soln_del(999);
    ASSERT_EQ(ret, -1);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("delete a non-existing object."));

    ct_resetStorage();
    ret = soln_del(0);
    ASSERT_EQ(ret, -1);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("delete a non-existing object."));
}

TEST(ct, new_solution)
{
    string name = "ohmech";
    int ref = soln_newSolution("h2o2.yaml", name.c_str(), "");

    int buflen = soln_name(ref, 0, 0) + 1; // include \0
    ASSERT_EQ(buflen, int(name.size() + 1));

    char* buf = new char[buflen];
    soln_name(ref, buflen, buf);
    string solName = buf;
    ASSERT_EQ(solName, name);
    delete[] buf;
}

TEST(ct, soln_objects)
{
    ct_resetStorage();

    thermo_newFromFile("gri30.yaml", "gri30");

    int ref = soln_newSolution("gri30.yaml", "gri30", "none");
    ASSERT_EQ(ref, 0);
    int ref2 = soln_newSolution("h2o2.yaml", "ohmech", "");
    ASSERT_EQ(ref2, 1);

    int thermo = soln_thermo(ref2);
    ASSERT_EQ(thermo, 2);
    ASSERT_EQ(thermo_nSpecies(thermo), 10u);

    int kin = soln_kinetics(ref2);
    ASSERT_EQ(kin, 1);
    ASSERT_EQ(kin_nReactions(kin), 29u);

    int trans = soln_kinetics(ref2);
    ASSERT_EQ(trans, 1);
    int buflen = trans_transportModel(trans, 0, 0);
    vector<char> buf(buflen);
    trans_transportModel(trans, buflen, buf.data());
    string trName(buf.data());
    ASSERT_EQ(trName, "mixture-averaged");

    soln_del(ref2);
    size_t nsp = thermo_nSpecies(thermo);
    ASSERT_EQ(nsp, npos);
    string err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    size_t nr = thermo_nSpecies(thermo);
    ASSERT_EQ(nr, npos);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    trans = soln_setTransportModel(ref, "mixture-averaged");
    ASSERT_EQ(trans, 2);
    buflen = trans_transportModel(trans, 0, 0);
    buf.resize(buflen);
    trans_transportModel(trans, buflen, buf.data());
    trName = buf.data();
    ASSERT_EQ(trName, "mixture-averaged");
}

TEST(ct, new_interface)
{
    ct_resetStorage();

    int sol = soln_newSolution("ptcombust.yaml", "gas", "");
    ASSERT_EQ(sol, 0);

    vector<int> adj{sol};
    int surf = soln_newInterface("ptcombust.yaml", "Pt_surf", 1, adj.data());
    ASSERT_EQ(surf, 1);

    int ph_surf = soln_thermo(surf);
    int buflen = soln_name(ph_surf, 0, 0) + 1; // include \0
    char* buf = new char[buflen];
    soln_name(ph_surf, buflen, buf);
    string solName = buf;
    ASSERT_EQ(solName, "Pt_surf");
    delete[] buf;

    int kin_surf = soln_kinetics(surf);
    buflen = kin_getType(kin_surf, 0, 0) + 1; // include \0
    buf = new char[buflen];
    kin_getType(ph_surf, buflen, buf);
    string kinType = buf;
    ASSERT_EQ(kinType, "surface");
    delete[] buf;
}

TEST(ct, new_interface_auto)
{
    ct_resetStorage();

    vector<int> adj;
    int surf = soln_newInterface("ptcombust.yaml", "Pt_surf", 0, adj.data());
    ASSERT_EQ(surf, 1);

    ASSERT_EQ(soln_nAdjacent(surf), 1u);
    int gas = soln_adjacent(surf, 0);
    ASSERT_EQ(gas, 0);

    int buflen = soln_name(gas, 0, 0) + 1; // include \0
    char* buf = new char[buflen];
    soln_name(gas, buflen, buf);
    string solName = buf;
    ASSERT_EQ(solName, "gas");
    delete[] buf;
}

TEST(ct, thermo)
{
    int ret;
    int thermo = thermo_newFromFile("gri30.yaml", "gri30");
    ASSERT_GE(thermo, 0);
    size_t nsp = thermo_nSpecies(thermo);
    ASSERT_EQ(nsp, 53u);

    ret = thermo_setTemperature(thermo, 500);
    ASSERT_EQ(ret, 0);
    ret = thermo_setPressure(thermo, 5 * 101325);
    ASSERT_EQ(ret, 0);
    ret = thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    ASSERT_EQ(ret, 0);

    ret = thermo_equilibrate(thermo, "HP", 0, 1e-9, 50000, 1000, 0);
    ASSERT_EQ(ret, 0);
    double T = thermo_temperature(thermo);
    ASSERT_GT(T, 2200);
    ASSERT_LT(T, 2300);
}

TEST(ct, kinetics)
{
    int thermo = thermo_newFromFile("gri30.yaml", "gri30");
    int kin = kin_newFromFile("gri30.yaml", "", thermo, -1, -1, -1, -1);
    ASSERT_GE(kin, 0);

    size_t nr = kin_nReactions(kin);
    ASSERT_EQ(nr, 325u);

    thermo_equilibrate(thermo, "HP", 0, 1e-9, 50000, 1000, 0);
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
    int thermo = thermo_newFromFile("gri30.yaml", "gri30");
    int tran = trans_newDefault(thermo, 0);
    ASSERT_GE(tran, 0);

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
