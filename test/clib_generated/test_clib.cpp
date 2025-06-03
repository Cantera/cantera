#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/ctkin.h"
#include "cantera_clib/cttrans.h"

using namespace Cantera;
using ::testing::HasSubstr;

string reportError()
{
    vector<char> output_buf;
    int32_t buflen = ct_getCanteraError(0, output_buf.data());
    output_buf.resize(buflen);
    ct_getCanteraError(buflen, output_buf.data());
    return string(output_buf.data());
}

TEST(ct, cabinet_exceptions)
{
    sol_newSolution("h2o2.yaml", "ohmech", "default");
    sol_name(999, 0, 0);

    string err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 999 out of range."));

    sol_thermo(998);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 998 out of range."));

    int32_t ret = sol_del(997);
    ASSERT_EQ(ret, -1);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 997 out of range."));

    int32_t ref = sol_newSolution("h2o2.yaml", "ohmech", "default");
    sol_del(ref);
    int32_t thermo = sol_thermo(ref);
    EXPECT_EQ(thermo, -2);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    ct_resetStorage();
    ret = sol_del(0);
    ASSERT_EQ(ret, -1);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("Index 0 out of range."));
}

TEST(ct, new_solution)
{
    ct_resetStorage();

    string name = "ohmech";
    int32_t ref = sol_newSolution("h2o2.yaml", name.c_str(), "default");
    ASSERT_EQ(ref, 0);

    ASSERT_EQ(sol_cabinetSize(), 1);
    ASSERT_EQ(thermo_cabinetSize(), 1);
    ASSERT_EQ(kin_cabinetSize(), 1);

    int32_t buflen = sol_name(ref, 0, 0); // includes \0
    ASSERT_EQ(buflen, static_cast<int32_t>(name.size() + 1));

    int32_t thermo = sol_thermo(ref);
    ASSERT_EQ(thermo_parentHandle(thermo), ref);

    vector<char> buf(buflen);
    sol_name(ref, buflen, buf.data());
    string solName(buf.data());
    ASSERT_EQ(solName, name);
}

TEST(ct, sol_objects)
{
    ct_resetStorage();

    int32_t ref = sol_newSolution("gri30.yaml", "gri30", "none");
    ASSERT_EQ(ref, 0);
    ASSERT_EQ(thermo_cabinetSize(), 1); // one ThermoPhase object

    int32_t ref2 = sol_newSolution("h2o2.yaml", "ohmech", "default");
    ASSERT_EQ(ref2, 1);
    ASSERT_EQ(thermo_cabinetSize(), 2); // two ThermoPhase objects

    int32_t thermo = sol_thermo(ref);
    ASSERT_EQ(thermo_parentHandle(thermo), ref);

    int32_t thermo2 = sol_thermo(ref2);
    ASSERT_EQ(thermo2, 1); // references stored object with index '1'
    ASSERT_EQ(thermo_nSpecies(thermo2), 10u);
    ASSERT_EQ(thermo_parentHandle(thermo2), ref2);

    int32_t kin = sol_kinetics(ref);

    int32_t kin2 = sol_kinetics(ref2);
    ASSERT_EQ(kin2, 1);
    ASSERT_EQ(kin_nReactions(kin2), 29u);
    ASSERT_EQ(kin_parentHandle(kin2), ref2);
    ASSERT_EQ(kin_parentHandle(kin), ref);

    int32_t trans = sol_transport(ref);
    ASSERT_EQ(trans_parentHandle(trans), ref);

    int32_t trans2 = sol_transport(ref2);
    ASSERT_EQ(trans2, 1);
    int32_t buflen = trans_transportModel(trans2, 0, 0);
    vector<char> buf(buflen);
    trans_transportModel(trans2, buflen, buf.data());
    string trName(buf.data());
    ASSERT_EQ(trName, "mixture-averaged");
    ASSERT_EQ(trans_parentHandle(trans2), ref2);

    sol_del(ref2);
    int32_t nsp = thermo_nSpecies(thermo2);
    ASSERT_EQ(nsp, ERR);
    string err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    nsp = thermo_nSpecies(thermo2);
    ASSERT_EQ(nsp, ERR);
    err = reportError();
    EXPECT_THAT(err, HasSubstr("has been deleted."));

    trans2 = sol_setTransportModel(ref, "mixture-averaged");
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

    int32_t sol = sol_newSolution("ptcombust.yaml", "gas", "none");
    ASSERT_EQ(sol, 0);

    vector<int32_t> adj{sol};
    int32_t surf = sol_newInterface("ptcombust.yaml", "Pt_surf", 1, adj.data());
    ASSERT_EQ(surf, 1);

    int32_t ph_surf = sol_thermo(surf);
    int32_t buflen = sol_name(ph_surf, 0, 0) + 1; // include \0
    vector<char> buf(buflen);
    sol_name(ph_surf, buflen, buf.data());
    string solName(buf.data());
    ASSERT_EQ(solName, "Pt_surf");

    int32_t kin_surf = sol_kinetics(surf);
    buflen = kin_getType(kin_surf, 0, 0) + 1; // include \0
    buf.resize(buflen);
    kin_getType(ph_surf, buflen, buf.data());
    string kinType(buf.data());
    ASSERT_EQ(kinType, "surface");
}

TEST(ct, new_interface_auto)
{
    ct_resetStorage();

    vector<int32_t> adj;
    int32_t surf = sol_newInterface("ptcombust.yaml", "Pt_surf", 0, adj.data());
    ASSERT_EQ(surf, 0);

    ASSERT_EQ(sol_nAdjacent(surf), 1u);
    int32_t gas = sol_adjacent(surf, 0);
    ASSERT_EQ(gas, 1);

    int32_t buflen = sol_name(gas, 0, 0) + 1; // include \0
    vector<char> buf(buflen);
    sol_name(gas, buflen, buf.data());
    string solName(buf.data());
    ASSERT_EQ(solName, "gas");

    buflen = sol_adjacentName(surf, 0, 0, 0) + 1;
    buf.resize(buflen);
    sol_adjacentName(surf, 0, buflen, buf.data());
    solName = buf.data();
    ASSERT_EQ(solName, "gas");
}

TEST(ct, transport)
{
    int32_t sol0 = sol_newSolution("gri30.yaml", "gri30", "default");
    int32_t thermo = sol_thermo(sol0);
    int32_t tran = sol_transport(sol0);

    size_t nsp = thermo_nSpecies(thermo);
    vector<double> c_dkm(nsp);
    int32_t ret = trans_getMixDiffCoeffs(tran, 53, c_dkm.data());
    ASSERT_EQ(ret, 0);

    vector<double> cpp_dkm(nsp);
    auto sol = newSolution("gri30.yaml", "gri30");
    auto transport = sol->transport();
    transport->getMixDiffCoeffs(cpp_dkm.data());

    for (size_t n = 0; n < nsp; n++) {
        ASSERT_NEAR(cpp_dkm[n], c_dkm[n], 1e-10);
    }
}

TEST(ct, constants)
{
    ASSERT_EQ(ct_Avogadro(), Avogadro);
    ASSERT_EQ(ct_GasConstant(), GasConstant);
    ASSERT_EQ(ct_OneAtm(), OneAtm);
}


int main(int argc, char** argv)
{
    printf("Running main() from test_clib3.cpp\n");
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
