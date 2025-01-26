#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/clib/clib_defs.h"
#include "cantera/clib_experimental/ct3.h"
#include "cantera/clib_experimental/ctsol3.h"
#include "cantera/clib_experimental/ctthermo3.h"
#include "cantera/clib_experimental/ctkin3.h"
#include "cantera/clib_experimental/ctrxn3.h"

using namespace Cantera;
using ::testing::HasSubstr;

string reportError();  // forward declaration


TEST(ctkin3, kinetics)
{
    int sol0 = sol3_newSolution("gri30.yaml", "gri30", "none");
    int thermo = sol3_thermo(sol0);
    int kin = sol3_kinetics(sol0);
    ASSERT_GE(kin, 0);

    size_t nr = kin3_nReactions(kin);
    ASSERT_EQ(nr, 325u);

    thermo3_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    double T = thermo3_temperature(thermo);
    thermo3_setTemperature(thermo, T - 200);

    auto sol = newSolution("gri30.yaml", "gri30", "none");
    auto phase = sol->thermo();
    auto kinetics = sol->kinetics();

    phase->equilibrate("HP");
    ASSERT_NEAR(T, phase->temperature(), 1e-2);
    phase->setTemperature(T - 200);

    vector<double> c_ropf(nr);
    kin3_getFwdRatesOfProgress(kin, 325, c_ropf.data());
    vector<double> cpp_ropf(nr);
    kinetics->getFwdRatesOfProgress(cpp_ropf.data());

    for (size_t n = 0; n < nr; n++) {
        ASSERT_NEAR(cpp_ropf[n], c_ropf[n], 1e-6);
    }
}

TEST(ctkin3, exceptions)
{
    int sol0 = sol3_newSolution("h2o2.yaml", "", "none");
    int kin = sol3_kinetics(sol0);

    kin3_reaction(kin, 998);
    string err = reportError();
    EXPECT_THAT(err, HasSubstr("IndexError: 998 outside valid range of 0 to 28."));
}

TEST(ctkin3, reaction)
{
    int sol0 = sol3_newSolution("h2o2.yaml", "", "none");
    int kin = sol3_kinetics(sol0);
    int rxn = kin3_reaction(kin, 0);

    auto sol = newSolution("h2o2.yaml", "", "none");
    auto kinetics = sol->kinetics();
    auto reaction = kinetics->reaction(0);

    int buflen = rxn3_type(rxn, 0, 0);
    vector<char> buf(buflen);
    rxn3_type(rxn, buflen, buf.data());
    string rxnType = buf.data();
    ASSERT_EQ(rxnType, reaction->type());
    ASSERT_EQ(rxnType, "three-body-Arrhenius");

    buflen = rxn3_equation(rxn, 0, 0);
    buf.resize(buflen);
    rxn3_equation(rxn, buflen, buf.data());
    string rxnEqn = buf.data();
    ASSERT_EQ(rxnEqn, reaction->equation());
    ASSERT_EQ(rxnEqn, "2 O + M <=> O2 + M");
    ASSERT_EQ(rxn3_usesThirdBody(rxn), 1);
}
