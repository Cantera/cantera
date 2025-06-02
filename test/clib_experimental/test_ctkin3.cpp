#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/clib/clib_defs.h"
#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"
#include "cantera_clib/ctkin.h"
#include "cantera_clib/ctrxn.h"
#include "cantera_clib/ctrdiag.h"

using namespace Cantera;
using ::testing::HasSubstr;
using ::testing::StartsWith;

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
    EXPECT_THAT(err, HasSubstr("IndexError: 998 outside valid range of 0 to"));
}

TEST(ctrxn3, reaction)
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
    string text = buf.data();
    ASSERT_EQ(text, reaction->type());
    ASSERT_EQ(text, "three-body-Arrhenius");

    buflen = rxn3_equation(rxn, 0, 0);
    buf.resize(buflen);
    rxn3_equation(rxn, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, reaction->equation());
    ASSERT_EQ(text, "2 O + M <=> O2 + M");
    ASSERT_EQ(rxn3_usesThirdBody(rxn), 1);
    ASSERT_EQ(rxn3_allowNonreactantOrders(rxn), 0);
    ASSERT_EQ(rxn3_allowNonreactantOrders(rxn),
              int(reaction->allow_nonreactant_orders));

    buflen = rxn3_id(rxn, 0, 0);
    buf.resize(buflen);
    rxn3_id(rxn, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, reaction->id);
    ASSERT_EQ(text, "");

    rxn3_setId(rxn, "spam");
    buflen = rxn3_id(rxn, 0, 0);
    buf.resize(buflen);
    rxn3_id(rxn, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, "spam");
}

TEST(ctrdiag3, diagram)
{
    int sol = sol3_newSolution("h2o2.yaml", "", "none");
    int thermo = sol3_thermo(sol);
    int kin = sol3_kinetics(sol);
    thermo3_set_TP(thermo, 1000., ct3_OneAtm());

    int diag = rdiag3_newReactionPathDiagram(kin, "H");
    ASSERT_EQ(rdiag3_boldThreshold(diag), 0.2);

    rdiag3_setFlowType(diag, "spam");
    string text = reportError();
    EXPECT_THAT(text, HasSubstr("Unknown flow type 'spam'"));
    int buflen = rdiag3_flowType(diag, 0, 0);
    vector<char> buf(buflen);
    rdiag3_flowType(diag, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, "NetFlow");
    rdiag3_setFlowType(diag, "OneWayFlow");
    buflen = rdiag3_flowType(diag, 0, 0);
    buf.resize(buflen);
    rdiag3_flowType(diag, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, "OneWayFlow");

    buflen = rdiag3_getLog(diag, 0, 0);
    buf.resize(buflen);
    rdiag3_getLog(diag, buflen, buf.data());
    text = buf.data();
    EXPECT_THAT(text, StartsWith("\nReaction 1:"));

    buflen = rdiag3_getData(diag, 0, 0);
    buf.resize(buflen);
    rdiag3_getData(diag, buflen, buf.data());
    text = buf.data();
    EXPECT_THAT(text, StartsWith("\nH H2 \nH H2"));

    buflen = rdiag3_getDot(diag, 0, 0);
    buf.resize(buflen);
    rdiag3_getDot(diag, buflen, buf.data());
    text = buf.data();
    EXPECT_THAT(text, StartsWith("digraph reaction_paths"));
}
