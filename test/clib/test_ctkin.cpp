#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/kinetics/Reaction.h"
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


TEST(ctkin, kinetics)
{
    int32_t sol0 = sol_newSolution("gri30.yaml", "gri30", "none");
    int32_t thermo = sol_thermo(sol0);
    int32_t kin = sol_kinetics(sol0);
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

TEST(ctkin, exceptions)
{
    int32_t sol0 = sol_newSolution("h2o2.yaml", "", "none");
    int32_t kin = sol_kinetics(sol0);

    kin_reaction(kin, 998);
    string err = reportError();
    EXPECT_THAT(err, HasSubstr("IndexError: 998 outside valid range of 0 to"));
}

TEST(ctrxn, reaction)
{
    int32_t sol0 = sol_newSolution("h2o2.yaml", "", "none");
    int32_t kin = sol_kinetics(sol0);
    int32_t rxn = kin_reaction(kin, 0);

    auto sol = newSolution("h2o2.yaml", "", "none");
    auto kinetics = sol->kinetics();
    auto reaction = kinetics->reaction(0);

    int32_t buflen = rxn_type(rxn, 0, 0);
    vector<char> buf(buflen);
    rxn_type(rxn, buflen, buf.data());
    string text = buf.data();
    ASSERT_EQ(text, reaction->type());
    ASSERT_EQ(text, "three-body-Arrhenius");

    buflen = rxn_equation(rxn, 0, 0);
    buf.resize(buflen);
    rxn_equation(rxn, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, reaction->equation());
    ASSERT_EQ(text, "2 O + M <=> O2 + M");
    ASSERT_EQ(rxn_usesThirdBody(rxn), 1);
    ASSERT_EQ(rxn_allowNonreactantOrders(rxn), 0);
    ASSERT_EQ(rxn_allowNonreactantOrders(rxn),
              static_cast<int32_t>(reaction->allow_nonreactant_orders));

    buflen = rxn_id(rxn, 0, 0);
    buf.resize(buflen);
    rxn_id(rxn, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, reaction->id);
    ASSERT_EQ(text, "");

    rxn_setId(rxn, "spam");
    buflen = rxn_id(rxn, 0, 0);
    buf.resize(buflen);
    rxn_id(rxn, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, "spam");
}

TEST(ctrdiag, diagram)
{
    int32_t sol = sol_newSolution("h2o2.yaml", "", "none");
    int32_t thermo = sol_thermo(sol);
    int32_t kin = sol_kinetics(sol);
    thermo_setState_TP(thermo, 1000., ct_OneAtm());

    int32_t diag = rdiag_newReactionPathDiagram(kin, "H");
    ASSERT_EQ(rdiag_boldThreshold(diag), 0.2);

    rdiag_setFlowType(diag, "spam");
    string text = reportError();
    EXPECT_THAT(text, HasSubstr("Unknown flow type 'spam'"));
    int32_t buflen = rdiag_flowType(diag, 0, 0);
    vector<char> buf(buflen);
    rdiag_flowType(diag, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, "NetFlow");
    rdiag_setFlowType(diag, "OneWayFlow");
    buflen = rdiag_flowType(diag, 0, 0);
    buf.resize(buflen);
    rdiag_flowType(diag, buflen, buf.data());
    text = buf.data();
    ASSERT_EQ(text, "OneWayFlow");

    buflen = rdiag_getLog(diag, 0, 0);
    buf.resize(buflen);
    rdiag_getLog(diag, buflen, buf.data());
    text = buf.data();
    EXPECT_THAT(text, StartsWith("\nReaction 1:"));

    buflen = rdiag_getData(diag, 0, 0);
    buf.resize(buflen);
    rdiag_getData(diag, buflen, buf.data());
    text = buf.data();
    EXPECT_THAT(text, StartsWith("\nH H2 \nH H2"));

    buflen = rdiag_getDot(diag, 0, 0);
    buf.resize(buflen);
    rdiag_getDot(diag, buflen, buf.data());
    text = buf.data();
    EXPECT_THAT(text, StartsWith("digraph reaction_paths"));
}
