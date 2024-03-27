#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/base/Solution.h"
#include "cantera/base/Interface.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/ReactionRateFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/Custom.h"
#include "cantera/kinetics/ElectronCollisionPlasmaRate.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/InterfaceRate.h"
#include "cantera/kinetics/PlogRate.h"
#include "cantera/kinetics/TwoTempPlasmaRate.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/Array.h"

using namespace Cantera;

TEST(ReactionRate, ModifyArrheniusRate)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    const auto& rr = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_TRUE(rr->allowNegativePreExponentialFactor());
    rr->setAllowNegativePreExponentialFactor(false);
    EXPECT_FALSE(rr->allowNegativePreExponentialFactor());
}

TEST(ReactionRate, ArrheniusUnits)
{
    AnyMap rxn1 = AnyMap::fromYamlString(
        "{rate-constant: [2.70000E+13 cm^3/mol/s, 0, 355 cal/mol]}");
    ASSERT_THROW(newReactionRate(rxn1), InputFileError);

    AnyMap rxn2 = AnyMap::fromYamlString(
        "{units: {quantity: mol, length: cm},\n"
         "nested: {rate-constant: [27.0, 0, 355 cal/mol]}}");
    rxn2.applyUnits();
    ASSERT_THROW(newReactionRate(rxn2["nested"].as<AnyMap>()), InputFileError);

    AnyMap rxn3 = AnyMap::fromYamlString(
        "{rate-constant: [2.70000E+13 m^3/kmol/s, 0, 355 cal/mol]}");
    auto rate = newReactionRate(rxn3);
}

TEST(Reaction, ElementaryFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_FALSE(R->usesThirdBody());
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->type(), "Arrhenius");
    EXPECT_EQ(R->rate()->type(), "Arrhenius");
    EXPECT_FALSE(R->allow_negative_orders);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_TRUE(rate->allowNegativePreExponentialFactor());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), -2.7e10);
    EXPECT_DOUBLE_EQ(rate->activationEnergy(), 355 * 4184.0);
}

TEST(Reaction, ThreeBodyFromYaml1)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O + M = O2 + M,"
        " type: three-body,"
        " rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
        " efficiencies: {AR: 0.83, H2O: 5}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);
    EXPECT_EQ(R->type(), "three-body-Arrhenius");
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiencies["H2O"], 5.0);
    EXPECT_DOUBLE_EQ(R->thirdBody()->default_efficiency, 1.0);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 1.2e11);
}

TEST(Reaction, ThreeBodyFromYaml2)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O + M = O2 + M,"
        " rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
        " efficiencies: {AR: 0.83, H2O: 5}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);
    EXPECT_EQ(R->type(), "three-body-Arrhenius");
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiencies["H2O"], 5.0);
    EXPECT_DOUBLE_EQ(R->thirdBody()->default_efficiency, 1.0);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 1.2e11);

    AnyMap input = R->parameters(false);
    EXPECT_FALSE(input.hasKey("type"));
    EXPECT_TRUE(input.hasKey("efficiencies"));

    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 2u);
    EXPECT_EQ(efficiencies["AR"], 0.83);
    EXPECT_EQ(efficiencies["H2O"], 5.);
}

TEST(Reaction, ThreeBodyFromYaml3)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + M <=> CH2(S) + M,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0}}");

    auto R = new Reaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->type(), "three-body-Arrhenius");
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->thirdBody()->name(), "M");

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 5.0e+9);
}

TEST(Reaction, ThreeBodyFromYaml4)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2(S) + O2,"
        " type: three-body,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->type(), "three-body-Arrhenius");
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->thirdBody()->name(), "O2");

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 5.0e+9);

    AnyMap input = R->parameters(false);
    EXPECT_EQ(input.getString("type", ""), "three-body");
    EXPECT_FALSE(input.hasKey("efficiencies"));
    EXPECT_FALSE(input.hasKey("default-efficiency"));
}

TEST(Reaction, ThreeBodyFromYaml5)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0},"
        " efficiencies: {O2: 1.}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->type(), "three-body-Arrhenius");
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->thirdBody()->name(), "O2");
    EXPECT_EQ(R->thirdBody()->default_efficiency, 0.);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 5.0e+9);

    AnyMap input = R->parameters(false);
    EXPECT_FALSE(input.hasKey("type"));
    EXPECT_TRUE(input.hasKey("efficiencies"));
    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 1u);
    EXPECT_EQ(efficiencies.begin()->first, "O2");
    EXPECT_FALSE(input.hasKey("default-efficiency"));
}

TEST(Reaction, ThreeBodyFromYaml6)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " type: three-body,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0},"
        " default-efficiency: 0.,"
        " efficiencies: {O2: 1.}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->type(), "three-body-Arrhenius");
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->thirdBody()->name(), "O2");
    EXPECT_EQ(R->thirdBody()->default_efficiency, 0.);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 5.0e+9);

    AnyMap input = R->parameters(false);
    EXPECT_EQ(input.getString("type", ""), "three-body");
    EXPECT_TRUE(input.hasKey("efficiencies"));
    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 1u);
    EXPECT_EQ(efficiencies.begin()->first, "O2");
    EXPECT_FALSE(input.hasKey("default-efficiency"));
}

TEST(Reaction, ThreeBodyFromYamlMissingM)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O <=> O2," // Missing "M" on each side of the equation
        " type: three-body,"
        " rate-constant: [1.20000E+17, -1, 0]}");

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYamlMultiple)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " type: three-body," // ambiguous valid explicit third bodies
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0}}");

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYamlIncompatible1)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0},"
        " default-efficiency: 0.,"
        " efficiencies: {AR: 1.}}"); // incompatible third body

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYamlIncompatible2)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2(S) + O2,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0},"
        " default-efficiency: 0.,"
        " efficiencies: {AR: 1.}}"); // incompatible single third body

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYamlIncompatible3)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0},"
        " default-efficiency: 1.," // Needs to be zero
        " efficiencies: {O2: 1.}}");

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYamlIncompatible4)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0},"
        " default-efficiency: 0.}"); // missing efficiencies field

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYamlIncompatible5)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2 + O2,"
        " type: three-body,"
        " rate-constant: {A: 5.0e+9, b: 0.0, Ea: 0.0}}"); // third body ambiguous

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, FalloffFromYaml1)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N2O (+M) = N2 + O (+ M),"
        " type: falloff,"
        " high-P-rate-constant: [7.91000E+10, 0, 56020],"
        " low-P-rate-constant: [6.37000E+14, 0, 56640],"
        " SRI: {A: 1.1, B: 700.0, C: 1234.0, D: 56.0, E: 0.7},"
        " efficiencies: {AR: 0.625}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->type(), "falloff-SRI");
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiency("AR"), 0.625);
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiency("N2"), 1.0);
    const auto rate = std::dynamic_pointer_cast<SriRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->highRate().preExponentialFactor(), 7.91E+10);
    EXPECT_DOUBLE_EQ(rate->lowRate().preExponentialFactor(), 6.37E+14);
    EXPECT_DOUBLE_EQ(rate->lowRate().activationEnergy(), 56640);
}

TEST(Reaction, FalloffFromYaml2)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: H + CH2 (+ N2) <=> CH3 (+N2),"
        " type: falloff,"
        " high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0],"
        " low-P-rate-constant: {A: 1.04000E+26 cm^6/mol^2/s, b: -2.76, Ea: 1600},"
        " Troe: {A: 0.562, T3: 91, T1: 5836},"
        " source: somewhere}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->type(), "falloff-Troe");
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiency("N2"), 1.0);
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiency("H2O"), 0.0);
    const auto rate = std::dynamic_pointer_cast<TroeRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->highRate().preExponentialFactor(), 6e11);
    EXPECT_DOUBLE_EQ(rate->lowRate().preExponentialFactor(), 1.04e20);
    EXPECT_DOUBLE_EQ(rate->lowRate().activationEnergy(), 1600);
    vector<double> params;
    rate->getFalloffCoeffs(params);
    ASSERT_EQ(params.size(), (size_t) 3);
    EXPECT_DOUBLE_EQ(params[0], 0.562);
    EXPECT_DOUBLE_EQ(params[1], 91.0);
    EXPECT_EQ(R->input["source"].asString(), "somewhere");
}

TEST(Reaction, FalloffFromYaml3)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: HCN (+ M) <=> H + CN (+ M),"
        " type: falloff,"
        " low-P-rate-constant: {A: 3.57e+26, b: -2.6, Ea: 1.249e+05},"
        " high-P-rate-constant: {A: 8.3e+17, b: -0.93, Ea: 1.238e+05},"
        " Tsang: {A: 0.95, B: -1.0e-04},"
        " efficiencies: {CO2: 1.6, H2O: 5.0, N2: 1.0, N2O: 5.0},"
        " source: ARL-TR-5088}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->type(), "falloff-Tsang");
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiency("N2"), 1.0);
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiency("H2O"), 5.0);
    const auto rate = std::dynamic_pointer_cast<TsangRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->highRate().preExponentialFactor(), 8.3e17);
    EXPECT_DOUBLE_EQ(rate->lowRate().preExponentialFactor(), 3.57e26);
    EXPECT_DOUBLE_EQ(rate->highRate().activationEnergy(), 123800.0);
    EXPECT_DOUBLE_EQ(rate->lowRate().activationEnergy(), 124900.0);
    vector<double> params(2);
    rate->getFalloffCoeffs(params);
    EXPECT_DOUBLE_EQ(params[0], 0.95);
    EXPECT_DOUBLE_EQ(params[1], -0.0001);
    EXPECT_EQ(R->input["source"].asString(), "ARL-TR-5088");
}

TEST(Reaction, ChemicallyActivatedFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH3 + OH (+M) <=> CH2O + H2 (+M),"
        " units: {length: cm, quantity: mol},"
        " type: chemically-activated,"
        " high-P-rate-constant: [5.88E-14, 6.721, -3022.227],"
        " low-P-rate-constant: [282320.078, 1.46878, -3270.56495]}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->type(), "chemically-activated-Lindemann");
    const auto& rate = std::dynamic_pointer_cast<LindemannRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->highRate().preExponentialFactor(), 5.88e-14);
    EXPECT_DOUBLE_EQ(rate->lowRate().preExponentialFactor(), 2.82320078e2);
    EXPECT_EQ(rate->nParameters(), (size_t) 0);
}

TEST(Reaction, PlogFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "equation: 'H + CH4 <=> H2 + CH3'\n"
        "units: {pressure: atm}\n"
        "type: pressure-dependent-Arrhenius\n"
        "rate-constants:\n"
        "- {P: 0.039474, A: 2.720000e+09 cm^3/mol/s, b: 1.2, Ea: 6834.0}\n"
        "- {P: 1.0 atm, A: 1.260000e+20, b: -1.83, Ea: 15003.0}\n"
        "- {P: 1.0 atm, A: 1.230000e+04, b: 2.68, Ea: 6335.0}\n"
        "- {P: 1.01325 MPa, A: 1.680000e+16, b: -0.6, Ea: 14754.0}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_FALSE(R->usesThirdBody());
    const auto& rateMap = std::dynamic_pointer_cast<PlogRate>(R->rate())->getRates();
    vector<pair<double, ArrheniusRate>> rates(rateMap.begin(), rateMap.end());
    EXPECT_EQ(rates.size(), (size_t) 4);
    EXPECT_NEAR(rates[0].first, 0.039474 * OneAtm, 1e-6);
    EXPECT_NEAR(rates[2].first, OneAtm, 1e-6);
    EXPECT_NEAR(rates[3].first, 10 * OneAtm, 1e-6);
    EXPECT_DOUBLE_EQ(rates[0].second.preExponentialFactor(), 2.72e6);
    EXPECT_DOUBLE_EQ(rates[3].second.preExponentialFactor(), 1.68e16);
    EXPECT_DOUBLE_EQ(rates[3].second.temperatureExponent(), -0.6);
}

TEST(Reaction, ChebyshevFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "equation: 'CH4 <=> CH3 + H'\n"
        "type: Chebyshev\n"
        "temperature-range: [290, 3000]\n"
        "pressure-range: [0.0098692326671601278 atm, 98.692326671601279 atm]\n"
        "data: [[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],\n"
        "       [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],\n"
        "       [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],\n"
        "       [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],\n"
        "       [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],\n"
        "       [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]]\n");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_FALSE(R->usesThirdBody());
    EXPECT_EQ(R->reactants.size(), (size_t) 1);
    const auto& rate = std::dynamic_pointer_cast<ChebyshevRate>(R->rate());
    double T = 1800;
    double P = 2e6;
    EXPECT_EQ(rate->data().nRows(), (size_t) 6);
    EXPECT_EQ(rate->data().nColumns(), (size_t) 4);
    EXPECT_DOUBLE_EQ(rate->Tmax(), 3000);
    EXPECT_DOUBLE_EQ(rate->Pmin(), 1000);
    EXPECT_NEAR(rate->eval(T, P), 130512.2773948636, 2e-9);
}

TEST(Reaction, BlowersMaselFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: O + H2 <=> H + OH,"
        " type: Blowers-Masel,"
        " rate-constant: [-3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_FALSE(R->usesThirdBody());
    EXPECT_EQ(R->reactants.at("H2"), 1);
    EXPECT_EQ(R->products.at("OH"), 1);

    double E_intrinsic = 6260 * 4184.0; // J/kmol
    double H_big = 5 * E_intrinsic;
    double H_small = -5 * E_intrinsic;
    double H_mid = 4 * E_intrinsic;
    double w = 1e9 * 4184.0; // J/kmol
    double vp = 2 * w * ((w + E_intrinsic) / (w - E_intrinsic));
    double Ea = (w + H_mid / 2) * (vp - 2 * w + H_mid) * (vp - 2 * w + H_mid)
        / (vp * vp - 4 * w * w + H_mid * H_mid );
    const auto& rate = std::dynamic_pointer_cast<BlowersMaselRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), -38.7);
    EXPECT_DOUBLE_EQ(rate->bondEnergy(), w);
    rate->setDeltaH(H_big);
    EXPECT_DOUBLE_EQ(rate->activationEnergy(), H_big);
    rate->setDeltaH(H_small);
    EXPECT_DOUBLE_EQ(rate->activationEnergy(), 0);
    rate->setDeltaH(H_mid);
    EXPECT_NEAR(rate->activationEnergy(), Ea, 5e-4);
    EXPECT_TRUE(rate->allowNegativePreExponentialFactor());
    EXPECT_FALSE(R->allow_negative_orders);
}

TEST(Reaction, ThreeBodyBlowersMaselFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH2 + O2 <=> CH2(S) + O2,"
        " type: three-body-Blowers-Masel,"
        " rate-constant: [3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol]}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->type(), "three-body-Blowers-Masel");
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->thirdBody()->name(), "O2");

    const auto& rate = std::dynamic_pointer_cast<BlowersMaselRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 38.7);

    AnyMap input = R->parameters(false);
    EXPECT_EQ(input.getString("type", ""), "three-body-Blowers-Masel");
}

TEST(Reaction, InvalidBlowersMaselFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: O + H2 <=> H + OH,"
        " type: three-body-Blowers-Masel,"
        " rate-constant: [3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol]}");

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, TwoTempPlasmaFromYaml)
{
    auto sol = newSolution("ET_test.yaml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: E + O2 + O2 => O2^- + O2,"
        " type: two-temperature-plasma,"
        " rate-constant: [1.523e+27 cm^6/mol^2/s, -1.0, -100 K, 700 K]}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_TRUE(R->usesThirdBody());
    EXPECT_EQ(R->thirdBody()->name(), "O2");
    EXPECT_EQ(R->reactants.at("O2"), 1);
    EXPECT_EQ(R->reactants.at("E"), 1);
    EXPECT_EQ(R->products.at("O2^-"), 1);

    const auto rate = std::dynamic_pointer_cast<TwoTempPlasmaRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 1.523e21);
    EXPECT_DOUBLE_EQ(rate->temperatureExponent(), -1.0);
    EXPECT_DOUBLE_EQ(rate->activationEnergy(), -100 * GasConstant);
    EXPECT_DOUBLE_EQ(rate->activationElectronEnergy(), 700 * GasConstant);
}

TEST(Reaction, ElectronCollisionPlasmaFromYaml)
{
    auto sol = newSolution("oxygen-plasma.yaml", "", "none");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: O2 + E => E + O2,"
        " type: electron-collision-plasma,"
        " energy-levels: [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0],"
        " cross-sections: [0.0, 5.97e-20, 6.45e-20, 6.74e-20, 6.93e-20, 7.2e-20, "
        " 7.52e-20, 7.86e-20, 8.21e-20, 8.49e-20, 8.8e-20]}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.at("O2"), 1);
    EXPECT_EQ(R->reactants.at("E"), 1);
    EXPECT_EQ(R->products.at("O2"), 1);
    EXPECT_EQ(R->products.at("E"), 1);

    const auto rate = std::dynamic_pointer_cast<ElectronCollisionPlasmaRate>(R->rate());

    for (size_t k = 0; k < rate->energyLevels().size(); k++) {
        EXPECT_DOUBLE_EQ(rate->energyLevels()[k], rxn["energy-levels"].asVector<double>()[k]);
        EXPECT_DOUBLE_EQ(rate->crossSections()[k], rxn["cross-sections"].asVector<double>()[k]);
    }
}

TEST(Reaction, PythonExtensibleRate)
{
    #ifdef CT_SKIP_PYTHON // Possibly set via test/SConscript
    GTEST_SKIP();
    #endif
    auto sol = newSolution("extensible-reactions.yaml");
    auto R = sol->kinetics()->reaction(0);
    EXPECT_EQ(R->type(), "square-rate");
    auto rate = R->rate();
    EXPECT_EQ(rate->type(), "square-rate");
    vector<double> kf(sol->kinetics()->nReactions());
    sol->kinetics()->getFwdRateConstants(kf.data());
    EXPECT_DOUBLE_EQ(kf[0], 3.14 * 300 * 300);
}

TEST(Kinetics, BulkKineticsFromYaml1)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode = infile["phases"].getMapWhere("name", "simple-kinetics");
    shared_ptr<ThermoPhase> thermo = newThermo(phaseNode, infile);
    auto kin = newKinetics({thermo}, phaseNode, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 2);
    const auto& R = kin->reaction(0);
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->id, "NOx-R1");
    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 2.7e10);
}

TEST(Kinetics, BulkKineticsFromYaml2)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode = infile["phases"].getMapWhere("name", "remote-kinetics");
    shared_ptr<ThermoPhase> thermo = newThermo(phaseNode, infile);
    auto kin = newKinetics({thermo}, phaseNode, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 3);
}

TEST(Kinetics, EfficienciesFromYaml)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode1 = infile["phases"].getMapWhere("name", "efficiency-error");
    shared_ptr<ThermoPhase> thermo1 = newThermo(phaseNode1, infile);
    // Reaction with efficiency defined for undeclared species "AR"
    EXPECT_THROW(newKinetics({thermo1}, phaseNode1, infile), CanteraError);

    auto& phaseNode2 = infile["phases"].getMapWhere("name", "efficiency-skip");
    shared_ptr<ThermoPhase> thermo2 = newThermo(phaseNode2, infile);
    auto kin = newKinetics({thermo2}, phaseNode2, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 1);
}

TEST(Kinetics, InterfaceKineticsFromYaml)
{
    auto soln = newInterface("surface-phases.yaml", "Pt-surf");
    auto kin = soln->kinetics();
    EXPECT_EQ(kin->nReactions(), (size_t) 3);
    EXPECT_EQ(kin->nTotalSpecies(), (size_t) 6);

    auto R1 = kin->reaction(0);
    EXPECT_DOUBLE_EQ(R1->orders["Pt(s)"], 1.0);
    const auto rate1 = std::dynamic_pointer_cast<InterfaceArrheniusRate>(R1->rate());
    EXPECT_DOUBLE_EQ(rate1->preExponentialFactor(), 4.4579e7);

    auto R2 = kin->reaction(1);
    const auto rate2 = std::dynamic_pointer_cast<InterfaceArrheniusRate>(R2->rate());
    EXPECT_DOUBLE_EQ(rate2->preExponentialFactor(), 3.7e20);
    AnyMap coverage_deps;
    rate2->getCoverageDependencies(coverage_deps);
    coverage_deps.applyUnits();
    auto& cov_map = coverage_deps["H(s)"].as<AnyMap>();
    EXPECT_DOUBLE_EQ(cov_map["E"].asDouble(), -6e6);


    auto R3 = kin->reaction(2);
    EXPECT_TRUE(std::dynamic_pointer_cast<StickingArrheniusRate>(R3->rate()));

    auto soln2 = newInterface("surface-phases.yaml", "Pt-multi-sites");
    auto kin2 = soln2->kinetics();
    auto R4 = kin2->reaction(3);
    vector<double> coeffs{1.0e6, 3.0e6, -7.0e7, 5.0e6};
    const auto rate4 = std::dynamic_pointer_cast<InterfaceArrheniusRate>(R4->rate());
    AnyMap coverage_deps2;
    rate4->getCoverageDependencies(coverage_deps2);
    coverage_deps2.applyUnits();
    auto& cov_map2 = coverage_deps2["O(s)"].as<AnyMap>();
    auto& poly_dep = cov_map2["E"].asVector<AnyValue>();
    for (size_t i = 0; i < poly_dep.size(); i++) {
        EXPECT_DOUBLE_EQ(poly_dep[i].asDouble(), coeffs[i]);
    }
}

TEST(Kinetics, BMInterfaceKineticsFromYaml)
{
    auto soln = newInterface("blowers-masel.yaml", "Pt_surf");
    auto kin = soln->kinetics();
    EXPECT_EQ(kin->nReactions(), (size_t) 6);
    EXPECT_EQ(kin->nTotalSpecies(), (size_t) 14);

    auto R1 = kin->reaction(5);
    EXPECT_DOUBLE_EQ(R1->orders["PT(s)"], 1.0);
    const auto rate = std::dynamic_pointer_cast<InterfaceBlowersMaselRate>(R1->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 4.4579e7);

    auto R2 = kin->reaction(0);
    const auto rate2 = std::dynamic_pointer_cast<InterfaceBlowersMaselRate>(R2->rate());
    EXPECT_DOUBLE_EQ(rate2->preExponentialFactor(), 3.7e20);
    AnyMap coverage_deps;
    rate2->getCoverageDependencies(coverage_deps);
    coverage_deps.applyUnits();
    auto& cov_map = coverage_deps["H(S)"].as<AnyMap>();
    EXPECT_DOUBLE_EQ(cov_map["E"].asDouble(), -6e6);

    auto R3 = kin->reaction(1);
    EXPECT_TRUE(std::dynamic_pointer_cast<StickingBlowersMaselRate>(R3->rate()));
}

TEST(Kinetics, ElectrochemFromYaml)
{
    auto soln = newInterface("surface-phases.yaml", "anode-surface");
    auto kin = soln->kinetics();
    soln->adjacent("graphite")->thermo()->setElectricPotential(0.4);
    vector<double> ropf(kin->nReactions()), ropr(kin->nReactions());
    kin->getFwdRatesOfProgress(ropf.data());
    kin->getRevRatesOfProgress(ropr.data());

    EXPECT_NEAR(ropf[0], 0.279762338, 1e-8);
    EXPECT_NEAR(ropr[0], 0.045559670, 1e-8);
}

TEST(KineticsFromYaml, NoKineticsModelOrReactionsField1)
{
    auto soln = newSolution("phase-reaction-spec1.yaml",
                            "nokinetics-noreactions");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "none");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, NoKineticsModelOrReactionsField2)
{
    auto soln = newSolution("phase-reaction-spec2.yaml",
                            "nokinetics-noreactions");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "none");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithReactionsNone1)
{
    auto soln = newSolution("phase-reaction-spec1.yaml",
                            "kinetics-reactions-none");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "bulk");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithReactionsNone2)
{
    auto soln = newSolution("phase-reaction-spec2.yaml",
                            "kinetics-reactions-none");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "bulk");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsSection1)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "kinetics-no-reaction-section1"),
                 InputFileError);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsSection2)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "kinetics-no-reaction-section2"),
                 InputFileError);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsSection3)
{
    auto soln = newSolution("phase-reaction-spec1.yaml",
                            "kinetics-no-reaction-section3");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "surface");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsSection4)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "kinetics-no-reaction-section4"),
                 InputFileError);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsField)
{
    auto soln = newSolution("phase-reaction-spec2.yaml",
                            "kinetics-noreactions");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "bulk");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 1);
}

TEST(KineticsFromYaml, ReactionsFieldWithoutKineticsModel1)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "nokinetics-reactions"),
                 InputFileError);
}

TEST(KineticsFromYaml, ReactionsFieldWithoutKineticsModel2)
{
    EXPECT_THROW(newSolution("phase-reaction-spec2.yaml",
                             "nokinetics-reactions"),
                 InputFileError);
}

TEST(KineticsFromYaml, InvalidExtension)
{
    AnyMap input = AnyMap::fromYamlFile("h2o2.yaml");
    newSolution(input["phases"].asVector<AnyMap>()[0], input);
    vector<AnyMap> extensions(1);
    extensions[0]["type"] = "nonexistent";
    extensions[0]["name"] = "fake";
    input["extensions"] = extensions;
    EXPECT_THROW(newSolution(input["phases"].asVector<AnyMap>()[0], input),
                 CanteraError);
}

class ReactionToYaml : public testing::Test
{
public:
    void duplicateReaction(size_t i) {
        auto kin = soln->kinetics();
        iOld = i;
        // Exclude the raw input data, to make sure this test actually relies on
        // the fields being populated by the Reaction types.
        AnyMap rdata1 = kin->reaction(iOld)->parameters(false);
        AnyMap rdata2 = AnyMap::fromYamlString(rdata1.toYamlString());
        duplicate = newReaction(rdata2, *kin);
        kin->addReaction(duplicate);
        iNew = kin->nReactions() - 1;
    }

    void compareReactions() {
        auto kin = soln->kinetics();
        EXPECT_EQ(kin->reaction(iOld)->equation(),
                  kin->reaction(iNew)->equation());
        EXPECT_EQ(kin->isReversible(iOld), kin->isReversible(iNew));

        vector<double> kf(kin->nReactions()), kr(kin->nReactions());
        vector<double> ropf(kin->nReactions()), ropr(kin->nReactions());
        kin->getFwdRateConstants(kf.data());
        kin->getRevRateConstants(kr.data());
        kin->getFwdRatesOfProgress(ropf.data());
        kin->getRevRatesOfProgress(ropr.data());
        EXPECT_NEAR(kf[iOld], kf[iNew], 1e-13 * kf[iOld]);
        EXPECT_NEAR(kr[iOld], kr[iNew], 1e-13 * kr[iOld]);
        EXPECT_NEAR(ropf[iOld], ropf[iNew], 1e-13 * ropf[iOld]);
        EXPECT_NEAR(ropr[iOld], ropr[iNew], 1e-13 * ropr[iOld]);
    }

    shared_ptr<Solution> soln;
    shared_ptr<Reaction> duplicate;

    size_t iOld;
    size_t iNew;
};

TEST_F(ReactionToYaml, elementary)
{
    soln = newSolution("h2o2.yaml", "", "none");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, O2:0.5, O:1e-8, OH:3e-8");
    duplicateReaction(2);
    EXPECT_EQ(duplicate->type(), "Arrhenius");
    EXPECT_EQ(duplicate->rate()->type(), "Arrhenius");
    compareReactions();
}

TEST_F(ReactionToYaml, threeBody)
{
    soln = newSolution("h2o2.yaml", "", "none");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, O2:0.5, O:1e-8, OH:3e-8, H:2e-7");
    duplicateReaction(1);
    EXPECT_TRUE(std::dynamic_pointer_cast<ArrheniusRate>(duplicate->rate()));
    compareReactions();
}

TEST_F(ReactionToYaml, TroeFalloff)
{
    soln = newSolution("h2o2.yaml", "", "none");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, O2:0.5, H2O2:1e-8, OH:3e-8");
    duplicateReaction(21);
    auto rate = std::dynamic_pointer_cast<TroeRate>(duplicate->rate());
    EXPECT_TRUE(rate);
    EXPECT_FALSE(rate->chemicallyActivated());
    compareReactions();
}

TEST_F(ReactionToYaml, SriFalloff)
{
    soln = newSolution("sri-falloff.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "R1A: 0.1, R1B:0.2, H: 0.2, R2:0.5");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<SriRate>(duplicate->rate()));
    compareReactions();
    duplicateReaction(1);
    compareReactions();
}

TEST_F(ReactionToYaml, TsangFalloff)
{
    soln = newSolution("tsang-falloff.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "NO:1.0, OH:1.0, H:1.0, CN:1.0");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<TsangRate>(duplicate->rate()));
    compareReactions();
    duplicateReaction(1);
    compareReactions();
}

TEST_F(ReactionToYaml, chemicallyActivated)
{
    soln = newSolution("chemically-activated-reaction.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, ch2o:0.1, ch3:1e-8, oh:3e-6");
    duplicateReaction(0);
    auto rate = std::dynamic_pointer_cast<TroeRate>(duplicate->rate());
    EXPECT_TRUE(rate);
    EXPECT_TRUE(rate->chemicallyActivated());
    compareReactions();
}

TEST_F(ReactionToYaml, pdepArrhenius)
{
    soln = newSolution("pdep-test.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "R2:1, H:0.1, P2A:2, P2B:0.3");
    duplicateReaction(1);
    EXPECT_TRUE(std::dynamic_pointer_cast<PlogRate>(duplicate->rate()));
    compareReactions();
    soln->thermo()->setState_TPY(1100, 1e3, "R2:1, H:0.2, P2A:2, P2B:0.3");
    compareReactions();
}

TEST_F(ReactionToYaml, Chebyshev)
{
    soln = newSolution("pdep-test.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "R6:1, P6A:2, P6B:0.3");
    duplicateReaction(5);
    EXPECT_TRUE(std::dynamic_pointer_cast<ChebyshevRate>(duplicate->rate()));
    compareReactions();
}

TEST_F(ReactionToYaml, surface)
{
    auto gas = newSolution("diamond.yaml", "gas");
    auto solid = newSolution("diamond.yaml", "diamond");
    soln = newSolution("diamond.yaml", "diamond_100", "none", {gas, solid});
    auto surf = std::dynamic_pointer_cast<SurfPhase>(soln->thermo());
    surf->setCoveragesByName("c6HH:0.1, c6H*:0.6, c6**:0.1");
    gas->thermo()->setMassFractionsByName("H2:0.7, CH4:0.3");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<InterfaceArrheniusRate>(duplicate->rate()));
    compareReactions();
}

TEST_F(ReactionToYaml, electrochemical)
{
    soln = newSolution("sofc.yaml", "tpb");
    auto oxide_surf = soln->adjacent("oxide_surface");
    auto oxide_bulk = oxide_surf->adjacent("oxide_bulk");

    // There should be only one underlying 'gas' object
    ASSERT_EQ(oxide_surf->adjacent("gas").get(),
              soln->adjacent("metal_surface")->adjacent("gas").get());

    auto ox_surf = std::dynamic_pointer_cast<SurfPhase>(oxide_surf->thermo());
    oxide_bulk->thermo()->setElectricPotential(-3.4);
    oxide_surf->thermo()->setElectricPotential(-3.4);
    ox_surf->setCoveragesByName("O''(ox):0.2, OH'(ox):0.3, H2O(ox):0.5");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<InterfaceArrheniusRate>(duplicate->rate()));
    compareReactions();
    compareReactions();
}

TEST_F(ReactionToYaml, unconvertible1)
{
    auto rate = make_shared<ArrheniusRate>(1e5, -1.0, 12.5);
    Reaction R({{"H2", 1}, {"OH", 1}}, {{"H2O", 1}, {"H", 1}}, rate);
    AnyMap params = R.parameters();
    UnitSystem U{"g", "cm", "mol"};
    params.setUnits(U);
    EXPECT_THROW(params.applyUnits(), CanteraError);
}

TEST_F(ReactionToYaml, unconvertible2)
{
    Array2D coeffs(2, 2, 1.0);
    auto rate = make_shared<ChebyshevRate>(273., 3000., 1.e2, 1.e7, coeffs);
    Reaction R({{"H2", 1}, {"OH", 1}}, {{"H2O", 1}, {"H", 1}}, rate);
    UnitSystem U{"g", "cm", "mol"};
    AnyMap params = R.parameters();
    params.setUnits(U);
    EXPECT_THROW(params.applyUnits(), CanteraError);
}

TEST_F(ReactionToYaml, unconvertible3)
{
    Reaction R(
        {{"H2", 1}, {"OH", 1}}, {{"H2O", 1}, {"H", 1}},
        make_shared<TroeRate>(
            ArrheniusRate(1e5, -1.0, 12.5), ArrheniusRate(1e5, -1.0, 12.5),
            vector<double>{0.562, 91.0, 5836.0, 8552.0}));
    AnyMap params = R.parameters();
    UnitSystem U{"g", "cm", "mol"};
    params.setUnits(U);
    EXPECT_THROW(params.applyUnits(), CanteraError);
}

TEST_F(ReactionToYaml, BlowersMaselRate)
{
    soln = newSolution("blowers-masel.yaml", "gas");
    soln->thermo()->setState_TPY(1100, 0.1 * OneAtm, "O:0.01, H2:0.8, O2:0.19");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<BlowersMaselRate>(duplicate->rate()));
    compareReactions();
}

TEST_F(ReactionToYaml, BlowersMaselInterface)
{
    auto gas = newSolution("blowers-masel.yaml", "gas");
    soln = newSolution("blowers-masel.yaml", "Pt_surf", "none", {gas});
    gas->thermo()->setState_TPY(1100, 0.1 * OneAtm, "O:0.01, H2:0.8, O2:0.19");
    soln->thermo()->setState_TP(1100, 0.1 * OneAtm);
    auto surf = std::dynamic_pointer_cast<SurfPhase>(soln->thermo());
    surf->setCoveragesByName("H(S):0.1, PT(S):0.8, H2O(S):0.1");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<InterfaceBlowersMaselRate>(duplicate->rate()));
    compareReactions();
}
