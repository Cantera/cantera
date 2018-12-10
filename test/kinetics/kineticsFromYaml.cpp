#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/IdealGasMix.h"

using namespace Cantera;

TEST(Reaction, ElementaryFromYaml)
{
    // @TODO: Use of XML input files in these tests of the YAML format needs to
    // be eliminated before we can deprecate the XML format.
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    UnitSystem U;
    auto R = newReaction(rxn, gas, U);
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->reaction_type, ELEMENTARY_RXN);

    auto ER = dynamic_cast<ElementaryReaction&>(*R);
    EXPECT_DOUBLE_EQ(ER.rate.preExponentialFactor(), -2.7e10);
    EXPECT_DOUBLE_EQ(ER.rate.activationEnergy_R(), 355 / GasConst_cal_mol_K);
    EXPECT_TRUE(ER.allow_negative_pre_exponential_factor);
    EXPECT_FALSE(ER.allow_negative_orders);
}

TEST(Reaction, ThreeBodyFromYaml1)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O + M <=> O2 + M,"
        " type: three-body,"
        " rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
        " efficiencies: {AR: 0.83, H2O: 5}}");

    UnitSystem U;
    auto R = newReaction(rxn, gas, U);
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);

    auto TBR = dynamic_cast<ThreeBodyReaction&>(*R);
    EXPECT_DOUBLE_EQ(TBR.rate.preExponentialFactor(), 1.2e11);
    EXPECT_DOUBLE_EQ(TBR.third_body.efficiencies["H2O"], 5.0);
    EXPECT_DOUBLE_EQ(TBR.third_body.default_efficiency, 1.0);
}

TEST(Reaction, ThreeBodyFromYaml2)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O <=> O2," // Missing "M" on each side of the equation
        " type: three-body,"
        " rate-constant: [1.20000E+17, -1, 0]}");

    UnitSystem U;
    EXPECT_THROW(newReaction(rxn, gas, U), CanteraError);
}

TEST(Reaction, FalloffFromYaml1)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N2O (+M) <=> N2 + O (+ M),"
        " type: falloff,"
        " high-P-rate-constant: [7.91000E+10, 0, 56020],"
        " low-P-rate-constant: [6.37000E+14, 0, 56640],"
        " SRI: {A: 1.1, B: 700.0, C: 1234.0, D: 56.0, E: 0.7},"
        " efficiencies: {AR: 0.625}}");

    UnitSystem U;
    auto R = newReaction(rxn, gas, U);
    auto FR = dynamic_cast<FalloffReaction&>(*R);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("AR"), 0.625);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("N2"), 1.0);
}

TEST(Reaction, FalloffFromYaml2)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: H + CH2 (+ N2) <=> CH3 (+N2),"
        " type: falloff,"
        " high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0],"
        " low-P-rate-constant: [1.04000E+26 cm^6/mol^2/s, -2.76, 1600],"
        " Troe: {A: 0.562, T3: 91, T1: 5836}}");

    UnitSystem U;
    auto R = newReaction(rxn, gas, U);
    auto FR = dynamic_cast<FalloffReaction&>(*R);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("N2"), 1.0);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("H2O"), 0.0);
    EXPECT_DOUBLE_EQ(FR.high_rate.preExponentialFactor(), 6e11);
    EXPECT_DOUBLE_EQ(FR.low_rate.preExponentialFactor(), 1.04e20);
    vector_fp params(4);
    FR.falloff->getParameters(params.data());
    EXPECT_DOUBLE_EQ(params[0], 0.562);
    EXPECT_DOUBLE_EQ(params[1], 91.0);
    EXPECT_DOUBLE_EQ(params[3], 0.0);
}

TEST(Reaction, ChemicallyActivatedFromYaml)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH3 + OH (+M) <=> CH2O + H2 (+M),"
        " type: chemically-activated,"
        " high-P-rate-constant: [5.88E-14, 6.721, -3022.227],"
        " low-P-rate-constant: [282320.078, 1.46878, -3270.56495]}");

    UnitSystem U;
    U.setDefaults({"cm", "mol"});
    auto R = newReaction(rxn, gas, U);
    auto CAR = dynamic_cast<ChemicallyActivatedReaction&>(*R);
    EXPECT_DOUBLE_EQ(CAR.high_rate.preExponentialFactor(), 5.88e-14);
    EXPECT_DOUBLE_EQ(CAR.low_rate.preExponentialFactor(), 2.82320078e2);
    EXPECT_EQ(CAR.falloff->nParameters(), (size_t) 0);
}

TEST(Reaction, PlogFromYaml)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "equation: 'H + CH4 <=> H2 + CH3'\n"
        "type: pressure-dependent-Arrhenius\n"
        "rate-constants:\n"
        "- [0.039474, [2.720000e+09 cm^3/mol/s, 1.2, 6834.0]]\n"
        "- [1.0 atm, [1.260000e+20, -1.83, 15003.0]]\n"
        "- [1.0 atm, [1.230000e+04, 2.68, 6335.0]]\n"
        "- [1.01325 MPa, [1.680000e+16, -0.6, 14754.0]]");

    UnitSystem U({"atm"});
    auto R = newReaction(rxn, gas, U);
    auto PR = dynamic_cast<PlogReaction&>(*R);
    const auto& rates = PR.rate.rates();
    EXPECT_EQ(rates.size(), (size_t) 4);
    EXPECT_NEAR(rates[0].first, 0.039474 * OneAtm, 1e-6);
    EXPECT_NEAR(rates[2].first, OneAtm, 1e-6);
    EXPECT_NEAR(rates[3].first, 10 * OneAtm, 1e-6);
    EXPECT_DOUBLE_EQ(rates[0].second.preExponentialFactor(), 2.72e6);
    EXPECT_DOUBLE_EQ(rates[3].second.preExponentialFactor(), 1.68e16);
}
