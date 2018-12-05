#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/IdealGasMix.h"

using namespace Cantera;

TEST(Reaction, ElementaryFromYaml)
{
    IdealGasMix gas("gri30.xml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
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
        " rate: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
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
        " rate: [1.20000E+17, -1, 0]}");

    UnitSystem U;
    EXPECT_THROW(newReaction(rxn, gas, U), CanteraError);
}
