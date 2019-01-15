#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Elements.h"

using namespace Cantera;

TEST(ThermoFromYaml, simpleIdealGas)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("simple"), infile);
    EXPECT_EQ(thermo->nSpecies(), (size_t) 3);
    EXPECT_DOUBLE_EQ(thermo->density(), 7.031763356741983);
    EXPECT_DOUBLE_EQ(thermo->cp_mass(), 1037.7632754708304);
}

TEST(ThermoFromYaml, failDuplicateSpecies)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    EXPECT_THROW(newPhase(*phaseNodes.at("duplicate-species"), infile),
                 CanteraError);
}

TEST(ThermoFromYaml, elementOverride)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("element-override"), infile);
    EXPECT_EQ(thermo->nElements(), (size_t) 3);
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(0), getElementWeight("N"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(1), getElementWeight("O"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(2), 36);
}

TEST(ThermoFromYaml, elementFromDifferentFile)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("element-remote"), infile);
    EXPECT_EQ(thermo->nElements(), (size_t) 3);
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(0), getElementWeight("N"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(1), getElementWeight("O"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(2), 38);
}

TEST(ThermoFromYaml, speciesFromDifferentFile)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("species-remote"), infile);
    EXPECT_EQ(thermo->nElements(), (size_t) 2);
    EXPECT_EQ(thermo->nSpecies(), (size_t) 4);
    EXPECT_EQ(thermo->species(0)->composition["O"], 2);
    EXPECT_EQ(thermo->species(3)->composition["O"], 1);
    EXPECT_EQ(thermo->species(2)->name, "NO2");
    EXPECT_DOUBLE_EQ(thermo->moleFraction(3), 0.3);
}

TEST(ThermoFromYaml, speciesAll)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("species-all"), infile);
    EXPECT_EQ(thermo->nElements(), (size_t) 3);
    EXPECT_EQ(thermo->nSpecies(), (size_t) 6);
    EXPECT_EQ(thermo->species(1)->name, "NO");
    EXPECT_EQ(thermo->species(2)->name, "N2");
}

TEST(ThermoFromYaml, StoichSubstance1)
{
    AnyMap infile = AnyMap::fromYamlFile("thermo-models.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("NaCl(s)"), infile);
    EXPECT_EQ(thermo->type(), "StoichSubstance");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    EXPECT_EQ(thermo->nElements(), (size_t) 2);
    EXPECT_DOUBLE_EQ(thermo->density(), 2165.0);
    EXPECT_DOUBLE_EQ(thermo->cp_mass(), 864.8437519457644); // Regression test based on XML
}

TEST(ThermoFromYaml, StoichSubstance2)
{
    AnyMap infile = AnyMap::fromYamlFile("thermo-models.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("KCl(s)"), infile);
    EXPECT_EQ(thermo->type(), "StoichSubstance");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    EXPECT_EQ(thermo->nElements(), (size_t) 2);
    EXPECT_NEAR(thermo->density(), 1980, 0.1);
}

TEST(ThermoFromYaml, WaterSSTP)
{
    AnyMap infile = AnyMap::fromYamlFile("thermo-models.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("liquid-water"), infile);
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    thermo->setState_TP(350, 2*OneAtm);
    // Regression tests based on XML
    EXPECT_NEAR(thermo->density(), 973.7736331, 1e-6);
    EXPECT_NEAR(thermo->enthalpy_mass(), -15649442.2898854, 1e-6);
}

TEST(ThermoFromYaml, FixedChemPot)
{
    AnyMap infile = AnyMap::fromYamlFile("thermo-models.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("Li-fixed"), infile);
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    double mu;
    thermo->getChemPotentials(&mu);
    EXPECT_DOUBLE_EQ(mu, -2.3e7);
}

TEST(ThermoFromYaml, Margules)
{
    AnyMap infile = AnyMap::fromYamlFile("thermo-models.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("molten-salt-Margules"), infile);
    EXPECT_EQ(thermo->type(), "Margules");

    // Regression test based on LiKCl_liquid.xml
    EXPECT_NEAR(thermo->density(), 2042.1165603245981, 1e-9);
    EXPECT_NEAR(thermo->gibbs_mass(), -9682981.421693124, 1e-5);
    EXPECT_NEAR(thermo->cp_mole(), 67478.48085733457, 1e-8);
}

TEST(ThermoFromYaml, IdealMolalSoln)
{
    AnyMap infile = AnyMap::fromYamlFile("thermo-models.yaml");
    auto phaseNodes = infile["phases"].asMap("name");
    auto thermo = newPhase(*phaseNodes.at("ideal-molal-aqueous"), infile);
    EXPECT_EQ(thermo->type(), "IdealMolalSoln");

    EXPECT_NEAR(thermo->enthalpy_mole(), 0.013282, 1e-6);
    EXPECT_NEAR(thermo->gibbs_mole(), -3.8986e7, 1e3);
    EXPECT_NEAR(thermo->density(), 12.058, 1e-3);
}
