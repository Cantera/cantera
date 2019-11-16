// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "cantera/base/YamlWriter.h"
#include "cantera/thermo.h"
#include "cantera/base/Solution.h"

using namespace Cantera;

TEST(YamlWriter, thermoDef)
{
    auto original = newSolution("ideal-gas.yaml", "simple");
    YamlWriter writer;
    writer.addPhase(original);
    writer.toYamlFile("generated-simple.yaml");

    auto duplicate = newSolution("generated-simple.yaml", "simple");

    auto thermo1 = original->thermo();
    auto thermo2 = duplicate->thermo();

    EXPECT_EQ(thermo1->type(), thermo2->type());
    EXPECT_EQ(thermo1->speciesNames(), thermo2->speciesNames());
    EXPECT_DOUBLE_EQ(thermo1->pressure(), thermo2->pressure());
    EXPECT_DOUBLE_EQ(thermo1->enthalpy_mole(), thermo2->enthalpy_mole());
}

TEST(YamlWriter, sharedSpecies)
{
    auto original1 = newSolution("ideal-gas.yaml", "simple");
    auto original2 = newSolution("ideal-gas.yaml", "species-remote");

    YamlWriter writer;
    writer.addPhase(original1);
    writer.addPhase(original2);
    writer.toYamlFile("generated-shared-species.yaml");

    auto duplicate = newSolution("generated-shared-species.yaml", "species-remote");
    auto thermo1 = original2->thermo();
    auto thermo2 = duplicate->thermo();

    EXPECT_EQ(thermo1->type(), thermo2->type());
    EXPECT_EQ(thermo1->speciesNames(), thermo2->speciesNames());
    EXPECT_DOUBLE_EQ(thermo1->pressure(), thermo2->pressure());
    EXPECT_DOUBLE_EQ(thermo1->enthalpy_mole(), thermo2->enthalpy_mole());
}

TEST(YamlWriter, duplicateName)
{
    auto original1 = newSolution("ideal-gas.yaml", "simple");
    auto original2 = newSolution("ideal-gas.yaml", "simple");
    YamlWriter writer;
    writer.addPhase(original1);
    EXPECT_THROW(writer.addPhase(original2), CanteraError);
}
