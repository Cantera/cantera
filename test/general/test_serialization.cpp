// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "cantera/base/YamlWriter.h"
#include "cantera/thermo.h"
#include "cantera/base/Solution.h"
#include "cantera/kinetics.h"

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

TEST(YamlWriter, reactions)
{
    auto original = newSolution("h2o2.yaml");
    YamlWriter writer;
    writer.addPhase(original);
    writer.toYamlFile("generated-h2o2.yaml");
    auto duplicate = newSolution("generated-h2o2.yaml");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector_fp kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_DOUBLE_EQ(kf1[i], kf2[i]) << "for reaction i = " << i;
    }

    AnyMap m = AnyMap::fromYamlFile("generated-h2o2.yaml");
    auto& reactions = m["reactions"].asVector<AnyMap>();
}

TEST(YamlWriter, multipleReactionSections)
{
    auto original1 = newSolution("h2o2.yaml");
    auto original2 = newSolution("h2o2.yaml");
    auto original3 = newSolution("h2o2.yaml");
    // this phase will require its own "reactions" section
    auto R = original3->kinetics()->reaction(3);
    R->duplicate = true;
    original3->kinetics()->addReaction(R);
    original2->setName("ohmech2");
    original3->setName("ohmech3");
    YamlWriter writer;
    writer.addPhase(original1);
    writer.addPhase(original2);
    writer.addPhase(original3);
    writer.toYamlFile("generated-multi-rxn-secs.yaml");

    auto duplicate1 = newSolution("generated-multi-rxn-secs.yaml", "ohmech");
    auto duplicate2 = newSolution("generated-multi-rxn-secs.yaml", "ohmech2");
    auto duplicate3 = newSolution("generated-multi-rxn-secs.yaml", "ohmech3");
    auto kin1 = duplicate1->kinetics();
    auto kin2 = duplicate2->kinetics();
    auto kin3 = duplicate3->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    ASSERT_EQ(kin2->nReactions() + 1, kin3->nReactions());
    ASSERT_EQ(kin2->reactionString(3),
              kin3->reactionString(kin3->nReactions() - 1));
}
