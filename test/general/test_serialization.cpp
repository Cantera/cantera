// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "cantera/base/YamlWriter.h"
#include "cantera/thermo.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/Solution.h"
#include "cantera/kinetics.h"
#include "cantera/transport/TransportData.h"

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

TEST(YamlWriter, userDefinedFields)
{
    auto original = newSolution("ideal-gas.yaml", "simple");
    YamlWriter writer;
    writer.addPhase(original);

    AnyMap input1 = AnyMap::fromYamlString(writer.toYamlString());
    auto thermo1 = newPhase(input1["phases"].getMapWhere("name", "simple"),
                            input1);

    // user-defined fields should be in place
    EXPECT_TRUE(thermo1->input()["custom-field"]["second"].is<vector_fp>());
    auto spec1 = thermo1->species("NO");
    EXPECT_EQ(spec1->input["extra-field"], "blue");
    EXPECT_EQ(spec1->thermo->input()["bonus-field"], "green");
    EXPECT_EQ(spec1->transport->input["bogus-field"], "red");

    writer.skipUserDefined();
    AnyMap input2 = AnyMap::fromYamlString(writer.toYamlString());
    auto thermo2 = newPhase(input2["phases"].getMapWhere("name", "simple"),
                            input2);
    // user-defined fields should have been removed
    EXPECT_FALSE(thermo2->input().hasKey("custom-field"));
    auto spec2 = thermo2->species("NO");
    EXPECT_FALSE(spec2->input.hasKey("extra-field"));
    EXPECT_FALSE(spec2->thermo->input().hasKey("bonus-field"));
    EXPECT_FALSE(spec2->transport->input.hasKey("bogus-field"));
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
    auto original = newSolution("h2o2.yaml", "", "None");
    YamlWriter writer;
    writer.addPhase(original);
    writer.setPrecision(14);
    writer.toYamlFile("generated-h2o2.yaml");
    auto duplicate = newSolution("generated-h2o2.yaml", "", "None");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector_fp kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }
}

TEST(YamlWriter, reaction_units_from_Yaml)
{
    auto original = newSolution("h2o2.yaml", "", "None");
    YamlWriter writer;
    writer.addPhase(original);
    writer.setPrecision(14);
    writer.setUnits({
        {"activation-energy", "K"},
        {"quantity", "mol"},
        {"length", "cm"}
    });
    writer.toYamlFile("generated-h2o2-outunits.yaml");
    auto duplicate = newSolution("generated-h2o2-outunits.yaml", "", "None");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector_fp kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }
}

TEST(YamlWriter, reaction_units_from_Xml)
{
    auto original = newSolution("h2o2.xml", "", "None");
    YamlWriter writer;
    writer.addPhase(original);
    writer.setPrecision(14);
    writer.setUnits({
        {"activation-energy", "K"},
        {"quantity", "mol"},
        {"length", "cm"}
    });

    writer.toYamlFile("generated-h2o2-from-xml.yaml");
    auto duplicate = newSolution("generated-h2o2-from-xml.yaml", "", "None");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector_fp kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }
}

TEST(YamlWriter, chebyshev_units_from_Yaml)
{
    auto original = newSolution("pdep-test.yaml");
    YamlWriter writer;
    writer.addPhase(original);
    writer.setPrecision(14);
    writer.setUnits({
        {"activation-energy", "K"},
        {"quantity", "mol"},
        {"length", "cm"},
        {"pressure", "atm"}
    });
    writer.toYamlFile("generated-pdep-test.yaml");
    auto duplicate = newSolution("generated-pdep-test.yaml");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector_fp kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }
}

TEST(YamlWriter, multipleReactionSections)
{
    auto original1 = newSolution("h2o2.yaml", "", "None");
    auto original2 = newSolution("h2o2.yaml", "", "None");
    auto original3 = newSolution("h2o2.yaml", "", "None");
    // this phase will require its own "reactions" section
    auto R = original3->kinetics()->reaction(3);
    R->duplicate = true;
    R->rate()->releaseEvaluator(); // unlink from kinetics object
    original3->kinetics()->addReaction(R);
    original2->setName("ohmech2");
    original3->setName("ohmech3");
    YamlWriter writer;
    writer.addPhase(original1);
    writer.addPhase(original2);
    writer.addPhase(original3);
    writer.toYamlFile("generated-multi-rxn-secs.yaml");

    auto duplicate1 = newSolution("generated-multi-rxn-secs.yaml", "ohmech", "None");
    auto duplicate2 = newSolution("generated-multi-rxn-secs.yaml", "ohmech2", "None");
    auto duplicate3 = newSolution("generated-multi-rxn-secs.yaml", "ohmech3", "None");
    auto kin1 = duplicate1->kinetics();
    auto kin2 = duplicate2->kinetics();
    auto kin3 = duplicate3->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    ASSERT_EQ(kin2->nReactions() + 1, kin3->nReactions());
    ASSERT_EQ(kin2->reactionString(3),
              kin3->reactionString(kin3->nReactions() - 1));
}

TEST(YamlWriter, Interface)
{
    shared_ptr<ThermoPhase> gas1(newPhase("ptcombust.yaml", "gas"));
    shared_ptr<ThermoPhase> surf1(newPhase("ptcombust.yaml", "Pt_surf"));
    std::vector<ThermoPhase*> phases1{surf1.get(), gas1.get()};
    shared_ptr<Kinetics> kin1 = newKinetics(phases1, "ptcombust.yaml", "Pt_surf");

    double T = 900;
    double P = OneAtm;
    surf1->setState_TPX(T, P, "PT(S): 0.5, H(S): 0.1, CO(S): 0.4");
    gas1->setState_TPY(T, P, "H2: 0.5, CH4:0.48, OH:0.005, H:0.005");

    YamlWriter writer;
    writer.addPhase(gas1);
    writer.addPhase(surf1, kin1);
    writer.setUnits({
        {"length", "mm"},
        {"quantity", "molec"},
        {"activation-energy", "K"}
    });
    writer.toYamlFile("generated-ptcombust.yaml");

    shared_ptr<ThermoPhase> gas2(newPhase("generated-ptcombust.yaml", "gas"));
    shared_ptr<ThermoPhase> surf2(newPhase("generated-ptcombust.yaml", "Pt_surf"));
    std::vector<ThermoPhase*> phases2{surf2.get(), gas2.get()};
    shared_ptr<Kinetics> kin2 = newKinetics(phases2, "generated-ptcombust.yaml", "Pt_surf");

    auto iface1 = std::dynamic_pointer_cast<SurfPhase>(surf1);
    auto iface2 = std::dynamic_pointer_cast<SurfPhase>(surf2);

    EXPECT_NEAR(iface1->siteDensity(), iface2->siteDensity(),
                1e-13 * iface2->siteDensity());

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector_fp kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }

    vector_fp wdot1(kin1->nTotalSpecies());
    vector_fp wdot2(kin2->nTotalSpecies());
    kin1->getNetProductionRates(wdot1.data());
    kin2->getNetProductionRates(wdot2.data());
    for (size_t i = 0; i < kin1->nTotalSpecies(); i++) {
        EXPECT_NEAR(wdot1[i], wdot2[i], 1e-13 * fabs(wdot1[i])) << "for species i = " << i;
    }
}

TEST(YamlWriter, sofc)
{
    shared_ptr<ThermoPhase> gas1(newPhase("sofc.yaml", "gas"));
    shared_ptr<ThermoPhase> metal1(newPhase("sofc.yaml", "metal"));
    shared_ptr<ThermoPhase> ox_bulk1(newPhase("sofc.yaml", "oxide_bulk"));
    shared_ptr<ThermoPhase> metal_surf1(newPhase("sofc.yaml", "metal_surface"));
    shared_ptr<ThermoPhase> oxide_surf1(newPhase("sofc.yaml", "oxide_surface"));
    shared_ptr<ThermoPhase> tpb1(newPhase("sofc.yaml", "tpb"));

    std::vector<ThermoPhase*> tpb_phases1{tpb1.get(), metal_surf1.get(), oxide_surf1.get(), metal1.get()};
    std::vector<ThermoPhase*> ox_phases1{oxide_surf1.get(), ox_bulk1.get(), gas1.get()};

    shared_ptr<Kinetics> tpb_kin1 = newKinetics(tpb_phases1, "sofc.yaml", "tpb");
    shared_ptr<Kinetics> ox_kin1 = newKinetics(ox_phases1, "sofc.yaml", "oxide_surface");

    YamlWriter writer;
    writer.addPhase(tpb1, tpb_kin1);
    writer.addPhase(metal_surf1);
    writer.addPhase(oxide_surf1, ox_kin1);
    writer.addPhase(metal1);
    writer.addPhase(gas1);
    writer.addPhase(ox_bulk1);
    writer.setUnits({
        {"length", "cm"},
        {"pressure", "atm"},
        {"activation-energy", "eV"}
    });
    writer.toYamlFile("generated-sofc.yaml");

    shared_ptr<ThermoPhase> gas2(newPhase("generated-sofc.yaml", "gas"));
    shared_ptr<ThermoPhase> metal2(newPhase("generated-sofc.yaml", "metal"));
    shared_ptr<ThermoPhase> ox_bulk2(newPhase("generated-sofc.yaml", "oxide_bulk"));
    shared_ptr<ThermoPhase> metal_surf2(newPhase("generated-sofc.yaml", "metal_surface"));
    shared_ptr<ThermoPhase> oxide_surf2(newPhase("generated-sofc.yaml", "oxide_surface"));
    shared_ptr<ThermoPhase> tpb2(newPhase("generated-sofc.yaml", "tpb"));

    std::vector<ThermoPhase*> tpb_phases2{tpb2.get(), metal_surf2.get(), oxide_surf2.get(), metal2.get()};
    std::vector<ThermoPhase*> ox_phases2{oxide_surf2.get(), ox_bulk2.get(), gas2.get()};

    shared_ptr<Kinetics> tpb_kin2 = newKinetics(tpb_phases2, "generated-sofc.yaml", "tpb");
    shared_ptr<Kinetics> ox_kin2 = newKinetics(ox_phases2, "generated-sofc.yaml", "oxide_surface");

    ASSERT_EQ(tpb_kin1->nReactions(), tpb_kin2->nReactions());
    vector_fp kf1(tpb_kin1->nReactions()), kf2(tpb_kin1->nReactions());
    tpb_kin1->getFwdRateConstants(kf1.data());
    tpb_kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < tpb_kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for tpb reaction i = " << i;
    }

    vector_fp wdot1(tpb_kin1->nTotalSpecies());
    vector_fp wdot2(tpb_kin2->nTotalSpecies());
    tpb_kin1->getNetProductionRates(wdot1.data());
    tpb_kin2->getNetProductionRates(wdot2.data());
    for (size_t i = 0; i < tpb_kin1->nTotalSpecies(); i++) {
        EXPECT_NEAR(wdot1[i], wdot2[i], 1e-13 * fabs(wdot1[i])) << "for species i = " << i;
    }

    ASSERT_EQ(ox_kin1->nReactions(), ox_kin2->nReactions());
    kf1.resize(ox_kin1->nReactions());
    kf2.resize(ox_kin1->nReactions());
    ox_kin1->getFwdRateConstants(kf1.data());
    ox_kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < ox_kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for ox reaction i = " << i;
    }

    wdot1.resize(ox_kin1->nTotalSpecies());
    wdot2.resize(ox_kin2->nTotalSpecies());
    ox_kin1->getNetProductionRates(wdot1.data());
    ox_kin2->getNetProductionRates(wdot2.data());
    for (size_t i = 0; i < ox_kin1->nTotalSpecies(); i++) {
        EXPECT_NEAR(wdot1[i], wdot2[i], 1e-13 * fabs(wdot1[i])) << "for ox species i = " << i;
    }
}
