// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "cantera/core.h"
#include "cantera/base/YamlWriter.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/transport/TransportData.h"
#include "cantera/base/Storage.h"
#include <fstream>

using namespace Cantera;
using namespace YAML;

TEST(YamlWriter, formatDouble)
{
    AnyMap m;
    int count;
    double delta;

    // check least significant digit
    // default precision is 15 (see AnyMap.cpp::getPrecision)
    delta = 1.e-15;
    count = 0;
    for (int i = -9; i < 10; i += 2) {
        m["a"] = 4. + i * delta;
        if (i < -5 || i > 4) {
            // round away from 4.
            EXPECT_NE(m.toYamlString(), "a: 4.0\n");
        } else {
            // round towards 4.
            EXPECT_EQ(m.toYamlString(), "a: 4.0\n");
            count++;
        }
    }
    EXPECT_EQ(count, 5); // only checking even multiples of delta

    // check edge cases
    m["a"] = -1061.793215682400;
    EXPECT_EQ(m.toYamlString(), "a: -1061.7932156824\n");

    m["a"] = 7.820059054328200e-111;
    EXPECT_EQ(m.toYamlString(), "a: 7.8200590543282e-111\n");
}

TEST(YamlWriter, formatDoubleExp)
{
    AnyMap m;
    int count;
    double delta;

    // check least significant digit
    // default precision is 15 (see AnyMap.cpp::getPrecision)
    delta = 1.e-5;
    count = 0;
    for (int i = -9; i < 10; i += 2) {
        m["a"] = 4.e10 + i * delta;
        if (i < -5 || i > 4) {
            // round away from 4.
            EXPECT_NE(m.toYamlString(), "a: 4.0e+10\n");
        } else {
            // round towards 4.
            EXPECT_EQ(m.toYamlString(), "a: 4.0e+10\n");
            count++;
        }
    }
    EXPECT_EQ(count, 5); // only checking even multiples of delta

    // check edge cases
    m["a"] = 1.629771953878800e+13;
    EXPECT_EQ(m.toYamlString(), "a: 1.6297719538788e+13\n");
}

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
    auto thermo1 = newThermo(input1["phases"].getMapWhere("name", "simple"),
                                  input1);

    // user-defined fields should be in place
    EXPECT_TRUE(thermo1->input()["custom-field"]["second"].is<vector<double>>());
    auto spec1 = thermo1->species("NO");
    EXPECT_EQ(spec1->input["extra-field"], "blue");
    EXPECT_EQ(spec1->thermo->input()["bonus-field"], "green");
    EXPECT_EQ(spec1->transport->input["bogus-field"], "red");

    writer.skipUserDefined();
    AnyMap input2 = AnyMap::fromYamlString(writer.toYamlString());
    auto thermo2 = newThermo(input2["phases"].getMapWhere("name", "simple"),
                                  input2);
    // user-defined fields should have been removed
    EXPECT_FALSE(thermo2->input().hasKey("custom-field"));
    auto spec2 = thermo2->species("NO");
    EXPECT_FALSE(spec2->input.hasKey("extra-field"));
    EXPECT_FALSE(spec2->thermo->input().hasKey("bonus-field"));
    EXPECT_FALSE(spec2->transport->input.hasKey("bogus-field"));
}

TEST(YamlWriter, literalStrings)
{
    auto original = newSolution("ideal-gas.yaml", "simple");
    YamlWriter writer;
    writer.addPhase(original);
    writer.toYamlFile("generated-literal.yaml");
    auto duplicate = newSolution("generated-literal.yaml");
    auto thermo1 = original->thermo();
    auto thermo2 = duplicate->thermo();

    EXPECT_EQ(thermo1->input()["literal-string"].asString(),
              thermo2->input()["literal-string"].asString());

    auto spec1 = thermo1->species("O2");
    auto spec2 = thermo2->species("O2");

    EXPECT_EQ(spec1->input["another-literal-string"].asString(),
              spec2->input["another-literal-string"].asString());
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
    auto original = newSolution("h2o2.yaml", "", "none");
    YamlWriter writer;
    writer.addPhase(original);
    writer.setPrecision(14);
    writer.toYamlFile("generated-h2o2.yaml");
    auto duplicate = newSolution("generated-h2o2.yaml", "", "none");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector<double> kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }
}

TEST(YamlWriter, reaction_units_from_Yaml)
{
    auto original = newSolution("h2o2.yaml", "", "none");
    YamlWriter writer;
    writer.addPhase(original);
    writer.setPrecision(14);
    auto units = UnitSystem();
    map<string, string> defaults{
        {"activation-energy", "K"},
        {"quantity", "mol"},
        {"length", "cm"}
    };
    units.setDefaults(defaults);
    writer.setUnitSystem(units);
    writer.toYamlFile("generated-h2o2-outunits.yaml");
    auto duplicate = newSolution("generated-h2o2-outunits.yaml", "", "none");

    auto kin1 = original->kinetics();
    auto kin2 = duplicate->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector<double> kf1(kin1->nReactions()), kf2(kin1->nReactions());
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
    vector<double> kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }
}

TEST(YamlWriter, multipleReactionSections)
{
    auto original1 = newSolution("h2o2.yaml", "", "none");
    auto original2 = newSolution("h2o2.yaml", "", "none");
    auto original3 = newSolution("h2o2.yaml", "", "none");
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

    auto duplicate1 = newSolution("generated-multi-rxn-secs.yaml", "ohmech", "none");
    auto duplicate2 = newSolution("generated-multi-rxn-secs.yaml", "ohmech2", "none");
    auto duplicate3 = newSolution("generated-multi-rxn-secs.yaml", "ohmech3", "none");
    auto kin1 = duplicate1->kinetics();
    auto kin2 = duplicate2->kinetics();
    auto kin3 = duplicate3->kinetics();

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    ASSERT_EQ(kin2->nReactions() + 1, kin3->nReactions());
    ASSERT_EQ(kin2->reaction(3)->equation(),
              kin3->reaction(kin3->nReactions() - 1)->equation());
}

TEST(YamlWriter, Interface)
{
    shared_ptr<ThermoPhase> gas1(newThermo("ptcombust.yaml", "gas"));
    shared_ptr<ThermoPhase> surf1(newThermo("ptcombust.yaml", "Pt_surf"));
    vector<shared_ptr<ThermoPhase>> phases1{surf1, gas1};
    shared_ptr<Kinetics> kin1 = newKinetics(phases1, "ptcombust.yaml");

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

    shared_ptr<ThermoPhase> gas2(newThermo("generated-ptcombust.yaml", "gas"));
    shared_ptr<ThermoPhase> surf2(newThermo("generated-ptcombust.yaml", "Pt_surf"));
    vector<shared_ptr<ThermoPhase>> phases2{surf2, gas2};
    shared_ptr<Kinetics> kin2 = newKinetics(phases2, "generated-ptcombust.yaml");

    auto iface1 = std::dynamic_pointer_cast<SurfPhase>(surf1);
    auto iface2 = std::dynamic_pointer_cast<SurfPhase>(surf2);

    EXPECT_NEAR(iface1->siteDensity(), iface2->siteDensity(),
                1e-13 * iface2->siteDensity());

    ASSERT_EQ(kin1->nReactions(), kin2->nReactions());
    vector<double> kf1(kin1->nReactions()), kf2(kin1->nReactions());
    kin1->getFwdRateConstants(kf1.data());
    kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for reaction i = " << i;
    }

    vector<double> wdot1(kin1->nTotalSpecies());
    vector<double> wdot2(kin2->nTotalSpecies());
    kin1->getNetProductionRates(wdot1.data());
    kin2->getNetProductionRates(wdot2.data());
    for (size_t i = 0; i < kin1->nTotalSpecies(); i++) {
        EXPECT_NEAR(wdot1[i], wdot2[i], 1e-13 * fabs(wdot1[i])) << "for species i = " << i;
    }
}

TEST(YamlWriter, sofc)
{
    auto tpb1 = newSolution("sofc.yaml", "tpb");
    auto tpb_kin1 = tpb1->kinetics();
    auto ox_kin1 = tpb1->adjacent("oxide_surface")->kinetics();

    YamlWriter writer;
    writer.addPhase(tpb1);
    writer.setUnits({
        {"length", "cm"},
        {"pressure", "atm"},
        {"activation-energy", "eV"}
    });
    writer.skipUserDefined();
    writer.toYamlFile("generated-sofc.yaml");

    auto tpb2 = newSolution("generated-sofc.yaml", "tpb");
    auto tpb_kin2 = tpb2->kinetics();
    auto ox_kin2 = tpb2->adjacent("oxide_surface")->kinetics();

    ASSERT_EQ(tpb_kin1->nReactions(), tpb_kin2->nReactions());
    vector<double> kf1(tpb_kin1->nReactions()), kf2(tpb_kin1->nReactions());
    tpb_kin1->getFwdRateConstants(kf1.data());
    tpb_kin2->getFwdRateConstants(kf2.data());
    for (size_t i = 0; i < tpb_kin1->nReactions(); i++) {
        EXPECT_NEAR(kf1[i], kf2[i], 1e-13 * kf1[i]) << "for tpb reaction i = " << i;
    }

    vector<double> wdot1(tpb_kin1->nTotalSpecies());
    vector<double> wdot2(tpb_kin2->nTotalSpecies());
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

TEST(YamlWriter, customHeader)
{
    auto original = newSolution("h2o2.yaml", "", "none");
    AnyMap header;
    header["description"] = "Copy of H2O2 mechanism";
    header["spam"] = "eggs";

    YamlWriter writer;
    writer.setHeader(header);
    writer.addPhase(original);
    writer.toYamlFile("generated-custom-header.yaml");

    auto soln = newSolution("generated-custom-header.yaml", "", "none");

    ASSERT_EQ(soln->header()["cantera-version"].asString(), CANTERA_VERSION);
    ASSERT_EQ(soln->header()["generator"].asString(), "YamlWriter");
    ASSERT_EQ(soln->header()["description"].asString(),
              "Copy of H2O2 mechanism");
    ASSERT_EQ(soln->header()["spam"].asString(), "eggs");
}

#if CT_USE_HDF5

TEST(Storage, groups)
{
    // testing Storage class outside of SolutionArray
    const string fname = "groups.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }

    // open in write mode
    auto file = unique_ptr<Storage>(new Storage(fname, true));

    // create group implicitly if permissive flag is set
    EXPECT_FALSE(file->hasGroup("one"));
    EXPECT_THROW(file->checkGroup("one"), CanteraError);
    EXPECT_FALSE(file->checkGroup("one", true));
    EXPECT_TRUE(file->hasGroup("one"));

    // nested groups
    EXPECT_FALSE(file->hasGroup("one/two"));
    EXPECT_FALSE(file->checkGroup("one/two", true));
    EXPECT_TRUE(file->hasGroup("one/two"));
    EXPECT_FALSE(file->hasGroup("two"));

    // open in read mode
    file = unique_ptr<Storage>(new Storage(fname, false));

    // can check group, but not create
    EXPECT_FALSE(file->hasGroup("three"));
    EXPECT_THROW(file->checkGroup("three"), CanteraError);
    EXPECT_FALSE(file->checkGroup("three", true));
    EXPECT_FALSE(file->hasGroup("three"));
}

TEST(Storage, writeAttributes)
{
    // testing Storage class outside of SolutionArray
    const string fname = "writeAttributes.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }

    // open in write mode
    auto file = unique_ptr<Storage>(new Storage(fname, true));

    // try to write without creating group first
    AnyMap attr;
    attr["description"] = "writing attributes";
    // cannot write to root location
    EXPECT_THROW(file->writeAttributes("", attr), CanteraError);
    // cannot write to group that is not created
    EXPECT_THROW(file->writeAttributes("one", attr), CanteraError);

    file->checkGroup("one", true); // implicitly creates group
    EXPECT_FALSE(file->hasAttribute("one", "description"));
    file->writeAttributes("one", attr);
    EXPECT_TRUE(file->hasAttribute("one", "description"));
    // attributes cannot be overwritten
    EXPECT_THROW(file->writeAttributes("one", attr), NotImplementedError);

    attr.clear(); // reset AnyMap
    attr["spam"] = "eggs"; // string
    attr["foo"]["bar]"] = 1.; // float; also nested
    file->writeAttributes("one", attr); // succeds; new attribute

    file->checkGroup("two/three", true); // multi-level group
    file->writeAttributes("two/three", attr);

    attr["hello-world"] = 5; // integer
    attr["animals"] = vector<string>({"dog", "cat"});
    attr["floats"] = vector<double>({1.1, 2.2, 3.3});
    attr["integers"] = vector<long int>({1, 2, 3});
    attr["booleans"] = vector<bool>({true, false, true});
    file->checkGroup("four", true);
    file->writeAttributes("four", attr);

    // unsupported types
    attr.clear();
    attr["any"] = AnyValue();
    EXPECT_THROW(file->writeAttributes("four", attr), NotImplementedError);
    attr.clear();
    attr["invalid"] = vector<int>({1, 2, 3}); // not supported (needs long int)
    EXPECT_THROW(file->writeAttributes("four", attr), NotImplementedError);
}

TEST(Storage, readAttributes)
{
    // testing Storage class outside of SolutionArray
    const string fname = "readAttributes.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }
    // open in write mode
    auto file = unique_ptr<Storage>(new Storage(fname, true));

    AnyMap attr;
    attr["description"] = "reading attributes";
    attr["spam"] = "eggs"; // string
    attr["foo"]["bar"] = 1.; // float
    attr["hello-world"] = 5; // integer
    attr["animals"] = vector<string>({"dog", "cat"});
    attr["floats"] = vector<double>({1.1, 2.2, 3.3});
    attr["integers"] = vector<long int>({1, 2, 3});
    attr["booleans"] = vector<bool>({true, false, true});

    file->checkGroup("test", true); // implicitly creates group
    file->writeAttributes("test", attr);

    // open in read mode
    file = unique_ptr<Storage>(new Storage(fname, false));
    auto data = file->readAttributes("test", false); // only one level
    EXPECT_EQ(data["spam"].asString(), attr["spam"].asString());
    EXPECT_FALSE(data.hasKey("foo"));

    data = file->readAttributes("test", true); // recursive
    EXPECT_EQ(data["spam"].asString(), attr["spam"].asString());
    EXPECT_TRUE(data.hasKey("foo"));
    EXPECT_EQ(data["foo"]["bar"].asDouble(), attr["foo"]["bar"].asDouble());

    // cannot write to file in read mode
    AnyMap any;
    any["hello"] = "world";
    EXPECT_TRUE(file->hasGroup("test"));
    EXPECT_THROW(file->writeAttributes("test", any), CanteraError);

    data = file->readAttributes("test", false); // recursive
    EXPECT_FALSE(data.hasKey("hello"));
}

TEST(Storage, writeData)
{
    // testing Storage class outside of SolutionArray
    const string fname = "writeData.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }
    auto file = unique_ptr<Storage>(new Storage(fname, true));
    file->checkGroup("test", true); // implicitly creates group

    AnyValue any;

    // scalars don't work
    any = 3.1415;
    EXPECT_THROW(file->writeData("test", "double", any), CanteraError);
    any = 42;
    EXPECT_THROW(file->writeData("test", "integer", any), CanteraError);
    any = "hello world!";
    EXPECT_THROW(file->writeData("test", "string", any), CanteraError);

    // vectors
    any = vector<double>({1.1, 2.2, 3.3});
    file->writeData("test", "double-vector", any);
    any = vector<long int>({1, 2, 3, 4});
    file->writeData("test", "integer-vector", any);
    any = vector<string>({"dog", "cat"});
    file->writeData("test", "string-vector", any);

    // writing to root doesn't work
    EXPECT_THROW(file->writeData("", "string-vector", any), CanteraError);
    // overwriting of existing data doesn't work
    EXPECT_THROW(file->writeData("test", "string-vector", any), CanteraError);

    // not all types are implemented
    any = vector<bool>({true, false, true});
    EXPECT_THROW(file->writeData("test", "invalid0", any), NotImplementedError);
    any = vector<int>({1, 2, 3});
    EXPECT_THROW(file->writeData("test", "invalid1", any), CanteraError);

    // matrices
    any = vector<vector<double>>({{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}});
    EXPECT_EQ(any.matrixShape().first, 2u);
    EXPECT_EQ(any.matrixShape().second, 3u);
    file->writeData("test", "double-matrix", any);
    any = vector<vector<long int>>({{1, 2}, {3, 4}, {5, 6}});
    file->writeData("test", "integer-matrix", any);
    any = vector<vector<string>>({{"foo", "bar"}, {"spam", "eggs"}});
    file->writeData("test", "string-matrix", any);

    // not all types are implemented
    any = vector<vector<int>>({{1, 2}, {3, 4}});
    EXPECT_THROW(file->writeData("test", "invalid2", any), CanteraError);

    // invalid vectors of vectors (irregular shape)
    any = vector<vector<double>>({{1.1, 2.2}, {3.3}});
    EXPECT_EQ(any.matrixShape().first, 2u);
    EXPECT_EQ(any.matrixShape().second, npos);
    EXPECT_THROW(file->writeData("test", "invalid3", any), CanteraError);
    any = vector<vector<long int>>({{1, 2}, {3, 4, 5}, {6}});
    EXPECT_THROW(file->writeData("test", "invalid4", any), CanteraError);
    any = vector<vector<string>>({{"foo", "bar"}, {"a", "b", "c"}});
    EXPECT_THROW(file->writeData("test", "invalid5", any), CanteraError);
}

TEST(Storage, writeCompressed)
{
    // testing Storage class outside of SolutionArray
    const string fname = "writeCompressed.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }
    auto file = unique_ptr<Storage>(new Storage(fname, true));
    file->checkGroup("test", true); // implicitly creates group

    AnyValue any;

    EXPECT_THROW(file->setCompressionLevel(-1), CanteraError);
    EXPECT_THROW(file->setCompressionLevel(10), CanteraError);
    file->setCompressionLevel(9);
    file->setCompressionLevel(0);
    file->setCompressionLevel(5);

    // matrices
    any = vector<vector<double>>({{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}});
    file->writeData("test", "double-matrix", any);
    any = vector<vector<long int>>({{1, 2}, {3, 4}, {5, 6}});
    file->writeData("test", "integer-matrix", any);
    any = vector<vector<string>>({{"foo", "bar"}, {"spam", "eggs"}});
    file->writeData("test", "string-matrix", any);
}

TEST(Storage, readData)
{
    // testing Storage class outside of SolutionArray
    const string fname = "writeData.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }
    auto file = unique_ptr<Storage>(new Storage(fname, true));
    file->checkGroup("test", true); // implicitly creates group

    AnyValue any;

    // vectors
    any = vector<double>({1.1, 2.2, 3.3});
    file->writeData("test", "double-vector", any);
    any = vector<long int>({1, 2, 3, 4});
    file->writeData("test", "integer-vector", any);
    any = vector<string>({"dog", "cat"});
    file->writeData("test", "string-vector", any);

    // matrices
    any = vector<vector<double>>({{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}});
    file->writeData("test", "double-matrix", any);
    any = vector<vector<long int>>({{1, 2}, {3, 4}, {5, 6}});
    file->writeData("test", "integer-matrix", any);
    any = vector<vector<string>>({{"foo", "bar"}, {"spam", "eggs"}});
    file->writeData("test", "string-matrix", any);

    file = unique_ptr<Storage>(new Storage(fname, false));

    // cannot write to file in read mode
    EXPECT_THROW(file->writeData("test", "string-matrix", any), CanteraError);

    auto data = file->readData("test", "double-vector", 3);
    ASSERT_TRUE(data.isVector<double>());
    ASSERT_EQ(data.asVector<double>().size(), 3u);
    // need to know length
    EXPECT_THROW(file->readData("test", "double-vector", 4), CanteraError);

    data = file->readData("test", "integer-vector", 4, 0);
    ASSERT_TRUE(data.isVector<long int>());
    ASSERT_EQ(data.asVector<long int>().size(), 4u);
    // invalid number of columns
    EXPECT_THROW(file->readData("test", "integer-vector", 4, 2), CanteraError);

    data = file->readData("test", "string-vector", 2);
    ASSERT_TRUE(data.isVector<string>());
    ASSERT_EQ(data.asVector<string>().size(), 2u);

    data = file->readData("test", "double-matrix", 2);
    ASSERT_TRUE(data.isMatrix<double>());
    // incorrect number of rows
    EXPECT_THROW(file->readData("test", "double-matrix", 3), CanteraError);

    data = file->readData("test", "integer-matrix", 3, 2);
    ASSERT_TRUE(data.isMatrix<long int>());
    // invalid number of columns
    EXPECT_THROW(file->readData("test", "integer-matrix", 3, 4), CanteraError);

    data = file->readData("test", "string-matrix", 2, 2);
    ASSERT_TRUE(data.isMatrix<string>());
}

TEST(Storage, readCompressed)
{
    // testing Storage class outside of SolutionArray
    const string fname = "writeCompressed.h5";
    if (std::ifstream(fname).good()) {
        std::remove(fname.c_str());
    }
    auto file = unique_ptr<Storage>(new Storage(fname, true));
    file->checkGroup("test", true); // implicitly creates group
    file->setCompressionLevel(5);

    // matrices
    AnyValue any;
    any = vector<vector<double>>({{1.1, 2.2, 3.3}, {4.4, 5.5, 6.6}});
    file->writeData("test", "double-matrix", any);
    any = vector<vector<long int>>({{1, 2}, {3, 4}, {5, 6}});
    file->writeData("test", "integer-matrix", any);
    any = vector<vector<string>>({{"foo", "bar"}, {"spam", "eggs"}});
    file->writeData("test", "string-matrix", any);

    file = unique_ptr<Storage>(new Storage(fname, false));

    auto data = file->readData("test", "double-matrix", 2, 3);
    ASSERT_TRUE(data.isMatrix<double>());

    data = file->readData("test", "integer-matrix", 3, 2);
    ASSERT_TRUE(data.isMatrix<long int>());

    data = file->readData("test", "string-matrix", 2, 2);
    ASSERT_TRUE(data.isMatrix<string>());
}

#else

TEST(Storage, noSupport)
{
    EXPECT_THROW(unique_ptr<Storage>(new Storage("something.h5", true)), CanteraError);
}

#endif
