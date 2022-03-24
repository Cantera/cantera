#include "gtest/gtest.h"
#include "cantera/base/Interface.h"

using namespace Cantera;

TEST(Interface, incompatible_phase)
{
    ASSERT_THROW(newInterface("h2o2.yaml", "ohmech"), CanteraError);
    AnyMap h2o2 = AnyMap::fromYamlFile("h2o2.yaml");
    auto& ohmech = h2o2["phases"].getMapWhere("name", "ohmech");
    ASSERT_THROW(newInterface(ohmech, h2o2), CanteraError);
}

TEST(Interface, from_string_and_Solution)
{
    auto gas = newSolution("ptcombust.yaml", "gas");
    auto surf = newInterface("ptcombust.yaml", "Pt_surf", {gas});
    ASSERT_EQ(gas.get(), surf->adjacent(0).get());
}

TEST(Interface, from_AnyMap_and_Solution)
{
    auto gas = newSolution("ptcombust.yaml", "gas");
    AnyMap root = AnyMap::fromYamlFile("ptcombust.yaml");
    auto& phase = root["phases"].getMapWhere("name", "Pt_surf");
    auto surf = newInterface(phase, root, {gas});
    ASSERT_EQ(gas.get(), surf->adjacent(0).get());
    ASSERT_EQ(surf->kinetics()->nReactions(), 24);
}
