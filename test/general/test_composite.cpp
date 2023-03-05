#include "gtest/gtest.h"
#include "cantera/core.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/transport/MixTransport.h"

using namespace Cantera;

TEST(Solution, valid_cast)
{
    auto soln = newSolution("h2o2.yaml", "ohmech");
    ASSERT_TRUE(soln->thermo<IdealGasPhase>());
    auto hRT_ref = soln->thermo<IdealGasPhase>()->enthalpy_RT_ref();
    ASSERT_EQ(hRT_ref.size(), 10);
    ASSERT_TRUE(soln->kinetics<GasKinetics>());
    soln->kinetics<GasKinetics>()->updateROP();
    ASSERT_EQ(soln->kinetics<GasKinetics>()->nReactions(), 29);
    ASSERT_TRUE(soln->transport<MixTransport>());
}

TEST(Solution, invalid_cast)
{
    auto soln = newSolution("h2o2.yaml", "ohmech", "none");
    ASSERT_THROW(soln->thermo<SurfPhase>(), CanteraError);
    ASSERT_THROW(soln->kinetics<InterfaceKinetics>(), CanteraError);
    ASSERT_THROW(soln->transport<MixTransport>(), CanteraError);
}

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
    ASSERT_EQ(surf->kinetics()->nReactions(), 24u);
}
