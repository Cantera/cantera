#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace Cantera;

shared_ptr<ThermoPhase> newThermo(const std::string& fileName,
                                  const std::string& phaseName)
{
    return shared_ptr<ThermoPhase>(newPhase(fileName, phaseName));
}

TEST(ThermoToYaml, simpleIdealGas)
{
    auto thermo = newThermo("ideal-gas.yaml", "simple");
    AnyMap data;
    thermo->setState_TP(1010, 2e5);
    double rho = thermo->density();
    thermo->getParameters(data);

    ASSERT_EQ(data["thermo"], "ideal-gas");
    ASSERT_EQ(data["state"]["T"], 1010);
    ASSERT_EQ(data["state"]["density"], rho);
}
