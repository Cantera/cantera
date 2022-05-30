#include "gtest/gtest.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/Solution.h"

using namespace std;

namespace Cantera
{

vector<AnyMap> getStates(const string& name) {
    static AnyMap cases = AnyMap::fromYamlFile("consistency-cases.yaml");
    return cases[name]["states"].asVector<AnyMap>();
}

AnyMap getSetup(const string& name) {
    static AnyMap cases = AnyMap::fromYamlFile("consistency-cases.yaml");
    return cases[name]["input"].as<AnyMap>();
}

class TestConsistency : public testing::TestWithParam<std::tuple<AnyMap, AnyMap>>
{
public:
    TestConsistency() {
        auto param = GetParam();
        AnyMap inp = get<0>(param);
        AnyMap state = get<1>(param);
        pair<string, string> key = {inp["file"].asString(), inp.getString("phase", "")};
        if (cache.count(key) == 0) {
            cache[key].reset(newPhase(key.first, key.second));
        }
        phase = cache[key];
        phase->setState(state);
    }
    static map<pair<string, string>, shared_ptr<ThermoPhase>> cache;
    shared_ptr<ThermoPhase> phase;
};

map<pair<string, string>, shared_ptr<ThermoPhase>> TestConsistency::cache = {};

TEST_P(TestConsistency, h_eq_u_plus_Pv) {
    double h = phase->enthalpy_mole();
    double u = phase->intEnergy_mole();
    double v = phase->molarVolume();
    double p = phase->pressure();
    EXPECT_NEAR(h, u + p*v, 1e-10*p*v);
}

TEST_P(TestConsistency, g_eq_h_minus_Ts) {
    double g = phase->gibbs_mole();
    double h = phase->enthalpy_mole();
    double T = phase->temperature();
    double s = phase->entropy_mole();
    EXPECT_NEAR(g, h - T*s, 1e-10*T*s);
}

INSTANTIATE_TEST_SUITE_P(IdealGas, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-gas-h2o2")),
        testing::ValuesIn(getStates("ideal-gas-h2o2")))
);

INSTANTIATE_TEST_SUITE_P(RedlichKwong, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("redlich-kwong")),
        testing::ValuesIn(getStates("redlich-kwong")))
);

}
