#include "gtest/gtest.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/Solution.h"
#include "cantera/base/utilities.h"
#include <regex>

using namespace std;

namespace Cantera
{

vector<AnyMap> getStates(const string& name) {
    static AnyMap cases = AnyMap::fromYamlFile("consistency-cases.yaml");
    return cases[name]["states"].asVector<AnyMap>();
}

AnyMap getSetup(const string& name) {
    static AnyMap cases = AnyMap::fromYamlFile("consistency-cases.yaml");
    return cases[name]["setup"].as<AnyMap>();
}

// For more informative output about failing test cases
std::ostream& operator<<(std::ostream& s, const AnyMap& m)
{
    if (m.hasKey("file")) {
        s << fmt::format("file: {}, phase: {}",
                         m["file"].asString(), m.getString("phase", "<default>"));
    } else {
        s << "\n" << m.toYamlString();
    }
    return s;
}

class TestConsistency : public testing::TestWithParam<std::tuple<AnyMap, AnyMap>>
{
public:
    TestConsistency() {
        auto param = GetParam();
        setup = get<0>(param);
        AnyMap state = get<1>(param);
        pair<string, string> key = {setup["file"].asString(), setup.getString("phase", "")};
        if (cache.count(key) == 0) {
            cache[key].reset(newPhase(key.first, key.second));
        }
        atol = setup.getDouble("atol", 1e-5);
        atol_v = setup.getDouble("atol_v", 1e-11);

        phase = cache[key];
        phase->setState(state);
        nsp = phase->nSpecies();
        p = phase->pressure();
        T = phase->temperature();
    }

    void SetUp() {
        if (setup.hasKey("known-failures")) {
            auto name = testing::UnitTest::GetInstance()->current_test_info()->name();
            for (auto& item : setup["known-failures"].asMap<string>()) {
                if (regex_search(name, regex(item.first))) {
                    GTEST_SKIP() << item.second;
                    break;
                }
            }
        }
    }

    static map<pair<string, string>, shared_ptr<ThermoPhase>> cache;

    AnyMap setup;
    shared_ptr<ThermoPhase> phase;
    size_t nsp;
    double T, p;
    double atol, atol_v;
};

map<pair<string, string>, shared_ptr<ThermoPhase>> TestConsistency::cache = {};

TEST_P(TestConsistency, h_eq_u_plus_Pv) {
    double h = phase->enthalpy_mole();
    double u = phase->intEnergy_mole();
    double v = phase->molarVolume();
    EXPECT_NEAR(h, u + p * v, atol);
}

TEST_P(TestConsistency, g_eq_h_minus_Ts) {
    double g = phase->gibbs_mole();
    double h = phase->enthalpy_mole();
    double s = phase->entropy_mole();
    EXPECT_NEAR(g, h - T * s, atol);
}

TEST_P(TestConsistency, hk_eq_uk_plus_P_times_vk)
{
    vector_fp hk(nsp), uk(nsp), vk(nsp);
    phase->getPartialMolarEnthalpies(hk.data());
    phase->getPartialMolarIntEnergies(uk.data());
    phase->getPartialMolarVolumes(vk.data());
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(hk[k], uk[k] + p * vk[k], atol) << "k = " << k;
    }
}

TEST_P(TestConsistency, gk_eq_hk_minus_T_times_sk)
{
    vector_fp gk(nsp), hk(nsp), sk(nsp);
    phase->getChemPotentials(gk.data());
    phase->getPartialMolarEnthalpies(hk.data());
    phase->getPartialMolarEntropies(sk.data());
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(gk[k], hk[k] - T * sk[k], atol) << "k = " << k;
    }
}

TEST_P(TestConsistency, h_eq_sum_hk_Xk)
{
    vector_fp hk(nsp);
    phase->getPartialMolarEnthalpies(hk.data());
    EXPECT_NEAR(phase->enthalpy_mole(), phase->mean_X(hk), atol);
}

TEST_P(TestConsistency, u_eq_sum_uk_Xk)
{
    vector_fp uk(nsp);
    phase->getPartialMolarIntEnergies(uk.data());
    EXPECT_NEAR(phase->intEnergy_mole(), phase->mean_X(uk), atol);
}

TEST_P(TestConsistency, g_eq_sum_gk_Xk)
{
    vector_fp gk(nsp);
    phase->getChemPotentials(gk.data());
    EXPECT_NEAR(phase->gibbs_mole(), phase->mean_X(gk), atol);
}

TEST_P(TestConsistency, s_eq_sum_sk_Xk)
{
    vector_fp sk(nsp);
    phase->getPartialMolarEntropies(sk.data());
    EXPECT_NEAR(phase->entropy_mole(), phase->mean_X(sk), atol);
}

TEST_P(TestConsistency, v_eq_sum_vk_Xk)
{
    vector_fp vk(nsp);
    phase->getPartialMolarVolumes(vk.data());
    EXPECT_NEAR(phase->molarVolume(), phase->mean_X(vk), atol_v);
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

INSTANTIATE_TEST_SUITE_P(PengRobinson, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("peng-robinson")),
        testing::ValuesIn(getStates("peng-robinson")))
);

INSTANTIATE_TEST_SUITE_P(IdealMolalSolution, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-molal-solution")),
        testing::ValuesIn(getStates("ideal-molal-solution")))
);

INSTANTIATE_TEST_SUITE_P(IdealSolidSolnPhase1, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-condensed-1")),
        testing::ValuesIn(getStates("ideal-condensed-1")))
);

INSTANTIATE_TEST_SUITE_P(IdealSolidSolnPhase2, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-condensed-2")),
        testing::ValuesIn(getStates("ideal-condensed-2")))
);

INSTANTIATE_TEST_SUITE_P(BinarySolutionTabulated, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("binary-solution-tabulated")),
        testing::ValuesIn(getStates("binary-solution-tabulated")))
);

INSTANTIATE_TEST_SUITE_P(ElectronCloud, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("electron-cloud")),
        testing::ValuesIn(getStates("electron-cloud")))
);

INSTANTIATE_TEST_SUITE_P(NitrogenPureFluid, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("nitrogen-purefluid")),
        testing::ValuesIn(getStates("nitrogen-purefluid")))
);

INSTANTIATE_TEST_SUITE_P(PlasmaPhase, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("plasma")),
        testing::ValuesIn(getStates("plasma")))
);

}
