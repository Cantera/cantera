// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/thermo/PlasmaPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/base/Solution.h"
#include "cantera/base/utilities.h"
#include <boost/algorithm/string.hpp>
#include <regex>

// This is a set of tests that check all ThermoPhase models implemented in Cantera
// for thermodynamic identities such as g = h - T*s and c_p = dh/dT at constant P.

// Below the definition of the TestConsistency test fixture class, the individual
// consistency tests are implemented in the functions declared using the TEST_P macro.
// Each of these tests starts with the 'phase' object created and set to the
// thermodynamic state that should be used for the test.

// Following the individual test definitions, a "test suite" is created for each phase
// model using the INSTANTIATE_TEST_SUITE_P macro. The string passed to the getSetup()
// and getStates() functions in the test suite instantiation corresponds to a top-level
// section of the file test/data/consistency-cases.yaml. Each of these sections defines
// a phase definition and a sequence of thermodynamic states to be used for each of the
// individual test functions. Generally, there is one test suite for each ThermoPhase
// model. However, there are a few phase models that have multiple test suites to allow
// testing the effects of parameters that can only be varied within the phase
// definition.

using namespace std;
using namespace Cantera;

// Helper functions to reduce code duplication in test suite instantiations
vector<AnyMap> getStates(const string& name) {
    static AnyMap cases = AnyMap::fromYamlFile("consistency-cases.yaml");
    return cases[name]["states"].asVector<AnyMap>();
}

AnyMap getSetup(const string& name) {
    static AnyMap cases = AnyMap::fromYamlFile("consistency-cases.yaml");
    return cases[name]["setup"].as<AnyMap>();
}

namespace Cantera {

// For more informative output about failing test cases
std::ostream& operator<<(std::ostream& s, const AnyMap& m)
{
    if (m.hasKey("phase")) {
        s << fmt::format("file: {}, phase: {}",
                         m["file"].asString(), m["phase"].asString());
    } else if (m.hasKey("file")) {
        s << fmt::format("file: {}", m["file"].asString());
    } else {
        AnyMap out = m;
        out.setFlowStyle();
        s << "state: "
          << boost::algorithm::replace_all_copy(out.toYamlString(), "\n", " ");
    }
    return s;
}

} // end namespace Cantera

class TestConsistency : public testing::TestWithParam<std::tuple<AnyMap, AnyMap>>
{
public:
    TestConsistency() {
        auto param = GetParam();
        setup = get<0>(param);
        AnyMap state = get<1>(param);

        if (setup.getBool("ignore-deprecations", false)) {
            suppress_deprecation_warnings();
        }

        // For efficiency, cache the instantiated phase object rather than recreating
        // it for every single test case.
        pair<string, string> key = {setup["file"].asString(), setup.getString("phase", "")};
        if (cache.count(key) == 0) {
            cache[key] = newThermo(key.first, key.second);
        }
        atol = setup.getDouble("atol", 1e-5);
        rtol_fd = setup.getDouble("rtol_fd", 1e-6);
        atol_v = setup.getDouble("atol_v", 1e-11);
        atol_c = setup.getDouble("atol_c", 1e-14);
        atol_e = setup.getDouble("atol_e", 1e-18);

        phase = cache[key];
        phase->setState(state);
        nsp = phase->nSpecies();
        p = phase->pressure();
        T = phase->temperature();
        Te = phase->electronTemperature();
        RTe = Te * GasConstant;
        RT = T * GasConstant;
        if (phase->type() == "plasma") {
            ke = dynamic_cast<PlasmaPhase&>(*phase).electronSpeciesIndex();
        } else {
            ke = npos;
        }
    }

    ~TestConsistency() {
        make_deprecation_warnings_fatal();
    }

    void SetUp() override {
        // See if we should skip this test specific test case
        if (setup.hasKey("known-failures")) {
            auto current = testing::UnitTest::GetInstance()->current_test_info()->name();
            for (auto& [pattern, reason] : setup["known-failures"].asMap<string>()) {
                if (regex_search(current, regex(pattern))) {
                    GTEST_SKIP() << reason;
                    break;
                }
            }
        }
    }

    static map<pair<string, string>, shared_ptr<ThermoPhase>> cache;

    AnyMap setup;
    shared_ptr<ThermoPhase> phase;
    size_t nsp;
    double T, p, RT, RTe, Te;
    size_t ke;
    double atol; // absolute tolerance for molar energy comparisons
    double atol_v; // absolute tolerance for molar volume comparisons
    double atol_c; // absolute tolerance for compressibility comparison
    double atol_e; // absolute tolerance for expansion comparison
    double rtol_fd; // relative tolerance for finite difference comparisons
};

map<pair<string, string>, shared_ptr<ThermoPhase>> TestConsistency::cache = {};

// --------------- Definitions for individual consistency tests ---------------

TEST_P(TestConsistency, h_eq_u_plus_Pv) {
    double h, u, v;
    try {
        h = phase->enthalpy_mole();
        u = phase->intEnergy_mole();
        v = phase->molarVolume();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(h, u + p * v, atol);
}

TEST_P(TestConsistency, g_eq_h_minus_Ts) {
    double g, h, s;
    try {
        g = phase->gibbs_mole();
        h = phase->enthalpy_mole();
        s = phase->entropy_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(g, h - T * s, atol);
}

TEST_P(TestConsistency, hk_eq_uk_plus_P_vk)
{
    vector<double> hk(nsp), uk(nsp), vk(nsp);
    try {
        phase->getPartialMolarEnthalpies(hk.data());
        phase->getPartialMolarIntEnergies(uk.data());
        phase->getPartialMolarVolumes(vk.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(hk[k], uk[k] + p * vk[k], atol) << "k = " << k;
        } // not applicable for electron
    }
}

TEST_P(TestConsistency, gk_eq_hk_minus_T_sk)
{
    vector<double> gk(nsp), hk(nsp), sk(nsp);
    try {
        phase->getChemPotentials(gk.data());
        phase->getPartialMolarEnthalpies(hk.data());
        phase->getPartialMolarEntropies(sk.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(gk[k], hk[k] - T * sk[k], atol) << "k = " << k;
        } // not applicable for electron
    }
}

TEST_P(TestConsistency, h_eq_sum_hk_Xk)
{
    vector<double> hk(nsp);
    try {
        phase->getPartialMolarEnthalpies(hk.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(phase->enthalpy_mole(), phase->mean_X(hk), atol);
}

TEST_P(TestConsistency, u_eq_sum_uk_Xk)
{
    vector<double> uk(nsp);
    double u;
    try {
        phase->getPartialMolarIntEnergies(uk.data());
        u = phase->intEnergy_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(u, phase->mean_X(uk), atol);
}

TEST_P(TestConsistency, g_eq_sum_gk_Xk)
{
    vector<double> gk(nsp);
    double g;
    try {
        phase->getChemPotentials(gk.data());
        g = phase->gibbs_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(g, phase->mean_X(gk), atol);
}

TEST_P(TestConsistency, s_eq_sum_sk_Xk)
{
    vector<double> sk(nsp);
    double s;
    try {
        phase->getPartialMolarEntropies(sk.data());
        s = phase->entropy_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(s, phase->mean_X(sk), atol);
}

TEST_P(TestConsistency, v_eq_sum_vk_Xk)
{
    vector<double> vk(nsp);
    try {
        phase->getPartialMolarVolumes(vk.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(phase->molarVolume(), phase->mean_X(vk), atol_v);
}

TEST_P(TestConsistency, cp_eq_sum_cpk_Xk)
{
    vector<double> cpk(nsp);
    double cp;
    try {
        phase->getPartialMolarCp(cpk.data());
        cp = phase->cp_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(cp, phase->mean_X(cpk), atol);
}

TEST_P(TestConsistency, cp_eq_dhdT)
{
    double h1, cp1;
    try {
        h1 = phase->enthalpy_mole();
        cp1 = phase->cp_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    phase->setState_TP(T1 + dT, phase->pressure());
    double h2 = phase->enthalpy_mole();
    double cp2 = phase->cp_mole();
    double cp_mid = 0.5 * (cp1 + cp2);
    double cp_fd = (h2 - h1) / dT;
    EXPECT_NEAR(cp_fd, cp_mid, max({rtol_fd * cp_mid, rtol_fd * cp_fd, atol}));
}

TEST_P(TestConsistency, cv_eq_dudT)
{
    double u1, cv1;
    try {
        u1 = phase->intEnergy_mole();
        cv1 = phase->cv_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    if (phase->isCompressible()) {
        phase->setState_TD(T1 + dT, phase->density());
    } else {
        phase->setTemperature(T1 + dT);
    }
    double u2 = phase->intEnergy_mole();
    double cv2 = phase->cv_mole();
    double cv_mid = 0.5 * (cv1 + cv2);
    double cv_fd = (u2 - u1) / dT;
    EXPECT_NEAR(cv_fd, cv_mid, max({rtol_fd * cv_mid, rtol_fd * cv_fd, atol}));
}

TEST_P(TestConsistency, cp_eq_dsdT_const_p_times_T)
{
    double s1, cp1;
    try {
        s1 = phase->entropy_mole();
        cp1 = phase->cp_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-4 * phase->temperature();
    phase->setState_TP(T1 + dT, phase->pressure());
    double s2 = phase->entropy_mole();
    double cp2 = phase->cp_mole();
    double cp_mid = 0.5 * (cp1 + cp2);
    double cp_fd = (T1 + dT/2) * (s2 - s1) / dT;
    EXPECT_NEAR(cp_fd, cp_mid, max({rtol_fd * cp_mid, rtol_fd * cp_fd, atol}));
}

TEST_P(TestConsistency, cv_eq_dsdT_const_v_times_T)
{
    double s1, cv1;
    try {
        s1 = phase->entropy_mole();
        cv1 = phase->cv_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-4 * phase->temperature();
    if (phase->isCompressible()) {
        phase->setState_TD(T1 + dT, phase->density());
    } else {
        phase->setTemperature(T1 + dT);
    }
    double s2 = phase->entropy_mole();
    double cv2 = phase->cv_mole();
    double cv_mid = 0.5 * (cv1 + cv2);
    double cv_fd = (T1 + dT/2) * (s2 - s1) / dT;
    EXPECT_NEAR(cv_fd, cv_mid, max({rtol_fd * cv_mid, rtol_fd * cv_fd, atol}));
}

TEST_P(TestConsistency, dsdP_const_T_eq_minus_dV_dT_const_P)
{
    double P0 = phase->pressure();
    double T0 = phase->temperature();
    double s1, v1;
    double dP = 1e-4 * P0;
    double dT = 1e-4 * T0;
    try {
        phase->setState_TP(T0, P0 - dP);
        s1 = phase->entropy_mole();
        phase->setState_TP(T0 - dT, P0);
        v1 = phase->molarVolume();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    phase->setState_TP(T0, P0 + dP);
    double s2 = phase->entropy_mole();
    double dsdP = (s2 - s1) / (2 * dP);

    phase->setState_TP(T0 + dT, P0);
    double v2 = phase->molarVolume();
    double dvdT = (v2 - v1) / (2 * dT);
    double tol = rtol_fd * std::max(std::abs(dsdP), std::abs(dvdT));
    EXPECT_NEAR(dsdP, -dvdT, tol);
}

TEST_P(TestConsistency, dSdv_const_T_eq_dPdT_const_V) {
    if (phase->isCompressible()) {
        double s1;
        try {
            s1 = phase->entropy_mass();
        } catch (NotImplementedError& err) {
            GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
        }
        double v1 = 1 / phase->density();
        double P1 = phase->pressure();
        double v2 = v1 * (1 + 1e-7);
        phase->setState_TD(T, 1 / v2);
        double s2 = phase->entropy_mass();
        double dsdv = (s2 - s1) / (v2 - v1);

        double T2 = T * (1 + 1e-7);
        phase->setState_TD(T2, 1 / v1);
        double P2 = phase->pressure();
        double dPdT = (P2 - P1) / (T2 - T);
        double tol = rtol_fd * std::max(std::abs(dPdT), std::abs(dsdv));
        EXPECT_NEAR(dsdv, dPdT, tol);
    } else {
        GTEST_SKIP() << "Undefined for incompressible phase";
    }
}

TEST_P(TestConsistency, betaT_eq_minus_dmv_dP_const_T_div_mv)
{
    double betaT1;
    try {
        betaT1 = phase->isothermalCompressibility();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }

    double T = phase->temperature();
    double P1 = phase->pressure();
    double mv1 = phase->molarVolume();

    double P2 = P1 * (1 + 1e-6);
    phase->setState_TP(T, P2);
    double betaT2 = phase->isothermalCompressibility();
    double mv2 = phase->molarVolume();

    double betaT_mid = 0.5 * (betaT1 + betaT2);
    double mv_mid = 0.5 * (mv1 + mv2);
    double betaT_fd = -1 / mv_mid * (mv2 - mv1) / (P2 - P1);

    EXPECT_NEAR(betaT_fd, betaT_mid,
                max({rtol_fd * betaT_mid, rtol_fd * betaT_fd, atol_c}));
}

TEST_P(TestConsistency, alphaV_eq_dmv_dT_const_P_div_mv)
{
    double alphaV1;
    try {
        alphaV1 = phase->thermalExpansionCoeff();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }

    double P = phase->pressure();
    double T1 = phase->temperature();
    double mv1 = phase->molarVolume();

    double T2 = T1 * (1 + 1e-6);
    phase->setState_TP(T2, P);
    double alphaV2 = phase->thermalExpansionCoeff();
    double mv2 = phase->molarVolume();

    double alphaV_mid = 0.5 * (alphaV1 + alphaV2);
    double mv_mid = 0.5 * (mv1 + mv2);
    double alphaV_fd = 1 / mv_mid * (mv2 - mv1) / (T2 - T1);

    EXPECT_NEAR(alphaV_fd, alphaV_mid,
                max({rtol_fd * alphaV_mid, rtol_fd * alphaV_fd, atol_e}));
}

TEST_P(TestConsistency, c_eq_sqrt_dP_drho_const_s)
{
    double c1;
    try {
        c1 = phase->soundSpeed();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }

    double P1 = phase->pressure();
    double rho1 = phase->density();

    double rho2 = rho1 * (1 + 1e-6);
    phase->setState_SV(phase->entropy_mass(), 1 / rho2);
    double c2 = phase->soundSpeed();
    double P2 = phase->pressure();

    double c_mid = 0.5 * (c1 + c2);
    double c_fd = sqrt((P2 - P1) / (rho2 - rho1));

    EXPECT_NEAR(c_fd, c_mid, max({rtol_fd * c_mid, rtol_fd * c_fd, atol}));
}

// ---------- Tests for consistency of standard state properties ---------------

TEST_P(TestConsistency, hk0_eq_uk0_plus_p_vk0)
{
    vector<double> h0(nsp), u0(nsp), v0(nsp);
    try {
        phase->getEnthalpy_RT(h0.data());
        phase->getIntEnergy_RT(u0.data());
        phase->getStandardVolumes(v0.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(h0[k] * RT, u0[k] * RT + p * v0[k], atol) << "k = " << k;
    }
}

TEST_P(TestConsistency, gk0_eq_hk0_minus_T_sk0)
{
    vector<double> g0(nsp), h0(nsp), s0(nsp);
    try {
        phase->getEnthalpy_RT(h0.data());
        phase->getGibbs_RT(g0.data());
        phase->getEntropy_R(s0.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(g0[k] * RT ,
                    h0[k] * RT - T * s0[k] * GasConstant, atol) << "k = " << k;
    }
}

TEST_P(TestConsistency, cpk0_eq_dhk0dT)
{
    vector<double> h1(nsp), h2(nsp), cp1(nsp), cp2(nsp);
    try {
        phase->getEnthalpy_RT(h1.data());
        phase->getCp_R(cp1.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    phase->setState_TP(T1 + dT, phase->pressure());
    phase->getEnthalpy_RT(h2.data());
    phase->getCp_R(cp2.data());
    for (size_t k = 0; k < nsp; k++) {
        double cp_mid = 0.5 * (cp1[k] + cp2[k]) * GasConstant;
        double cp_fd = (h2[k] * (T1 + dT) - h1[k] * T1) / dT * GasConstant;
        double tol = max({rtol_fd * std::abs(cp_mid), rtol_fd * std::abs(cp_fd), atol});
        EXPECT_NEAR(cp_fd, cp_mid, tol) << "k = " << k;
    }
}

TEST_P(TestConsistency, standard_gibbs_nondim)
{
    vector<double> g0_RT(nsp), mu0(nsp);
    try {
        phase->getGibbs_RT(g0_RT.data());
        phase->getStandardChemPotentials(mu0.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(g0_RT[k] * RT , mu0[k], atol) << "k = " << k;
        } else {
            EXPECT_NEAR(g0_RT[k] * RTe , mu0[k], atol) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, chem_potentials_to_activities) {
    vector<double> mu0(nsp), mu(nsp), a(nsp);
    try {
        phase->getChemPotentials(mu.data());
        phase->getStandardChemPotentials(mu0.data());
        phase->getActivities(a.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            double a_from_mu = exp((mu[k] - mu0[k]) / RT);
            double scale = std::max(std::abs(a[k]), std::abs(a_from_mu));
            EXPECT_NEAR(a_from_mu, a[k], 1e-9 * scale + 1e-14) << "k = " << k;
        } else {
            double a_from_mu = exp((mu[k] - mu0[k]) / RTe);
            double scale = std::max(std::abs(a[k]), std::abs(a_from_mu));
            EXPECT_NEAR(a_from_mu, a[k], 1e-9 * scale + 1e-14) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, activity_coeffs) {
    vector<double> a(nsp), gamma(nsp), X(nsp);
    try {
        phase->getActivities(a.data());
        phase->getActivityCoefficients(gamma.data());
        phase->getMoleFractions(X.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        double scale = std::max(std::abs(a[k]), std::abs(gamma[k] * X[k]));
        EXPECT_NEAR(a[k], gamma[k] * X[k], 1e-10 * scale + 1e-20) << "k = " << k;
    }
}

TEST_P(TestConsistency, activity_concentrations) {
    vector<double> a(nsp), Cact(nsp);
    try {
        phase->getActivities(a.data());
        phase->getActivityConcentrations(Cact.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        double C0k = phase->standardConcentration(k);
        EXPECT_NEAR(a[k], Cact[k] / C0k, 1e-9 * std::abs(a[k])) << "k = " << k;
    }
}

TEST_P(TestConsistency, log_standard_concentrations) {
    for (size_t k = 0; k < nsp; k++) {
        double c0 = phase->standardConcentration(k);
        EXPECT_NEAR(c0, exp(phase->logStandardConc(k)), 1e-10 * c0) << "k = " << k;
    }
}

TEST_P(TestConsistency, log_activity_coeffs) {
    vector<double> gamma(nsp), log_gamma(nsp);
    try {
        phase->getActivityCoefficients(gamma.data());
        phase->getLnActivityCoefficients(log_gamma.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(gamma[k], exp(log_gamma[k]), 1e-10 * gamma[k]) << "k = " << k;
    }
}

// -------------------- Tests for reference state properties -------------------

TEST_P(TestConsistency, hRef_eq_uRef_plus_P_vRef)
{
    vector<double> hRef(nsp), uRef(nsp), vRef(nsp);
    try {
        phase->getEnthalpy_RT_ref(hRef.data());
        phase->getIntEnergy_RT_ref(uRef.data());
        phase->getStandardVolumes_ref(vRef.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(hRef[k] * RT, uRef[k] * RT + OneAtm * vRef[k], atol) << "k = " << k;
        } else {
            EXPECT_NEAR(hRef[k] * RTe, uRef[k] * RTe + OneAtm * vRef[k], atol) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, gRef_eq_hRef_minus_T_sRef)
{
    vector<double> hRef(nsp), gRef_RT(nsp), gRef(nsp), sRef(nsp);
    try {
        phase->getEnthalpy_RT_ref(hRef.data());
        phase->getGibbs_ref(gRef.data());
        phase->getGibbs_RT_ref(gRef_RT.data());
        phase->getEntropy_R_ref(sRef.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(gRef[k], gRef_RT[k] * RT, atol) << "k = " << k;
            EXPECT_NEAR(gRef[k], hRef[k] * RT - T * sRef[k] * GasConstant,
                        atol) << "k = " << k;
        } else {
            EXPECT_NEAR(gRef[k], gRef_RT[k] * RTe, atol) << "k = " << k;
            EXPECT_NEAR(gRef[k], hRef[k] * RTe - Te * sRef[k] * GasConstant,
                        atol) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, cpRef_eq_dhRefdT)
{
    vector<double> h1(nsp), h2(nsp), cp1(nsp), cp2(nsp);
    try {
        phase->getEnthalpy_RT_ref(h1.data());
        phase->getCp_R_ref(cp1.data());
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    phase->setState_TP(T1 + dT, phase->pressure());
    phase->getEnthalpy_RT_ref(h2.data());
    phase->getCp_R_ref(cp2.data());
    for (size_t k = 0; k < nsp; k++) {
        double cp_mid = 0.5 * (cp1[k] + cp2[k]) * GasConstant;
        double cp_fd = (h2[k] * (T1 + dT) - h1[k] * T1) / dT * GasConstant;
        double tol = max({rtol_fd * std::abs(cp_mid), rtol_fd * std::abs(cp_fd), atol});
        EXPECT_NEAR(cp_fd, cp_mid, tol) << "k = " << k;
    }
}

// ------- Instantiation of test suites defined in consistency-cases.yaml ------

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

INSTANTIATE_TEST_SUITE_P(DebyeHuckelDilute, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-dilute")),
        testing::ValuesIn(getStates("debye-huckel-dilute")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckelDilute_IAPWS, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-dilute-IAPWS")),
        testing::ValuesIn(getStates("debye-huckel-dilute-IAPWS")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_b_dot_ak, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-B-dot-ak")),
        testing::ValuesIn(getStates("debye-huckel-B-dot-ak")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_b_dot_ak_IAPWS, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-B-dot-ak-IAPWS")),
        testing::ValuesIn(getStates("debye-huckel-B-dot-ak-IAPWS")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_b_dot_a, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-B-dot-a")),
        testing::ValuesIn(getStates("debye-huckel-B-dot-a")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_b_dot_a_IAPWS, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-B-dot-a-IAPWS")),
        testing::ValuesIn(getStates("debye-huckel-B-dot-a-IAPWS")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_pitzer_beta_ij, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-pitzer-beta_ij")),
        testing::ValuesIn(getStates("debye-huckel-pitzer-beta_ij")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_pitzer_beta_ij_IAPWS, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-pitzer-beta_ij-IAPWS")),
        testing::ValuesIn(getStates("debye-huckel-pitzer-beta_ij-IAPWS")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_beta_ij, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-beta_ij")),
        testing::ValuesIn(getStates("debye-huckel-beta_ij")))
);

INSTANTIATE_TEST_SUITE_P(DebyeHuckel_beta_ij_IAPWS, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("debye-huckel-beta_ij-IAPWS")),
        testing::ValuesIn(getStates("debye-huckel-beta_ij-IAPWS")))
);

INSTANTIATE_TEST_SUITE_P(Margules, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("margules")),
        testing::ValuesIn(getStates("margules")))
);

INSTANTIATE_TEST_SUITE_P(FixedStoichiometry, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("fixed-stoichiometry")),
        testing::ValuesIn(getStates("fixed-stoichiometry")))
);

INSTANTIATE_TEST_SUITE_P(IdealSurface, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-surface")),
        testing::ValuesIn(getStates("ideal-surface")))
);

INSTANTIATE_TEST_SUITE_P(IdealEdge, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-edge")),
        testing::ValuesIn(getStates("ideal-edge")))
);

INSTANTIATE_TEST_SUITE_P(CoverageDependentSurface, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("coverage-dependent-surface")),
        testing::ValuesIn(getStates("coverage-dependent-surface")))
);

INSTANTIATE_TEST_SUITE_P(LiquidWaterIapws95, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("liquid-water-IAPWS95")),
        testing::ValuesIn(getStates("liquid-water-IAPWS95")))
);

INSTANTIATE_TEST_SUITE_P(IdealSolnVPSS_simple, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-solution-VPSS-simple")),
        testing::ValuesIn(getStates("ideal-solution-VPSS-simple")))
);

INSTANTIATE_TEST_SUITE_P(IdealSolnVPSS_HKFT, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("ideal-solution-VPSS-HKFT")),
        testing::ValuesIn(getStates("ideal-solution-VPSS-HKFT")))
);

INSTANTIATE_TEST_SUITE_P(RedlichKister_LiC6, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("Redlich-Kister-LiC6")),
        testing::ValuesIn(getStates("Redlich-Kister-LiC6")))
);

INSTANTIATE_TEST_SUITE_P(RedlichKister_complex, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("Redlich-Kister-complex")),
        testing::ValuesIn(getStates("Redlich-Kister-complex")))
);

INSTANTIATE_TEST_SUITE_P(HMWSoln, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("HMW-electrolyte")),
        testing::ValuesIn(getStates("HMW-electrolyte")))
);

int main(int argc, char** argv)
{
    printf("Running main() from consistency.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    Cantera::make_deprecation_warnings_fatal();
    Cantera::CanteraError::setStackTraceDepth(0);
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
