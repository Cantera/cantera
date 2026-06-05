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
#include <limits>
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
    if (phase->type() == "plasma") {
        // Whenever a temperature can be defined, the following relation holds:
        //   g_k = h_k - T_k * s_k
        // When multiplying by mole fraction and summing over all species, we get:
        //   sum_k X_k * g_k = sum_k X_k * h_k - sum_k X_k * T_k * s_k
        // The left side is simply g. The first term on the right side is h.
        // The second term on the right side can be separated into electron and
        // non-electron contributions (for a 2 temperature plasma model):
        //   sum_k X_k * T_k * s_k = T * sum_{k!=e} X_k * s_k + Te * X_e * s_e
        //                         = T * (s - X_e * s_e) + Te * X_e * s_e
        //                         = T * s + X_e * s_e * (Te - T)
        // Rearranging gives:
        //   g = h - T * s - X_e * s_e * (Te - T)

        // Get partial molar entropy of electron species.
        vector<double> sk(nsp);
        phase->getPartialMolarEntropies(sk);
        double se = sk[ke];
        // Get mole fraction of electron species.
        double Xe = phase->moleFraction(ke);

        EXPECT_NEAR(g, h - T * s - Xe * se * (Te - T), atol);
    }
    else {
        EXPECT_NEAR(g, h - T * s, atol);
    }
}

TEST_P(TestConsistency, hk_eq_uk_plus_P_vk)
{
    vector<double> hk(nsp), uk(nsp), vk(nsp);
    try {
        phase->getPartialMolarEnthalpies(hk);
        phase->getPartialMolarIntEnergies(uk);
        phase->getPartialMolarVolumes(vk);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(hk[k], uk[k] + p * vk[k], atol) << "k = " << k;
    }
}

TEST_P(TestConsistency, utilde_eq_uk_minus_piT_vk)
{
    vector<double> utilde(nsp), uk(nsp), vk(nsp);
    double piT;
    try {
        phase->getPartialMolarIntEnergies_TV(utilde);
        phase->getPartialMolarIntEnergies(uk);
        phase->getPartialMolarVolumes(vk);
        piT = phase->internalPressure();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    if (!std::isfinite(piT)) {
        GTEST_SKIP() << "internalPressure is not finite at this state";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            double expected = uk[k] - piT * vk[k];
            double tol = max({atol, rtol_fd * abs(utilde[k]), rtol_fd * abs(expected)});
            EXPECT_NEAR(utilde[k], expected, tol) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, piT_eq_dudv_const_T)
{
    if (!phase->isCompressible()) {
        GTEST_SKIP() << "Undefined for incompressible phase";
    }

    double piT;
    double T0 = phase->temperature();
    double v0 = 1.0 / phase->density(); // specific volume [m^3/kg]
    // Use a moderate perturbation to avoid cancellation in u(v) differences.
    double dv = 1e-4 * v0;
    double v1 = v0 - dv;
    double v2 = v0 + dv;

    try {
        piT = phase->internalPressure();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }

    phase->setState_TD(T0, 1.0 / v1);
    double u1 = phase->intEnergy_mass();
    phase->setState_TD(T0, 1.0 / v2);
    double u2 = phase->intEnergy_mass();

    double piT_fd = (u2 - u1) / (v2 - v1);
    double tol = max({rtol_fd * abs(piT), rtol_fd * abs(piT_fd), atol});
    EXPECT_NEAR(piT_fd, piT, tol);
}

TEST_P(TestConsistency, gk_eq_hk_minus_T_sk)
{
    vector<double> gk(nsp), hk(nsp), sk(nsp);
    try {
        phase->getChemPotentials(gk);
        phase->getPartialMolarEnthalpies(hk);
        phase->getPartialMolarEntropies(sk);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(gk[k], hk[k] - T * sk[k], atol) << "k = " << k;
        }
        else {
            EXPECT_NEAR(gk[k], hk[k] - Te * sk[k], atol) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, h_eq_sum_hk_Xk)
{
    vector<double> hk(nsp);
    try {
        phase->getPartialMolarEnthalpies(hk);
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
        phase->getPartialMolarIntEnergies(uk);
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
        phase->getChemPotentials(gk);
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
        phase->getPartialMolarEntropies(sk);
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
        phase->getPartialMolarVolumes(vk);
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
        phase->getPartialMolarCp(cpk);
        cp = phase->cp_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    EXPECT_NEAR(cp, phase->mean_X(cpk), atol);
}

TEST_P(TestConsistency, gibbs_duhem_const_T_P)
{
    // Gibbs-Duhem at constant T and P requires sum_k X_k * dmu_k = 0 for any
    // change in composition. This is the differential complement to the Euler
    // relation g = sum_k X_k * mu_k (tested by g_eq_sum_gk_Xk): together they
    // are equivalent to the fundamental relation
    // dU = T dS - P dV + sum_k mu_k dN_k. A model can report chemical potentials
    // that satisfy the Euler relation at every individual state yet still violate
    // this differential identity if the mu_k are not consistent partial molar
    // Gibbs energies of a single underlying G(T, P, N).
    vector<double> X0(nsp), mu_plus(nsp), mu_minus(nsp);
    double T0 = phase->temperature();
    double P0 = phase->pressure();
    phase->getMoleFractions(X0);

    // Only species that are actually present can be perturbed.
    vector<size_t> active;
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke && X0[k] > 1e-6) {
            active.push_back(k);
        }
    }
    if (active.size() < 2) {
        GTEST_SKIP() << "Fewer than two species available to perturb composition";
    }

    // Test each independent composition direction by transferring a small amount
    // of material from species active[0] to each of the other active species,
    // using a centered difference about the original composition X0. The step is
    // scaled to the smaller of the two mole fractions so that the relative
    // perturbation (and hence the finite-difference truncation error) is the same
    // regardless of how dilute the perturbed species are.
    double eps = 1e-4;
    size_t i = active[0];
    for (size_t n = 1; n < active.size(); n++) {
        size_t j = active[n];
        double dx = eps * std::min(X0[i], X0[j]);
        vector<double> Xp = X0, Xm = X0;
        Xp[i] -= dx; Xp[j] += dx;
        Xm[i] += dx; Xm[j] -= dx;

        try {
            phase->setMoleFractions(Xp);
            phase->setState_TP(T0, P0);
            phase->getChemPotentials(mu_plus);
            phase->setMoleFractions(Xm);
            phase->setState_TP(T0, P0);
            phase->getChemPotentials(mu_minus);
        } catch (NotImplementedError& err) {
            GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
        }

        double gd = 0.0;
        double scale = 0.0;
        for (size_t k = 0; k < nsp; k++) {
            double term = X0[k] * (mu_plus[k] - mu_minus[k]);
            gd += term;
            scale = std::max(scale, std::abs(term));
        }
        EXPECT_NEAR(gd, 0.0, rtol_fd * scale + atol)
            << "perturbing species " << i << " <-> " << j;
    }

    // Restore the original state for any subsequent use of the cached phase.
    phase->setMoleFractions(X0);
    phase->setState_TP(T0, P0);
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
    if (phase->type() == "plasma") {
        phase->setElectronTemperature(T1 + dT);
    }
    phase->setState_TP(T1 + dT, phase->pressure());
    double h2 = phase->enthalpy_mole();
    double cp2 = phase->cp_mole();
    double cp_mid = 0.5 * (cp1 + cp2);
    double cp_fd = (h2 - h1) / dT;
    EXPECT_NEAR(cp_fd, cp_mid, max({rtol_fd * cp_mid, rtol_fd * cp_fd, atol}));
}

TEST_P(TestConsistency, cv_eq_dudT)
{
    if (!phase->isCompressible()) {
        // For nearly-incompressible phases, we can't use setState_TD to evaluate
        // ∂u/∂T at constant V using finite differences to compare with cv.
        GTEST_SKIP() << "Not meaningful for incompressible phases";
    }
    double u1, cv1;
    try {
        u1 = phase->intEnergy_mole();
        cv1 = phase->cv_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    if (phase->type() == "plasma") {
        phase->setElectronTemperature(T1 + dT);
    }
    phase->setState_TD(T1 + dT, phase->density());
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
    double pressure = phase->pressure();
    if (phase->type() == "plasma") {
        phase->setElectronTemperature(T1 + dT);
    }
    phase->setState_TP(T1 + dT, pressure);
    double s2 = phase->entropy_mole();
    double cp2 = phase->cp_mole();
    double cp_mid = 0.5 * (cp1 + cp2);
    double cp_fd = (T1 + dT/2) * (s2 - s1) / dT;
    EXPECT_NEAR(cp_fd, cp_mid, max({rtol_fd * cp_mid, rtol_fd * cp_fd, atol}));
}

TEST_P(TestConsistency, cv_eq_dsdT_const_v_times_T)
{
    if (!phase->isCompressible()) {
        // For nearly-incompressible phases, we can't use setState_TD to evaluate ∂S/∂T
        // at constant V by finite difference to compare with cv.
        GTEST_SKIP() << "Not meaningful for incompressible phases";
    }
    double s1, cv1;
    try {
        s1 = phase->entropy_mole();
        cv1 = phase->cv_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    double T1 = phase->temperature();
    double dT = 1e-4 * phase->temperature();
    if (phase->type() == "plasma") {
        phase->setElectronTemperature(T1 + dT);
    }
    phase->setState_TD(T1 + dT, phase->density());
    double s2 = phase->entropy_mole();
    double cv2 = phase->cv_mole();
    double cv_mid = 0.5 * (cv1 + cv2);
    double cv_fd = (T1 + dT/2) * (s2 - s1) / dT;
    EXPECT_NEAR(cv_fd, cv_mid, max({rtol_fd * cv_mid, rtol_fd * cv_fd, atol}));
}

TEST_P(TestConsistency, cv_eq_cp_minus_T_V_beta2_over_kappa_T)
{
    double cv, cp, beta, kappa_T, V;
    try {
        cv = phase->cv_mole();
        cp = phase->cp_mole();
        beta = phase->thermalExpansionCoeff();
        kappa_T = phase->isothermalCompressibility();
        V = phase->molarVolume();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    if (kappa_T == 0.0) {
        // Truly incompressible: cv = cp exactly
        EXPECT_NEAR(cv, cp, atol);
    } else {
        EXPECT_NEAR(cv, cp - T * V * beta * beta / kappa_T, max({rtol_fd * cv, atol}));
    }
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
        if (phase->type() == "plasma") {
            phase->setElectronTemperature(T0 - dT);
        }
        phase->setState_TP(T0 - dT, P0);
        v1 = phase->molarVolume();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }

    if (phase->type() == "plasma") {
            phase->setElectronTemperature(T0);
    }
    phase->setState_TP(T0, P0 + dP);
    double s2 = phase->entropy_mole();
    double dsdP = (s2 - s1) / (2 * dP);

    if (phase->type() == "plasma") {
        phase->setElectronTemperature(T0 + dT);
    }
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
        if (phase->type() == "plasma") {
            phase->setElectronTemperature(T2);
        }
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
        max({rtol_fd * std::abs(alphaV_mid), rtol_fd * std::abs(alphaV_fd), atol_e}));
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

TEST_P(TestConsistency, sk_eq_minus_dmu_k_dT_const_P_X)
{
    // Composition<->T Maxwell relation: partial molar entropy s_k equals
    // -(d mu_k / dT) at constant P and composition. The s_k are read at the
    // midpoint T0 and compared with a central finite difference of mu_k.
    vector<double> s_k(nsp), mu1(nsp), mu2(nsp);
    double T0 = phase->temperature();
    double P0 = phase->pressure();
    double dT = 1e-5 * T0;
    try {
        phase->getPartialMolarEntropies(s_k);
        phase->setState_TP(T0 - dT, P0);
        phase->getChemPotentials(mu1);
        phase->setState_TP(T0 + dT, P0);
        phase->getChemPotentials(mu2);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    phase->setState_TP(T0, P0); // restore cached state
    for (size_t k = 0; k < nsp; k++) {
        if (k == ke) {
            continue; // electron species handled by two-temperature identities
        }
        double dmu_dT = (mu2[k] - mu1[k]) / (2 * dT);
        EXPECT_NEAR(-dmu_dT, s_k[k],
            max({rtol_fd * s_k[k], rtol_fd * std::abs(dmu_dT), atol}))
            << "for species " << k;
    }
}

TEST_P(TestConsistency, vk_eq_dmu_k_dP_const_T_X)
{
    // Composition<->P Maxwell relation: partial molar volume v_k equals
    // (d mu_k / dP) at constant T and composition.
    vector<double> v_k(nsp), mu1(nsp), mu2(nsp);
    double T0 = phase->temperature();
    double P0 = phase->pressure();
    double dP = 1e-4 * P0;
    try {
        phase->getPartialMolarVolumes(v_k);
        phase->setState_TP(T0, P0 - dP);
        phase->getChemPotentials(mu1);
        phase->setState_TP(T0, P0 + dP);
        phase->getChemPotentials(mu2);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    phase->setState_TP(T0, P0); // restore cached state
    for (size_t k = 0; k < nsp; k++) {
        if (k == ke) {
            continue;
        }
        double dmu_dP = (mu2[k] - mu1[k]) / (2 * dP);
        EXPECT_NEAR(dmu_dP, v_k[k],
            max({rtol_fd * std::abs(v_k[k]), rtol_fd * std::abs(dmu_dP), atol_v}))
            << "for species " << k;
    }
}

TEST_P(TestConsistency, v_eq_dgdP_const_T)
{
    // dG pressure coefficient: molar volume equals (dG/dP) at constant T.
    double T0 = phase->temperature();
    double P0 = phase->pressure();
    double v, g1, g2;
    double dP = 1e-4 * P0;
    try {
        v = phase->molarVolume();
        phase->setState_TP(T0, P0 - dP);
        g1 = phase->gibbs_mole();
        phase->setState_TP(T0, P0 + dP);
        g2 = phase->gibbs_mole();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    phase->setState_TP(T0, P0); // restore cached state
    double dgdP = (g2 - g1) / (2 * dP);
    double v_expected = v;
    if (phase->type() == "plasma" && ke != npos) {
        // For a two-temperature plasma, the electron species contributes an extra
        // term because the pressure term in mu_e uses the electron temperature Te,
        // giving (dG/dP)_T = V + X_e * R * (Te - T) / P (compare with h_eq_u_plus_Pv).
        vector<double> X(nsp);
        phase->getMoleFractions(X);
        v_expected += X[ke] * (RTe - RT) / P0;
    }
    EXPECT_NEAR(dgdP, v_expected,
        max({rtol_fd * std::abs(v_expected), rtol_fd * std::abs(dgdP), atol_v}));
}

TEST_P(TestConsistency, dhdP_const_T_eq_v_minus_T_dvdT_const_P)
{
    // dH pressure coefficient: (dH/dP)_T = V - T (dV/dT)_P. Both the dH/dP and
    // the dV/dT derivatives use centered differences about (T0, P0).
    double T0 = phase->temperature();
    double P0 = phase->pressure();
    double dP = 1e-4 * P0;
    double dT = 1e-4 * T0;
    double v, h1, v_lo;
    try {
        v = phase->molarVolume();
        phase->setState_TP(T0, P0 - dP);
        h1 = phase->enthalpy_mole();
        phase->setState_TP(T0 - dT, P0);
        v_lo = phase->molarVolume();
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    phase->setState_TP(T0, P0 + dP);
    double h2 = phase->enthalpy_mole();
    double dhdP = (h2 - h1) / (2 * dP);

    phase->setState_TP(T0 + dT, P0);
    double v_hi = phase->molarVolume();
    double dvdT = (v_hi - v_lo) / (2 * dT);

    phase->setState_TP(T0, P0); // restore cached state
    double rhs = v - T0 * dvdT;
    EXPECT_NEAR(dhdP, rhs,
        max({rtol_fd * std::abs(dhdP), rtol_fd * std::abs(rhs), atol_v}));
}

TEST_P(TestConsistency, dmu_k_dNj_eq_dmu_j_dNk_const_T_P)
{
    // Symmetry of the composition Hessian of G: d mu_k / dN_j = d mu_j / dN_k.
    // This is the pure-composition Maxwell relation, equivalent to the mu_k
    // being the gradient of a single Gibbs energy G(T, P, N). It is only
    // nontrivial with at least three independently variable species; for a
    // binary solution, the single composition direction is already covered by
    // gibbs_duhem_const_T_P. Mole fractions are treated as mole numbers with
    // total N = 1; setMoleFractions renormalizes, so perturbing one entry is
    // equivalent to adding that species to an open system.
    //
    // Each mu_k carries a large, composition-independent standard-state term,
    // so differencing mu_k has a roundoff floor of order |mu|_max * eps_machine.
    // Divided by the perturbation dN this gives a noise floor on each computed
    // derivative, which is added to the comparison tolerance. Species too dilute to
    // perturb meaningfully thus are given a large floor to avoid producing false
    // positives.
    vector<double> X0(nsp);
    double T0 = phase->temperature();
    double P0 = phase->pressure();
    phase->getMoleFractions(X0);

    vector<size_t> active;
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke && X0[k] > 1e-6) {
            active.push_back(k);
        }
    }
    if (active.size() < 3) {
        GTEST_SKIP() << "Fewer than three species: Hessian symmetry is trivial";
    }

    double eps = 1e-4;
    double mu_mag = 0.0;
    vector<double> dn(active.size());

    // cols[a][k] = d mu_k / dN_{active[a]} via centered difference in N_{active[a]}
    auto dmu_dN = [&](size_t a, vector<double>& col) {
        size_t j = active[a];
        dn[a] = eps * X0[j];
        vector<double> np = X0, nm = X0, mu_p(nsp), mu_m(nsp);
        np[j] += dn[a];
        nm[j] -= dn[a];
        phase->setMoleFractions(np);
        phase->setState_TP(T0, P0);
        phase->getChemPotentials(mu_p);
        phase->setMoleFractions(nm);
        phase->setState_TP(T0, P0);
        phase->getChemPotentials(mu_m);
        for (size_t k = 0; k < nsp; k++) {
            col[k] = (mu_p[k] - mu_m[k]) / (2 * dn[a]);
            mu_mag = std::max({mu_mag, std::abs(mu_p[k]), std::abs(mu_m[k])});
        }
    };

    try {
        vector<vector<double>> cols(active.size(), vector<double>(nsp));
        for (size_t a = 0; a < active.size(); a++) {
            dmu_dN(a, cols[a]);
        }
        double eps_mach = std::numeric_limits<double>::epsilon();
        for (size_t a = 0; a < active.size(); a++) {
            for (size_t b = a + 1; b < active.size(); b++) {
                double kj = cols[b][active[a]]; // d mu_{active[a]} / dN_{active[b]}
                double jk = cols[a][active[b]]; // d mu_{active[b]} / dN_{active[a]}
                double scale = std::max(std::abs(kj), std::abs(jk));
                double noise = 10 * mu_mag * eps_mach
                               * std::max(1.0 / dn[a], 1.0 / dn[b]);
                EXPECT_NEAR(kj, jk, rtol_fd * scale + noise + atol)
                    << "species " << active[a] << " <-> " << active[b];
            }
        }
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }

    phase->setMoleFractions(X0);
    phase->setState_TP(T0, P0); // restore cached state
}

// ---------- Tests for consistency of standard state properties ---------------

TEST_P(TestConsistency, hk0_eq_uk0_plus_p_vk0)
{
    vector<double> h0(nsp), u0(nsp), v0(nsp);
    try {
        phase->getEnthalpy_RT(h0);
        phase->getIntEnergy_RT(u0);
        phase->getStandardVolumes(v0);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    for (size_t k = 0; k < nsp; k++) {
        if (k != ke) {
            EXPECT_NEAR(h0[k] * RT, u0[k] * RT + p * v0[k], atol) << "k = " << k;
        } else {
            EXPECT_NEAR(h0[k] * RTe, u0[k] * RTe + p * v0[k], atol) << "k = " << k;
        }
    }
}

TEST_P(TestConsistency, gk0_eq_hk0_minus_T_sk0)
{
    vector<double> g0(nsp), h0(nsp), s0(nsp);
    try {
        phase->getEnthalpy_RT(h0);
        phase->getGibbs_RT(g0);
        phase->getEntropy_R(s0);
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
        phase->getEnthalpy_RT(h1);
        phase->getCp_R(cp1);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    // Perturb temperature.
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    // For plasma phases, also perturb the electron temperature.
    double Te = phase->electronTemperature();
    double dTe = 1e-5 * Te;

    phase->setState_TP(T1 + dT, phase->pressure());
    if (phase->type() == "plasma") {
        phase->setElectronTemperature(Te + dTe);
    }
    phase->getEnthalpy_RT(h2);
    phase->getCp_R(cp2);

    for (size_t k = 0; k < nsp; k++) {
        // Determine effective temperature and perturbation for species k.
        double T_eff = (k == ke) ? Te : T1;
        double dT_eff = (k == ke) ? dTe : dT;

        double cp_mid = 0.5 * (cp1[k] + cp2[k]) * GasConstant;
        double cp_fd = (h2[k] * (T_eff + dT_eff) - h1[k] * T_eff) / dT_eff * GasConstant;
        double tol = max({rtol_fd * std::abs(cp_mid), rtol_fd * std::abs(cp_fd), atol});
        EXPECT_NEAR(cp_fd, cp_mid, tol) << "k = " << k;
    }
}

TEST_P(TestConsistency, standard_gibbs_nondim)
{
    vector<double> g0_RT(nsp), mu0(nsp);
    try {
        phase->getGibbs_RT(g0_RT);
        phase->getStandardChemPotentials(mu0);
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
        phase->getChemPotentials(mu);
        phase->getStandardChemPotentials(mu0);
        phase->getActivities(a);
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
        phase->getActivities(a);
        phase->getActivityCoefficients(gamma);
        phase->getMoleFractions(X);
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
        phase->getActivities(a);
        phase->getActivityConcentrations(Cact);
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
        phase->getActivityCoefficients(gamma);
        phase->getLnActivityCoefficients(log_gamma);
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
        phase->getEnthalpy_RT_ref(hRef);
        phase->getIntEnergy_RT_ref(uRef);
        phase->getStandardVolumes_ref(vRef);
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
        phase->getEnthalpy_RT_ref(hRef);
        phase->getGibbs_ref(gRef);
        phase->getGibbs_RT_ref(gRef_RT);
        phase->getEntropy_R_ref(sRef);
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
        phase->getEnthalpy_RT_ref(h1);
        phase->getCp_R_ref(cp1);
    } catch (NotImplementedError& err) {
        GTEST_SKIP() << err.getMethod() << " threw NotImplementedError";
    }
    // Perturb temperature.
    double T1 = phase->temperature();
    double dT = 1e-5 * phase->temperature();
    // For plasma phases, also perturb the electron temperature.
    double Te = phase->electronTemperature();
    double dTe = 1e-5 * Te;

    phase->setState_TP(T1 + dT, phase->pressure());
    if (phase->type() == "plasma") {
        phase->setElectronTemperature(Te + dTe);
    }
    phase->getEnthalpy_RT_ref(h2);
    phase->getCp_R_ref(cp2);

    for (size_t k = 0; k < nsp; k++) {
        // Determine effective temperature and perturbation for species k.
        double T_eff = (k == ke) ? Te : T1;
        double dT_eff = (k == ke) ? dTe : dT;

        double cp_mid = 0.5 * (cp1[k] + cp2[k]) * GasConstant;
        double cp_fd = (h2[k] * (T_eff + dT_eff) - h1[k] * T_eff) / dT_eff * GasConstant;
        double tol = max({rtol_fd * std::abs(cp_mid), rtol_fd * std::abs(cp_fd), atol});
        EXPECT_NEAR(cp_fd, cp_mid, tol) << "k = " << k << "(ke = " << ke << ")";
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

INSTANTIATE_TEST_SUITE_P(MargulesSSVol, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("margules-SSVol")),
        testing::ValuesIn(getStates("margules-SSVol")))
);

INSTANTIATE_TEST_SUITE_P(MargulesWithExcessVolume, TestConsistency,
    testing::Combine(
        testing::Values(getSetup("margules-with-excess-volume")),
        testing::ValuesIn(getStates("margules-with-excess-volume")))
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
    Cantera::addDataDirectory("test/data");
    Cantera::addDataDirectory("data");
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
