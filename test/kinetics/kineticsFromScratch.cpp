#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/BulkKinetics.h"
#include "cantera/kinetics/Reaction.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/Arrhenius.h"
#include "cantera/kinetics/ChebyshevRate.h"
#include "cantera/kinetics/Custom.h"
#include "cantera/kinetics/Falloff.h"
#include "cantera/kinetics/InterfaceRate.h"
#include "cantera/kinetics/PlogRate.h"
#include "cantera/kinetics/TwoTempPlasmaRate.h"
#include "cantera/base/Array.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

class KineticsFromScratch : public testing::Test
{
public:
    KineticsFromScratch()
        : pp(newThermo("../data/kineticsfromscratch.yaml"))
        , pp_ref(newThermo("../data/kineticsfromscratch.yaml"))
    {
        vector<shared_ptr<ThermoPhase>> th;
        th.push_back(pp_ref);
        kin_ref = newKinetics(th, "../data/kineticsfromscratch.yaml");

        kin.addThermo(pp);
        kin.init();
    }

    shared_ptr<ThermoPhase> pp;
    shared_ptr<ThermoPhase> pp_ref;
    BulkKinetics kin;
    shared_ptr<Kinetics> kin_ref;

    //! iRef is the index of the corresponding reaction in the reference mech
    void check_rates(int iRef) {
        ASSERT_EQ((size_t) 1, kin.nReactions());

        string X = "O:0.02 H2:0.2 O2:0.5 H:0.03 OH:0.05 H2O:0.1 HO2:0.01";
        pp->setState_TPX(1200, 5*OneAtm, X);
        pp_ref->setState_TPX(1200, 5*OneAtm, X);

        vector<double> k(1), k_ref(kin_ref->nReactions());

        kin.getFwdRateConstants(&k[0]);
        kin_ref->getFwdRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);

        kin.getRevRateConstants(&k[0]);
        kin_ref->getRevRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);
    }
};

TEST_F(KineticsFromScratch, add_elementary_reaction1)
{
    // reaction 0:
    //     equation: O + H2 <=> H + OH  # Reaction 1
    //     rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0}
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 2.619184e+07);
    auto R = make_shared<Reaction>(reac, prod, rate);

    kin.addReaction(R);
    check_rates(0);
}

TEST_F(KineticsFromScratch, add_elementary_reaction2)
{
    // reaction 0:
    //     equation: O + H2 <=> H + OH  # Reaction 1
    //     rate-constant: {A: 38.7, b: 2.7, Ea: 6260.0}
    string equation = "O + H2 <=> H + OH";
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 2.619184e+07);
    auto R = make_shared<Reaction>(equation, rate);

    kin.addReaction(R);
    check_rates(0);
}

TEST_F(KineticsFromScratch, add_three_body_reaction1)
{
    // reaction 1:
    //     equation: 2 O + M <=> O2 + M  # Reaction 2
    //     type: three-body
    //     rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0}
    //     efficiencies: {AR: 0.83, H2: 2.4, H2O: 15.4}
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("AR:0.83 H2:2.4 H2O:15.4");
    auto R = make_shared<Reaction>(reac, prod, rate, tbody);

    kin.addReaction(R);
    check_rates(1);

    reac = parseCompString("O:2, M:1");
    prod = parseCompString("O2:1, M:1");
    ASSERT_THROW(make_shared<Reaction>(reac, prod, rate, tbody), CanteraError);

}

TEST_F(KineticsFromScratch, add_three_body_reaction2)
{
    // reaction 1:
    //     equation: 2 O + M <=> O2 + M  # Reaction 2
    //     type: three-body
    //     rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0}
    //     efficiencies: {AR: 0.83, H2: 2.4, H2O: 15.4}
    string equation = "2 O + M <=> O2 + M";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("AR:0.83 H2:2.4 H2O:15.4");
    auto R = make_shared<Reaction>(equation, rate, tbody);
    auto reac = R->reactants;
    EXPECT_EQ(reac.count("M"), (size_t) 0);

    kin.addReaction(R);
    check_rates(1);
}

TEST_F(KineticsFromScratch, add_three_body_reaction3)
{
    string equation = "2 O + M <=> O2 + M";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto R = make_shared<Reaction>(equation, rate);
    EXPECT_TRUE(R->usesThirdBody());
}

TEST_F(KineticsFromScratch, multiple_third_bodies1)
{
    string equation = "2 H + 2 O2 <=> H2 + 2 O2";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto R = make_shared<Reaction>(equation, rate);
    EXPECT_FALSE(R->usesThirdBody());
}

TEST_F(KineticsFromScratch, multiple_third_bodies2)
{
    string equation = "2 H + 2 M <=> H2 + 2 M";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    ASSERT_THROW(Reaction(equation, rate), CanteraError);
}

TEST_F(KineticsFromScratch, multiple_third_bodies3)
{
    string equation = "2 H + O2 + M <=> H2 + O2 + M";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto R = make_shared<Reaction>(equation, rate);
    EXPECT_TRUE(R->usesThirdBody());
}

TEST_F(KineticsFromScratch, multiple_third_bodies4)
{
    string equation = "H2 + O2 => H2 + O2";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>("O2");
    auto R = make_shared<Reaction>(equation, rate, tbody);
    EXPECT_EQ(R->thirdBody()->name(), "O2");
    EXPECT_EQ(R->thirdBody()->default_efficiency, 0.);
    EXPECT_EQ(R->reactants.count("H2"), (size_t) 1);
    EXPECT_EQ(R->reactants.count("O2"), (size_t) 0);
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);

    AnyMap input = R->parameters(false);
    EXPECT_FALSE(input.hasKey("type"));
    EXPECT_TRUE(input.hasKey("efficiencies"));
    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 1u);
    EXPECT_EQ(efficiencies.begin()->first, "O2");
    EXPECT_FALSE(input.hasKey("default-efficiency"));
}

TEST_F(KineticsFromScratch, multiple_third_bodies5)
{
    string equation = "H2 + O2 => H2 + O2";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>("AR"); // incompatible third body
    ASSERT_THROW(Reaction(equation, rate, tbody), CanteraError);
}

TEST_F(KineticsFromScratch, multiple_third_bodies6)
{
    Composition reac = parseCompString("H2:1");
    Composition prod = parseCompString("H2:1");
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>("O2");
    auto R = make_shared<Reaction>(reac, prod, rate, tbody);
    EXPECT_EQ(R->thirdBody()->name(), "O2");
    EXPECT_EQ(R->thirdBody()->default_efficiency, 0.);
    EXPECT_EQ(R->reactants.count("H2"), (size_t) 1);
    EXPECT_EQ(R->reactants.count("O2"), (size_t) 0);
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);

    AnyMap input = R->parameters(false);
    EXPECT_FALSE(input.hasKey("type"));
    EXPECT_TRUE(input.hasKey("efficiencies"));
    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 1u);
    EXPECT_EQ(efficiencies.begin()->first, "O2");
    EXPECT_FALSE(input.hasKey("default-efficiency"));
}

TEST_F(KineticsFromScratch, multiple_third_bodies7)
{
    string equation = "CH2OCH + M <=> CH2CHO + M";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>("O2");
    auto R = make_shared<Reaction>(equation, rate, tbody);
    EXPECT_EQ(R->thirdBody()->name(), "M");

    AnyMap input = R->parameters(false);
    EXPECT_FALSE(input.hasKey("type"));
    EXPECT_TRUE(input.hasKey("efficiencies"));
    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 1u);
    EXPECT_EQ(efficiencies.begin()->first, "O2");
    EXPECT_TRUE(input.hasKey("default-efficiency"));
    EXPECT_EQ(input["default-efficiency"].asDouble(), 0.);
}

TEST_F(KineticsFromScratch, multiple_third_bodies8)
{
    string equation = "CH2OCH + M <=> CH2CHO + M";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("O2:1");
    tbody->default_efficiency = 0.;
    auto R = make_shared<Reaction>(equation, rate, tbody);
    EXPECT_EQ(R->thirdBody()->name(), "M");

    AnyMap input = R->parameters(false);
    EXPECT_FALSE(input.hasKey("type"));
    EXPECT_TRUE(input.hasKey("efficiencies"));
    auto efficiencies = input["efficiencies"].asMap<double>();
    EXPECT_EQ(efficiencies.size(), 1u);
    EXPECT_EQ(efficiencies.begin()->first, "O2");
    EXPECT_TRUE(input.hasKey("default-efficiency"));
    EXPECT_EQ(input["default-efficiency"].asDouble(), 0.);
}

TEST_F(KineticsFromScratch, multiple_third_bodies9)
{
    Composition reac = parseCompString("H2:1, O2:1");
    Composition prod = parseCompString("H2:1, O:2");
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>("O2");
    auto R = make_shared<Reaction>(reac, prod, rate, tbody);
    EXPECT_EQ(R->thirdBody()->name(), "O2");
    EXPECT_EQ(R->thirdBody()->default_efficiency, 0.);
    EXPECT_EQ(R->reactants.count("H2"), (size_t) 1);
    EXPECT_EQ(R->reactants.count("O2"), (size_t) 1);
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);

    reac = parseCompString("H2:1, O2:1");
    prod = parseCompString("H2:1, O2:1");
    ASSERT_THROW(make_shared<Reaction>(reac, prod, rate, tbody), CanteraError);
}

TEST_F(KineticsFromScratch, add_two_temperature_plasma)
{
    string equation = "O + H => O + H";
    auto rate = make_shared<TwoTempPlasmaRate>(17283, -3.1, -5820000, 1081000);
    auto R = make_shared<Reaction>(equation, rate);
    EXPECT_FALSE(R->usesThirdBody());
}

TEST_F(KineticsFromScratch, undefined_third_body)
{
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("H2:0.1 CO2:0.83");
    auto R = make_shared<Reaction>(reac, prod, rate, tbody);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
}

TEST_F(KineticsFromScratch, skip_undefined_third_body)
{
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("H2:0.1 CO2:0.83");
    auto R = make_shared<Reaction>(reac, prod, rate, tbody);

    kin.skipUndeclaredThirdBodies(true);
    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch, skip_explicit_third_body)
{
    string equation = "2 O + CO2 <=> O2 + CO2";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto R = make_shared<Reaction>(equation, rate);
    EXPECT_EQ(R->thirdBody()->name(), "CO2");

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    kin.skipUndeclaredThirdBodies(true);
    kin.addReaction(R);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch, third_body_composition)
{
    string equation = "2 O + H2O <=> O2 + H2O";
    auto rate = make_shared<ArrheniusRate>(1.2e11, -1.0, 0.0);
    auto R = make_shared<Reaction>(equation, rate);
    EXPECT_EQ(R->thirdBody()->name(), "H2O");
    EXPECT_TRUE(R->thirdBody()->explicit_3rd);

    Composition reac = R->reactants;
    EXPECT_EQ(reac.count("H2O"), (size_t) 0);
    EXPECT_EQ(reac.count("M"), (size_t) 0);
    Composition prod = R->products;
    EXPECT_EQ(prod.count("H2O"), (size_t) 0);
    EXPECT_EQ(reac.count("M"), (size_t) 0);
}

TEST_F(KineticsFromScratch, add_falloff_reaction1)
{
    // reaction 2:
    //     equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
    //     duplicate: true
    //     type: falloff
    //     low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0}
    //     high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0}
    //     Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
    //     efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
    Composition reac = parseCompString("OH:2");
    Composition prod = parseCompString("H2O2:1");
    ArrheniusRate high_rate(7.4e10, -0.37, 0.0);
    ArrheniusRate low_rate(2.3e12, -0.9, -7112800.0);
    vector<double> falloff_params { 0.7346, 94.0, 1756.0, 5182.0 };
    auto rate = make_shared<TroeRate>(low_rate, high_rate, falloff_params);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("AR:0.7 H2:2.0 H2O:6.0");
    auto R = make_shared<Reaction>(reac, prod, rate, tbody);

    kin.addReaction(R);
    check_rates(2);
}

TEST_F(KineticsFromScratch, add_falloff_reaction2)
{
    // reaction 2:
    //     equation: 2 OH (+ M) <=> H2O2 (+ M)  # Reaction 3
    //     duplicate: true
    //     type: falloff
    //     low-P-rate-constant: {A: 2.3e+12, b: -0.9, Ea: -1700.0}
    //     high-P-rate-constant: {A: 7.4e+10, b: -0.37, Ea: 0.0}
    //     Troe: {A: 0.7346, T3: 94.0, T1: 1756.0, T2: 5182.0}
    //     efficiencies: {AR: 0.7, H2: 2.0, H2O: 6.0}
    string equation = "2 OH (+ M) <=> H2O2 (+ M)";
    ArrheniusRate high_rate(7.4e10, -0.37, 0.0);
    ArrheniusRate low_rate(2.3e12, -0.9, -7112800.0);
    vector<double> falloff_params { 0.7346, 94.0, 1756.0, 5182.0 };
    auto rate = make_shared<TroeRate>(low_rate, high_rate, falloff_params);
    auto tbody = make_shared<ThirdBody>();
    tbody->efficiencies = parseCompString("AR:0.7 H2:2.0 H2O:6.0");
    auto R = make_shared<Reaction>(equation, rate, tbody);

    kin.addReaction(R);
    check_rates(2);
}

TEST_F(KineticsFromScratch, missing_third_body)
{
    string equation = "2 OH <=> H2O2";
    ArrheniusRate high_rate(7.4e10, -0.37, 0.0);
    ArrheniusRate low_rate(2.3e12, -0.9, -7112800.0);
    vector<double> falloff_params { 0.7346, 94.0, 1756.0, 5182.0 };
    auto rate = make_shared<TroeRate>(low_rate, high_rate, falloff_params);
    ASSERT_THROW(Reaction(equation, rate), CanteraError);
}

TEST_F(KineticsFromScratch, add_plog_reaction)
{
    // reaction 3:
    //     equation: H2 + O2 <=> 2 OH  # Reaction 4
    //     type: pressure-dependent-Arrhenius
    //     rate-constants:
    //     - {P: 0.01 atm, A: 1.2124e+16, b: -0.5779, Ea: 1.08727e+04}
    //     - {P: 1.0 atm, A: 4.9108e+31, b: -4.8507, Ea: 2.47728e+04}
    //     - {P: 10.0 atm, A: 1.2866e+47, b: -9.0246, Ea: 3.97965e+04}
    //     - {P: 100.0 atm, A: 5.9632e+56, b: -11.529, Ea: 5.25996e+04}
    Composition reac = parseCompString("H2:1, O2:1");
    Composition prod = parseCompString("OH:2");
    std::multimap<double, ArrheniusRate> rates {
        { 0.01*101325, ArrheniusRate(1.212400e+16, -0.5779, 10872.7 * 4184.0) },
        { 1.0*101325, ArrheniusRate(4.910800e+31, -4.8507, 24772.8 * 4184.0) },
        { 10.0*101325, ArrheniusRate(1.286600e+47, -9.0246, 39796.5 * 4184.0) },
        { 100.0*101325, ArrheniusRate(5.963200e+56, -11.529, 52599.6 * 4184.0) }
    };

    auto R = make_shared<Reaction>(reac, prod, make_shared<PlogRate>(rates));
    kin.addReaction(R);
    check_rates(3);
}

TEST_F(KineticsFromScratch, plog_invalid_rate)
{
    Composition reac = parseCompString("H2:1, O2:1");
    Composition prod = parseCompString("OH:2");
    std::multimap<double, ArrheniusRate> rates {
        { 0.01*101325, ArrheniusRate(1.2124e+16, -0.5779, 10872.7 * 4184.0) },
        { 10.0*101325, ArrheniusRate(1e15, -1, 10000 * 4184.0) },
        { 10.0*101325, ArrheniusRate(-2e20, -2.0, 20000 * 4184.0) },
        { 100.0*101325, ArrheniusRate(5.9632e+56, -11.529, 52599.6 * 4184.0) }
    };

    auto R = make_shared<Reaction>(reac, prod, make_shared<PlogRate>(rates));
    ASSERT_THROW(kin.addReaction(R), CanteraError);
}

TEST_F(KineticsFromScratch, add_chebyshev_reaction)
{
    // reaction 4:
    //     equation: HO2 <=> OH + O  # Reaction 5
    //     type: Chebyshev
    //     temperature-range: [290.0, 3000.0]
    //     pressure-range: [9.869232667160128e-03 atm, 98.69232667160128 atm]
    //     data:
    //     - [8.2883, -1.1397, -0.12059, 0.016034]
    //     - [1.9764, 1.0037, 7.2865e-03, -0.030432]
    //     - [0.3177, 0.26889, 0.094806, -7.6385e-03]
    Composition reac = parseCompString("HO2:1");
    Composition prod = parseCompString("OH:1 O:1");
    Array2D coeffs(3, 4);
    coeffs(0,0) = 8.2883e+00;
    coeffs(0,1) = -1.1397e+00;
    coeffs(0,2) = -1.2059e-01;
    coeffs(0,3) = 1.6034e-02;
    coeffs(1,0) = 1.9764e+00;
    coeffs(1,1) = 1.0037e+00;
    coeffs(1,2) = 7.2865e-03;
    coeffs(1,3) = -3.0432e-02;
    coeffs(2,0) = 3.1770e-01;
    coeffs(2,1) = 2.6889e-01;
    coeffs(2,2) = 9.4806e-02;
    coeffs(2,3) = -7.6385e-03;
    auto rate = make_shared<ChebyshevRate>(290., 3000., 1000.0, 10000000.0, coeffs);

    auto R = make_shared<Reaction>(reac, prod, rate);
    kin.addReaction(R);
    check_rates(4);
}

TEST_F(KineticsFromScratch, undeclared_species)
{
    Composition reac = parseCompString("CO:1 OH:1");
    Composition prod = parseCompString("CO2:1 H:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch, skip_undeclared_species)
{
    Composition reac = parseCompString("CO:1 OH:1");
    Composition prod = parseCompString("CO2:1 H:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);

    kin.skipUndeclaredSpecies(true);
    kin.addReaction(R);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch, negative_A_error)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(-3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    ASSERT_THROW(make_shared<Reaction>(reac, prod, rate), CanteraError);
}

TEST_F(KineticsFromScratch, allow_negative_A)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(-3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    rate->setAllowNegativePreExponentialFactor(true);
    auto R = make_shared<Reaction>(reac, prod, rate);

    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch, invalid_reversible_with_orders)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    R->orders["H2"] = 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch, negative_order_override)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    R->reversible = false;
    R->allow_negative_orders = true;
    R->orders["H2"] = - 0.5;

    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch, invalid_negative_orders)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    R->reversible = false;
    R->orders["H2"] = - 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch, nonreactant_order_override)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    R->reversible = false;
    R->allow_nonreactant_orders = true;
    R->orders["OH"] = 0.5;

    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch, invalid_nonreactant_order)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    R->reversible = false;
    R->orders["OH"] = 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch, unbalanced_neutral)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 H2O:1");
    auto rate = make_shared<ArrheniusRate>(1.0, 0.0, 0.0);
    auto R = make_shared<Reaction>(reac, prod, rate);
    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST(KineticsFromScratch2, unbalanced_positive_ion)
{
    auto soln = newSolution("kineticsfromscratch.yaml", "ions");
    Composition reac = parseCompString("H3O+:1");
    Composition prod = parseCompString("H2O:1 H:1");
    auto rate = make_shared<ArrheniusRate>(1.0, 0.0, 0.0);
    auto R = make_shared<Reaction>(reac, prod, rate);
    ASSERT_THROW(soln->kinetics()->addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, soln->kinetics()->nReactions());
}

TEST(KineticsFromScratch2, unbalanced_negative_ion)
{
    auto soln = newSolution("kineticsfromscratch.yaml", "ions");
    Composition reac = parseCompString("O2-:1");
    Composition prod = parseCompString("O2:1 E:2");
    auto rate = make_shared<ArrheniusRate>(1.0, 0.0, 0.0);
    auto R = make_shared<Reaction>(reac, prod, rate);
    ASSERT_THROW(soln->kinetics()->addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, soln->kinetics()->nReactions());
}

class InterfaceKineticsFromScratch : public testing::Test
{
public:
    InterfaceKineticsFromScratch()
        : gas(newThermo("sofc.yaml", "gas"))
        , gas_ref(newThermo("sofc.yaml", "gas"))
        , surf(newThermo("sofc.yaml", "metal_surface"))
        , surf_ref(newThermo("sofc.yaml", "metal_surface"))
    {
        kin_ref = newKinetics({surf_ref, gas_ref}, "sofc.yaml");
        kin.addThermo(surf);
        kin.addThermo(gas);
    }

    shared_ptr<ThermoPhase> gas;
    shared_ptr<ThermoPhase> gas_ref;
    shared_ptr<ThermoPhase> surf;
    shared_ptr<ThermoPhase> surf_ref;
    InterfaceKinetics kin;
    shared_ptr<Kinetics> kin_ref;

    //! iRef is the index of the corresponding reaction in the reference mech
    void check_rates(int iRef) {
        ASSERT_EQ((size_t) 1, kin.nReactions());

        string X = "H2:0.2 O2:0.5 H2O:0.1 N2:0.2";
        string Xs = "H(m):0.1 O(m):0.2 OH(m):0.3 (m):0.4";
        gas->setState_TPX(1200, 5*OneAtm, X);
        gas_ref->setState_TPX(1200, 5*OneAtm, X);
        surf->setState_TP(1200, 5*OneAtm);
        surf_ref->setState_TP(1200, 5*OneAtm);
        std::dynamic_pointer_cast<SurfPhase>(surf)->setCoveragesByName(Xs);
        std::dynamic_pointer_cast<SurfPhase>(surf_ref)->setCoveragesByName(Xs);

        vector<double> k(1), k_ref(kin_ref->nReactions());

        kin.getFwdRateConstants(&k[0]);
        kin_ref->getFwdRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);

        kin.getRevRateConstants(&k[0]);
        kin_ref->getRevRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);
    }
};

TEST_F(InterfaceKineticsFromScratch, add_surface_reaction)
{
    // Reaction 3 on the metal surface
    //     equation: H(m) + O(m) <=> OH(m) + (m)  # Reaction 4
    //     id: metal-rxn4
    //     rate-constant: {A: 5.0e+22, b: 0, Ea: 100.0}
    Composition reac = parseCompString("H(m):1 O(m):1");
    Composition prod = parseCompString("OH(m):1 (m):1");
    auto rate = make_shared<InterfaceArrheniusRate>(5e21, 0, 100.0e6);
    auto R = make_shared<Reaction>(reac, prod, rate);
    kin.addReaction(R);
    check_rates(3);
}

TEST_F(InterfaceKineticsFromScratch, add_sticking_reaction)
{
    // Reaction 0 on the metal surface
    //     equation: H2 + (m) + (m) <=> H(m) + H(m)  # Reaction 1
    //     sticking-coefficient: {A: 0.1, b: 0, Ea: 0}
    //     id: metal-rxn1
    Composition reac = parseCompString("H2:1 (m):2");
    Composition prod = parseCompString("H(m):2");
    auto rate = make_shared<StickingArrheniusRate>(0.1, 0, 0.0);
    auto R = make_shared<Reaction>(reac, prod, rate);
    kin.addReaction(R);
    check_rates(0);
}

TEST_F(InterfaceKineticsFromScratch, unbalanced_sites)
{
    Composition reac = parseCompString("H(m):1 O(m):1");
    Composition prod = parseCompString("OH(m):1");
    auto rate = make_shared<InterfaceArrheniusRate>(5e21, 0, 100.0e6);
    auto R = make_shared<Reaction>(reac, prod, rate);
    ASSERT_THROW(kin.addReaction(R), CanteraError);
}

class KineticsAddSpecies : public testing::Test
{
public:
    KineticsAddSpecies()
        : pp_ref(newThermo("../data/kineticsfromscratch.yaml"))
    {
        p = make_shared<IdealGasPhase>();
        vector<shared_ptr<ThermoPhase>> th;
        th.push_back(pp_ref);
        kin_ref = newKinetics(th, "../data/kineticsfromscratch.yaml");
        kin.addThermo(p);

        vector<shared_ptr<Species>> S = getSpecies(
            AnyMap::fromYamlFile("h2o2.yaml")["species"]);
        for (auto sp : S) {
            species[sp->name] = sp;
        }
        reactions = getReactions(
            AnyMap::fromYamlFile("../data/kineticsfromscratch.yaml")["reactions"],
            *kin_ref);
    }

    shared_ptr<ThermoPhase> p;
    shared_ptr<ThermoPhase> pp_ref;
    BulkKinetics kin;
    shared_ptr<Kinetics> kin_ref;
    vector<shared_ptr<Reaction>> reactions;
    map<string, shared_ptr<Species>> species;

    void check_rates(size_t N, const string& X) {
        for (size_t i = 0; i < kin_ref->nReactions(); i++) {
            if (i >= N) {
                kin_ref->setMultiplier(i, 0);
            } else {
                kin_ref->setMultiplier(i, 1);
            }
        }
        p->setState_TPX(1200, 5*OneAtm, X);
        pp_ref->setState_TPX(1200, 5*OneAtm, X);

        // need to invalidate cache to force update
        kin_ref->invalidateCache();

        vector<double> k(kin.nReactions()), k_ref(kin_ref->nReactions());
        vector<double> w(kin.nTotalSpecies()), w_ref(kin_ref->nTotalSpecies());

        kin.getFwdRateConstants(k.data());
        kin_ref->getFwdRateConstants(k_ref.data());
        for (size_t i = 0; i < kin.nReactions(); i++) {
            EXPECT_DOUBLE_EQ(k_ref[i], k[i]) << "i = " << i << "; N = " << N;
        }

        kin.getFwdRatesOfProgress(k.data());
        kin_ref->getFwdRatesOfProgress(k_ref.data());
        for (size_t i = 0; i < kin.nReactions(); i++) {
            EXPECT_DOUBLE_EQ(k_ref[i], k[i]) << "i = " << i << "; N = " << N;
        }

        kin.getRevRateConstants(k.data());
        kin_ref->getRevRateConstants(k_ref.data());
        for (size_t i = 0; i < kin.nReactions(); i++) {
            EXPECT_NEAR(k_ref[i], k[i], k_ref[i]*1e-12) << "i = " << i << "; N = " << N;
        }

        kin.getRevRatesOfProgress(k.data());
        kin_ref->getRevRatesOfProgress(k_ref.data());
        for (size_t i = 0; i < kin.nReactions(); i++) {
            EXPECT_NEAR(k_ref[i], k[i], k_ref[i]*1e-12) << "i = " << i << "; N = " << N;
        }

        kin.getCreationRates(w.data());
        kin_ref->getCreationRates(w_ref.data());
        for (size_t i = 0; i < kin.nTotalSpecies(); i++) {
            size_t iref = pp_ref->speciesIndex(p->speciesName(i));
            EXPECT_NEAR(w_ref[iref], w[i], w_ref[iref]*1e-12) << "sp = " << p->speciesName(i) << "; N = " << N;
        }
    }
};

TEST_F(KineticsAddSpecies, add_species_sequential)
{
    ASSERT_EQ((size_t) 0, kin.nReactions());

    for (auto s : {"AR", "O", "H2", "H", "OH"}) {
        p->addSpecies(species[s]);
    }
    kin.addReaction(reactions[0]);
    ASSERT_EQ(5, (int) kin.nTotalSpecies());
    check_rates(1, "O:0.001, H2:0.1, H:0.005, OH:0.02, AR:0.88");

    p->addSpecies(species["O2"]);
    p->addSpecies(species["H2O"]);
    kin.addReaction(reactions[1]);
    ASSERT_EQ(7, (int) kin.nTotalSpecies());
    ASSERT_EQ(2, (int) kin.nReactions());
    check_rates(2, "O:0.001, H2:0.1, H:0.005, OH:0.02, O2:0.5, AR:0.38");

    p->addSpecies(species["H2O2"]);
    kin.addReaction(reactions[2]);
    kin.addReaction(reactions[3]);
    check_rates(4, "O:0.001, H2:0.1, H:0.005, OH:0.02, O2:0.5, AR:0.38"); // no change
    check_rates(4, "O:0.001, H2:0.1, H:0.005, OH:0.02, O2:0.5, AR:0.35, H2O2:0.03");

    p->addSpecies(species["HO2"]);
    kin.addReaction(reactions[4]);
    check_rates(5, "O:0.01, H2:0.1, H:0.02, OH:0.03, O2:0.4, AR:0.3, H2O2:0.03, HO2:0.01");
}

TEST_F(KineticsAddSpecies, add_species_err_first)
{
    for (auto s : {"AR", "O", "H2", "H"}) {
        p->addSpecies(species[s]);
    }
    ASSERT_THROW(kin.addReaction(reactions[0]), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());

    p->addSpecies(species["OH"]);
    kin.addReaction(reactions[0]);
    ASSERT_EQ(5, (int) kin.nTotalSpecies());
    ASSERT_EQ((size_t) 1, kin.nReactions());
    check_rates(1, "O:0.001, H2:0.1, H:0.005, OH:0.02, AR:0.88");
}
