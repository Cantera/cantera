#include "gtest/gtest.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/thermo/Species.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/kinetics/InterfaceKinetics.h"
#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/base/Array.h"
#include "cantera/base/stringUtils.h"

using namespace Cantera;

class KineticsFromScratch3 : public testing::Test
{
public:
    KineticsFromScratch3()
        : p("../data/kineticsfromscratch.yaml")
        , p_ref("../data/kineticsfromscratch.yaml")
    {
        std::vector<ThermoPhase*> th;
        th.push_back(&p_ref);
        kin_ref = newKinetics(th, "../data/kineticsfromscratch.yaml", "ohmech");
        kin.addPhase(p);
        kin.init();
    }

    IdealGasPhase p;
    IdealGasPhase p_ref;
    GasKinetics kin;
    unique_ptr<Kinetics> kin_ref;

    //! iRef is the index of the corresponding reaction in the reference mech
    void check_rates(int iRef) {
        ASSERT_EQ((size_t) 1, kin.nReactions());

        std::string X = "O:0.02 H2:0.2 O2:0.5 H:0.03 OH:0.05 H2O:0.1 HO2:0.01";
        p.setState_TPX(1200, 5*OneAtm, X);
        p_ref.setState_TPX(1200, 5*OneAtm, X);

        vector_fp k(1), k_ref(kin_ref->nReactions());

        kin.getFwdRateConstants(&k[0]);
        kin_ref->getFwdRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);

        kin.getRevRateConstants(&k[0]);
        kin_ref->getRevRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);
    }
};

TEST_F(KineticsFromScratch3, add_elementary_reaction)
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

TEST_F(KineticsFromScratch3, add_three_body_reaction)
{
    // reaction 1:
    //     equation: 2 O + M <=> O2 + M  # Reaction 2
    //     type: three-body
    //     rate-constant: {A: 1.2e+11, b: -1.0, Ea: 0.0}
    //     efficiencies: {AR: 0.83, H2: 2.4, H2O: 15.4}
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    ArrheniusRate rate(1.2e11, -1.0, 0.0);
    ThirdBody tbody;
    tbody.efficiencies = parseCompString("AR:0.83 H2:2.4 H2O:15.4");
    auto R = make_shared<ThreeBodyReaction3>(reac, prod, rate, tbody);

    kin.addReaction(R);
    check_rates(1);
}

TEST_F(KineticsFromScratch3, undefined_third_body)
{
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    ArrheniusRate rate(1.2e11, -1.0, 0.0);
    ThirdBody tbody;
    tbody.efficiencies = parseCompString("H2:0.1 CO2:0.83");
    auto R = make_shared<ThreeBodyReaction3>(reac, prod, rate, tbody);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
}

TEST_F(KineticsFromScratch3, skip_undefined_third_body)
{
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    ArrheniusRate rate(1.2e11, -1.0, 0.0);
    ThirdBody tbody;
    tbody.efficiencies = parseCompString("H2:0.1 CO2:0.83");
    auto R = make_shared<ThreeBodyReaction3>(reac, prod, rate, tbody);

    kin.skipUndeclaredThirdBodies(true);
    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch3, add_falloff_reaction)
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
    vector_fp falloff_params { 0.7346, 94.0, 1756.0, 5182.0 };
    TroeRate rate(low_rate, high_rate, falloff_params);
    ThirdBody tbody;
    tbody.efficiencies = parseCompString("AR:0.7 H2:2.0 H2O:6.0");
    auto R = make_shared<FalloffReaction3>(reac, prod, rate, tbody);
    kin.addReaction(R);
    check_rates(2);
}

TEST_F(KineticsFromScratch3, add_plog_reaction)
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

TEST_F(KineticsFromScratch3, plog_invalid_rate)
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

TEST_F(KineticsFromScratch3, add_chebyshev_reaction)
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

TEST_F(KineticsFromScratch3, undeclared_species)
{
    Composition reac = parseCompString("CO:1 OH:1");
    Composition prod = parseCompString("CO2:1 H:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, skip_undeclared_species)
{
    Composition reac = parseCompString("CO:1 OH:1");
    Composition prod = parseCompString("CO2:1 H:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);

    kin.skipUndeclaredSpecies(true);
    kin.addReaction(R);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, negative_A_error)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(-3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, allow_negative_A)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(-3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    auto rr = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    rr->setAllowNegativePreExponentialFactor(true);

    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch3, invalid_reversible_with_orders)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    auto rate = make_shared<ArrheniusRate>(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<Reaction>(reac, prod, rate);
    R->orders["H2"] = 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, negative_order_override)
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

TEST_F(KineticsFromScratch3, invalid_negative_orders)
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

TEST_F(KineticsFromScratch3, nonreactant_order_override)
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

TEST_F(KineticsFromScratch3, invalid_nonreactant_order)
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

class InterfaceKineticsFromScratch3 : public testing::Test
{
public:
    InterfaceKineticsFromScratch3()
        : gas("sofc.yaml", "gas")
        , gas_ref("sofc.yaml", "gas")
        , surf("sofc.yaml", "metal_surface")
        , surf_ref("sofc.yaml", "metal_surface")
    {
        kin_ref = newKinetics({&surf_ref, &gas_ref}, "sofc.yaml", "metal_surface");
        kin.addPhase(surf);
        kin.addPhase(gas);
    }

    IdealGasPhase gas;
    IdealGasPhase gas_ref;
    SurfPhase surf;
    SurfPhase surf_ref;
    InterfaceKinetics kin;
    unique_ptr<Kinetics> kin_ref;

    //! iRef is the index of the corresponding reaction in the reference mech
    void check_rates(int iRef) {
        ASSERT_EQ((size_t) 1, kin.nReactions());

        std::string X = "H2:0.2 O2:0.5 H2O:0.1 N2:0.2";
        std::string Xs = "H(m):0.1 O(m):0.2 OH(m):0.3 (m):0.4";
        gas.setState_TPX(1200, 5*OneAtm, X);
        gas_ref.setState_TPX(1200, 5*OneAtm, X);
        surf.setState_TP(1200, 5*OneAtm);
        surf_ref.setState_TP(1200, 5*OneAtm);
        surf.setCoveragesByName(Xs);
        surf_ref.setCoveragesByName(Xs);

        vector_fp k(1), k_ref(kin_ref->nReactions());

        kin.getFwdRateConstants(&k[0]);
        kin_ref->getFwdRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);

        kin.getRevRateConstants(&k[0]);
        kin_ref->getRevRateConstants(&k_ref[0]);
        EXPECT_DOUBLE_EQ(k_ref[iRef], k[0]);
    }
};

TEST_F(InterfaceKineticsFromScratch3, add_surface_reaction)
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

TEST_F(InterfaceKineticsFromScratch3, add_sticking_reaction)
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

TEST_F(InterfaceKineticsFromScratch3, unbalanced_sites)
{
    Composition reac = parseCompString("H(m):1 O(m):1");
    Composition prod = parseCompString("OH(m):1");
    auto rate = make_shared<InterfaceArrheniusRate>(5e21, 0, 100.0e6);
    auto R = make_shared<Reaction>(reac, prod, rate);
    ASSERT_THROW(kin.addReaction(R), CanteraError);
}

class KineticsAddSpecies3 : public testing::Test
{
public:
    KineticsAddSpecies3()
        : p_ref("../data/kineticsfromscratch.yaml")
    {
        std::vector<ThermoPhase*> th;
        th.push_back(&p_ref);
        kin_ref = newKinetics(th, "../data/kineticsfromscratch.yaml", "ohmech");
        kin.addPhase(p);

        std::vector<shared_ptr<Species>> S = getSpecies(
            AnyMap::fromYamlFile("h2o2.yaml")["species"]);
        for (auto sp : S) {
            species[sp->name] = sp;
        }
        reactions = getReactions(
            AnyMap::fromYamlFile("../data/kineticsfromscratch.yaml")["reactions"],
            *kin_ref);
    }

    IdealGasPhase p;
    IdealGasPhase p_ref;
    GasKinetics kin;
    unique_ptr<Kinetics> kin_ref;
    std::vector<shared_ptr<Reaction>> reactions;
    std::map<std::string, shared_ptr<Species>> species;

    void check_rates(size_t N, const std::string& X) {
        for (size_t i = 0; i < kin_ref->nReactions(); i++) {
            if (i >= N) {
                kin_ref->setMultiplier(i, 0);
            } else {
                kin_ref->setMultiplier(i, 1);
            }
        }
        p.setState_TPX(1200, 5*OneAtm, X);
        p_ref.setState_TPX(1200, 5*OneAtm, X);

        // need to invalidate cache to force update
        kin_ref->invalidateCache();

        vector_fp k(kin.nReactions()), k_ref(kin_ref->nReactions());
        vector_fp w(kin.nTotalSpecies()), w_ref(kin_ref->nTotalSpecies());

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
            size_t iref = p_ref.speciesIndex(p.speciesName(i));
            EXPECT_NEAR(w_ref[iref], w[i], w_ref[iref]*1e-12) << "sp = " << p.speciesName(i) << "; N = " << N;
        }
    }
};

TEST_F(KineticsAddSpecies3, add_species_sequential)
{
    ASSERT_EQ((size_t) 0, kin.nReactions());

    for (auto s : {"AR", "O", "H2", "H", "OH"}) {
        p.addSpecies(species[s]);
    }
    kin.addReaction(reactions[0]);
    ASSERT_EQ(5, (int) kin.nTotalSpecies());
    check_rates(1, "O:0.001, H2:0.1, H:0.005, OH:0.02, AR:0.88");

    p.addSpecies(species["O2"]);
    p.addSpecies(species["H2O"]);
    kin.addReaction(reactions[1]);
    ASSERT_EQ(7, (int) kin.nTotalSpecies());
    ASSERT_EQ(2, (int) kin.nReactions());
    check_rates(2, "O:0.001, H2:0.1, H:0.005, OH:0.02, O2:0.5, AR:0.38");

    p.addSpecies(species["H2O2"]);
    kin.addReaction(reactions[2]);
    kin.addReaction(reactions[3]);
    check_rates(4, "O:0.001, H2:0.1, H:0.005, OH:0.02, O2:0.5, AR:0.38"); // no change
    check_rates(4, "O:0.001, H2:0.1, H:0.005, OH:0.02, O2:0.5, AR:0.35, H2O2:0.03");

    p.addSpecies(species["HO2"]);
    kin.addReaction(reactions[4]);
    check_rates(5, "O:0.01, H2:0.1, H:0.02, OH:0.03, O2:0.4, AR:0.3, H2O2:0.03, HO2:0.01");
}

TEST_F(KineticsAddSpecies3, add_species_err_first)
{
    for (auto s : {"AR", "O", "H2", "H"}) {
        p.addSpecies(species[s]);
    }
    ASSERT_THROW(kin.addReaction(reactions[0]), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());

    p.addSpecies(species["OH"]);
    kin.addReaction(reactions[0]);
    ASSERT_EQ(5, (int) kin.nTotalSpecies());
    ASSERT_EQ((size_t) 1, kin.nReactions());
    check_rates(1, "O:0.001, H2:0.1, H:0.005, OH:0.02, AR:0.88");
}
