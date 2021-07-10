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
    // reaction('O + H2 <=> H + OH', [3.870000e+01, 2.7, 6260.0])
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    ArrheniusRate rate(3.87e1, 2.7, 2.619184e+07);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);

    kin.addReaction(R);
    check_rates(0);
}

TEST_F(KineticsFromScratch3, add_three_body_reaction)
{
    // reaction 1:
    // three_body_reaction('2 O + M <=> O2 + M', [1.200000e+11, -1.0, 0.0],
    //                     efficiencies='AR:0.83 H2:2.4 H2O:15.4')
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

TEST_F(KineticsFromScratch3, add_plog_reaction)
{
    // reaction 3:
    // pdep_arrhenius('H2 + O2 <=> 2 OH',
    //                [(0.01, 'atm'), 1.212400e+16, -0.5779, 10872.7],
    //                [(1.0, 'atm'), 4.910800e+31, -4.8507, 24772.8],
    //                [(10.0, 'atm'), 1.286600e+47, -9.0246, 39796.5],
    //                [(100.0, 'atm'), 5.963200e+56, -11.529, 52599.6])
    Composition reac = parseCompString("H2:1, O2:1");
    Composition prod = parseCompString("OH:2");
    std::vector<std::pair<double, Arrhenius>> rates {
        { 0.01*101325, Arrhenius(1.212400e+16, -0.5779, 10872.7 / GasConst_cal_mol_K) },
        { 1.0*101325, Arrhenius(4.910800e+31, -4.8507, 24772.8 / GasConst_cal_mol_K) },
        { 10.0*101325, Arrhenius(1.286600e+47, -9.0246, 39796.5 / GasConst_cal_mol_K) },
        { 100.0*101325, Arrhenius(5.963200e+56, -11.529, 52599.6 / GasConst_cal_mol_K) }
    };

    auto R = make_shared<PlogReaction3>(reac, prod, PlogRate(rates));
    kin.addReaction(R);
    check_rates(3);
}

TEST_F(KineticsFromScratch3, plog_invalid_rate)
{
    Composition reac = parseCompString("H2:1, O2:1");
    Composition prod = parseCompString("OH:2");
    std::vector<std::pair<double, Arrhenius>> rates {
        { 0.01*101325, Arrhenius(1.2124e+16, -0.5779, 10872.7 / GasConst_cal_mol_K) },
        { 10.0*101325, Arrhenius(1e15, -1, 10000 / GasConst_cal_mol_K) },
        { 10.0*101325, Arrhenius(-2e20, -2.0, 20000 / GasConst_cal_mol_K) },
        { 100.0*101325, Arrhenius(5.9632e+56, -11.529, 52599.6 / GasConst_cal_mol_K) }
    };

    auto R = make_shared<PlogReaction3>(reac, prod, PlogRate(rates));
    ASSERT_THROW(kin.addReaction(R), CanteraError);
}

TEST_F(KineticsFromScratch3, add_chebyshev_reaction)
{
    // reaction 4:
    // chebyshev_reaction(
    //     'HO2 <=> OH + O',
    //     Tmin=290.0, Tmax=3000.0,
    //     Pmin=(0.0098692326671601278, 'atm'), Pmax=(98.692326671601279, 'atm'),
    //     coeffs=[[ 8.2883e+00, -1.1397e+00, -1.2059e-01,  1.6034e-02],
    //             [ 1.9764e+00,  1.0037e+00,  7.2865e-03, -3.0432e-02],
    //             [ 3.1770e-01,  2.6889e-01,  9.4806e-02, -7.6385e-03]])
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
    ChebyshevRate3 rate(290, 3000, 1000.0, 10000000.0, coeffs);

    auto R = make_shared<ChebyshevReaction3>(reac, prod, rate);
    kin.addReaction(R);
    check_rates(4);
}

TEST_F(KineticsFromScratch3, undeclared_species)
{
    Composition reac = parseCompString("CO:1 OH:1");
    Composition prod = parseCompString("CO2:1 H:1");
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, skip_undeclared_species)
{
    Composition reac = parseCompString("CO:1 OH:1");
    Composition prod = parseCompString("CO2:1 H:1");
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);

    kin.skipUndeclaredSpecies(true);
    kin.addReaction(R);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, negative_A_error)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    ArrheniusRate rate(-3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, allow_negative_A)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    ArrheniusRate rate(-3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);
    auto rr = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    rr->allow_negative_pre_exponential_factor = true;

    kin.addReaction(R);
    ASSERT_EQ((size_t) 1, kin.nReactions());
}

TEST_F(KineticsFromScratch3, invalid_reversible_with_orders)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);
    R->orders["H2"] = 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, negative_order_override)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);
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
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);
    R->reversible = false;
    R->orders["H2"] = - 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
}

TEST_F(KineticsFromScratch3, nonreactant_order_override)
{
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);
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
    ArrheniusRate rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    auto R = make_shared<ElementaryReaction3>(reac, prod, rate);
    R->reversible = false;
    R->orders["OH"] = 0.5;

    ASSERT_THROW(kin.addReaction(R), CanteraError);
    ASSERT_EQ((size_t) 0, kin.nReactions());
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
            EXPECT_DOUBLE_EQ(k_ref[i], k[i]) << "i = " << i << "; N = " << N;
        }

        kin.getRevRatesOfProgress(k.data());
        kin_ref->getRevRatesOfProgress(k_ref.data());
        for (size_t i = 0; i < kin.nReactions(); i++) {
            EXPECT_DOUBLE_EQ(k_ref[i], k[i]) << "i = " << i << "; N = " << N;
        }

        kin.getCreationRates(w.data());
        kin_ref->getCreationRates(w_ref.data());
        for (size_t i = 0; i < kin.nTotalSpecies(); i++) {
            size_t iref = p_ref.speciesIndex(p.speciesName(i));
            EXPECT_DOUBLE_EQ(w_ref[iref], w[i]) << "sp = " << p.speciesName(i) << "; N = " << N;
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
