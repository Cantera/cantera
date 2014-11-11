#include "gtest/gtest.h"
#include "cantera/kinetics/importKinetics.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/kinetics/GasKinetics.h"

using namespace Cantera;

class KineticsFromScratch : public testing::Test
{
public:
    KineticsFromScratch()
        : p("../data/kineticsfromscratch.cti")
        , p_ref("../data/kineticsfromscratch.cti")
    {
        std::vector<ThermoPhase*> th;
        th.push_back(&p_ref);
        importKinetics(p_ref.xml(), th, &kin_ref);
        kin.addPhase(p);
        kin.init();
    }

    IdealGasPhase p;
    IdealGasPhase p_ref;
    GasKinetics kin;
    GasKinetics kin_ref;

    //! iRef is the index of the corresponding reaction in the reference mech
    void check_rates(int iRef) {
        ASSERT_EQ((size_t) 1, kin.nReactions());

        std::string X = "O:0.02 H2:0.2 O2:0.7 H:0.03 OH:0.05";
        p.setState_TPX(1200, 5*OneAtm, X);
        p_ref.setState_TPX(1200, 5*OneAtm, X);

        vector_fp k(1), k_ref(kin_ref.nReactions());

        kin.getFwdRateConstants(&k[0]);
        kin_ref.getFwdRateConstants(&k_ref[0]);
        EXPECT_FLOAT_EQ(k_ref[iRef], k[0]);

        kin.getRevRateConstants(&k[0]);
        kin_ref.getRevRateConstants(&k_ref[0]);
        EXPECT_FLOAT_EQ(k_ref[iRef], k[0]);
    }
};

TEST_F(KineticsFromScratch, add_elementary_reaction)
{
    // reaction 0:
    // reaction('O + H2 <=> H + OH', [3.870000e+01, 2.7, 6260.0])
    Composition reac = parseCompString("O:1 H2:1");
    Composition prod = parseCompString("H:1 OH:1");
    Arrhenius rate(3.87e1, 2.7, 6260.0 / GasConst_cal_mol_K);
    shared_ptr<ElementaryReaction> R(new ElementaryReaction(reac, prod, rate));

    kin.addReaction(R);
    kin.finalize();
    check_rates(0);
}

TEST_F(KineticsFromScratch, add_three_body_reaction)
{
    // reaction 1:
    // three_body_reaction('2 O + M <=> O2 + M', [1.200000e+11, -1.0, 0.0],
    //                     efficiencies='AR:0.83 H2:2.4 H2O:15.4')
    Composition reac = parseCompString("O:2");
    Composition prod = parseCompString("O2:1");
    Arrhenius rate(1.2e11, -1.0, 0.0);
    ThirdBody tbody;
    tbody.efficiencies = parseCompString("AR:0.83 H2:2.4 H2O:15.4");
    shared_ptr<ThirdBodyReaction> R(new ThirdBodyReaction(reac, prod, rate, tbody));

    kin.addReaction(R);
    kin.finalize();
    check_rates(1);
}

TEST_F(KineticsFromScratch, add_falloff_reaction)
{
    // reaction 2:
    // falloff_reaction('2 OH (+ M) <=> H2O2 (+ M)',
    //                  kf=[7.400000e+10, -0.37, 0.0],
    //                  kf0=[2.300000e+12, -0.9, -1700.0],
    //                  efficiencies='AR:0.7 H2:2.0 H2O:6.0',
    //                  falloff=Troe(A=0.7346, T3=94.0, T1=1756.0, T2=5182.0))
    Composition reac = parseCompString("OH:2");
    Composition prod = parseCompString("H2O2:1");
    Arrhenius high_rate(7.4e10, -0.37, 0.0);
    Arrhenius low_rate(2.3e12, -0.9, -1700.0 / GasConst_cal_mol_K);
    vector_fp falloff_params;
    falloff_params.push_back(0.7346);
    falloff_params.push_back(94.0);
    falloff_params.push_back(1756.0);
    falloff_params.push_back(5182.0);
    ThirdBody tbody;
    tbody.efficiencies = parseCompString("AR:0.7 H2:2.0 H2O:6.0");
    shared_ptr<FalloffReaction> R(new FalloffReaction(reac, prod, low_rate,
                                                      high_rate, tbody, TROE_FALLOFF,
                                                      falloff_params));
    kin.addReaction(R);
    kin.finalize();
    check_rates(2);
}

TEST_F(KineticsFromScratch, add_plog_reaction)
{
    // reaction 3:
    // pdep_arrhenius('H2 + O2 <=> 2 OH',
    //                [(0.01, 'atm'), 1.212400e+16, -0.5779, 10872.7],
    //                [(1.0, 'atm'), 4.910800e+31, -4.8507, 24772.8],
    //                [(10.0, 'atm'), 1.286600e+47, -9.0246, 39796.5],
    //                [(100.0, 'atm'), 5.963200e+56, -11.529, 52599.6])
    Composition reac = parseCompString("H2:1, O2:1");
    Composition prod = parseCompString("OH:2");
    std::multimap<double, Arrhenius> rates;
    typedef std::multimap<double, Arrhenius>::value_type item;
    rates.insert(item(0.01*101325, Arrhenius(1.212400e+16, -0.5779, 10872.7 / GasConst_cal_mol_K)));
    rates.insert(item(1.0*101325, Arrhenius(4.910800e+31, -4.8507, 24772.8 / GasConst_cal_mol_K)));
    rates.insert(item(10.0*101325, Arrhenius(1.286600e+47, -9.0246, 39796.5 / GasConst_cal_mol_K)));
    rates.insert(item(100.0*101325, Arrhenius(5.963200e+56, -11.529, 52599.6 / GasConst_cal_mol_K)));

    shared_ptr<PlogReaction> R(new PlogReaction(reac, prod, Plog(rates)));
    kin.addReaction(R);
    kin.finalize();
    check_rates(3);
}
