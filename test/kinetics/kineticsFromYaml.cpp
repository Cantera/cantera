#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/base/Solution.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/kinetics/ReactionFactory.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/Array.h"

using namespace Cantera;

TEST(ReactionRate, ModifyArrheniusRate)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    auto ER = dynamic_cast<ElementaryReaction3&>(*R);

    auto rr = std::dynamic_pointer_cast<ArrheniusRate>(ER.rate());
    EXPECT_TRUE(rr->allow_negative_pre_exponential_factor);
    rr->allow_negative_pre_exponential_factor = false;
    EXPECT_FALSE(rr->allow_negative_pre_exponential_factor);
}

TEST(Reaction, ElementaryFromYaml3)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->type(), "elementary");
    EXPECT_FALSE(R->allow_negative_orders);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_TRUE(rate->allow_negative_pre_exponential_factor);
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), -2.7e10);
    EXPECT_DOUBLE_EQ(rate->activationEnergy_R(), 355 / GasConst_cal_mol_K);
}

TEST(Reaction, ElementaryFromYaml2)
{
    auto sol = newSolution("gri30.yaml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " type: elementary-legacy,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->type(), "elementary-legacy");

    auto ER = dynamic_cast<ElementaryReaction2&>(*R);
    EXPECT_DOUBLE_EQ(ER.rate.preExponentialFactor(), -2.7e10);
    EXPECT_DOUBLE_EQ(ER.rate.activationEnergy_R(), 355 / GasConst_cal_mol_K);
    EXPECT_TRUE(ER.allow_negative_pre_exponential_factor);
    EXPECT_FALSE(ER.allow_negative_orders);
}

TEST(Reaction, ThreeBodyFromYaml3)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O + M = O2 + M,"
        " type: three-body,"
        " rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
        " efficiencies: {AR: 0.83, H2O: 5}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);
    EXPECT_EQ(R->type(), "three-body");
    EXPECT_DOUBLE_EQ(R->thirdBody()->efficiencies["H2O"], 5.0);
    EXPECT_DOUBLE_EQ(R->thirdBody()->default_efficiency, 1.0);

    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(R->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 1.2e11);
}

TEST(Reaction, ThreeBodyFromYamlMissingM)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O <=> O2," // Missing "M" on each side of the equation
        " type: three-body,"
        " rate-constant: [1.20000E+17, -1, 0]}");

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, ThreeBodyFromYaml2)
{
    auto sol = newSolution("gri30.yaml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O + M = O2 + M,"
        " type: three-body-legacy,"
        " rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
        " efficiencies: {AR: 0.83, H2O: 5}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);
    EXPECT_EQ(R->type(), "three-body-legacy");

    auto TBR = dynamic_cast<ThreeBodyReaction2&>(*R);
    EXPECT_DOUBLE_EQ(TBR.rate.preExponentialFactor(), 1.2e11);
    EXPECT_DOUBLE_EQ(TBR.third_body.efficiencies["H2O"], 5.0);
    EXPECT_DOUBLE_EQ(TBR.third_body.default_efficiency, 1.0);
}

TEST(Reaction, FalloffFromYaml1)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N2O (+M) = N2 + O (+ M),"
        " type: falloff,"
        " high-P-rate-constant: [7.91000E+10, 0, 56020],"
        " low-P-rate-constant: [6.37000E+14, 0, 56640],"
        " SRI: {A: 1.1, B: 700.0, C: 1234.0, D: 56.0, E: 0.7},"
        " efficiencies: {AR: 0.625}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    auto FR = dynamic_cast<FalloffReaction&>(*R);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("AR"), 0.625);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("N2"), 1.0);
}

TEST(Reaction, FalloffFromYaml2)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: H + CH2 (+ N2) <=> CH3 (+N2),"
        " type: falloff,"
        " high-P-rate-constant: [6.00000E+14 cm^3/mol/s, 0, 0],"
        " low-P-rate-constant: {A: 1.04000E+26 cm^6/mol^2/s, b: -2.76, Ea: 1600},"
        " Troe: {A: 0.562, T3: 91, T1: 5836},"
        " source: somewhere}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    auto FR = dynamic_cast<FalloffReaction&>(*R);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("N2"), 1.0);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("H2O"), 0.0);
    EXPECT_DOUBLE_EQ(FR.high_rate.preExponentialFactor(), 6e11);
    EXPECT_DOUBLE_EQ(FR.low_rate.preExponentialFactor(), 1.04e20);
    EXPECT_DOUBLE_EQ(FR.low_rate.activationEnergy_R(), 1600 / GasConstant);
    vector_fp params(4);
    FR.falloff->getParameters(params.data());
    EXPECT_DOUBLE_EQ(params[0], 0.562);
    EXPECT_DOUBLE_EQ(params[1], 91.0);
    EXPECT_DOUBLE_EQ(params[3], 0.0);
    EXPECT_EQ(R->input["source"].asString(), "somewhere");
}

TEST(Reaction, FalloffFromYaml3)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: HCN (+ M) <=> H + CN (+ M),"
        " type: falloff,"
        " low-P-rate-constant: {A: 3.57e+26, b: -2.6, Ea: 1.249e+05},"
        " high-P-rate-constant: {A: 8.3e+17, b: -0.93, Ea: 1.238e+05},"
        " Tsang: {A: 0.95, B: -1.0e-04},"
        " efficiencies: {CO2: 1.6, H2O: 5.0, N2: 1.0, N2O: 5.0},"
        " source: ARL-TR-5088}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    auto FR = dynamic_cast<FalloffReaction&>(*R);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("N2"), 1.0);
    EXPECT_DOUBLE_EQ(FR.third_body.efficiency("H2O"), 5.0);
    EXPECT_DOUBLE_EQ(FR.high_rate.preExponentialFactor(), 8.3e17);
    EXPECT_DOUBLE_EQ(FR.low_rate.preExponentialFactor(), 3.57e26);
    EXPECT_DOUBLE_EQ(FR.high_rate.activationEnergy_R(), 123800.0 / GasConstant);
    EXPECT_DOUBLE_EQ(FR.low_rate.activationEnergy_R(), 124900.0 / GasConstant);
    vector_fp params(2);
    FR.falloff->getParameters(params.data());
    EXPECT_DOUBLE_EQ(params[0], 0.95);
    EXPECT_DOUBLE_EQ(params[1], -0.0001);
    EXPECT_EQ(R->input["source"].asString(), "ARL-TR-5088");
}

TEST(Reaction, ChemicallyActivatedFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: CH3 + OH (+M) <=> CH2O + H2 (+M),"
        " units: {length: cm, quantity: mol},"
        " type: chemically-activated,"
        " high-P-rate-constant: [5.88E-14, 6.721, -3022.227],"
        " low-P-rate-constant: [282320.078, 1.46878, -3270.56495]}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    auto CAR = dynamic_cast<ChemicallyActivatedReaction&>(*R);
    EXPECT_DOUBLE_EQ(CAR.high_rate.preExponentialFactor(), 5.88e-14);
    EXPECT_DOUBLE_EQ(CAR.low_rate.preExponentialFactor(), 2.82320078e2);
    EXPECT_EQ(CAR.falloff->nParameters(), (size_t) 0);
}

TEST(Reaction, PlogFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "equation: 'H + CH4 <=> H2 + CH3'\n"
        "units: {pressure: atm}\n"
        "type: pressure-dependent-Arrhenius\n"
        "rate-constants:\n"
        "- {P: 0.039474, A: 2.720000e+09 cm^3/mol/s, b: 1.2, Ea: 6834.0}\n"
        "- {P: 1.0 atm, A: 1.260000e+20, b: -1.83, Ea: 15003.0}\n"
        "- {P: 1.0 atm, A: 1.230000e+04, b: 2.68, Ea: 6335.0}\n"
        "- {P: 1.01325 MPa, A: 1.680000e+16, b: -0.6, Ea: 14754.0}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    const auto& rateMap = std::dynamic_pointer_cast<PlogRate>(R->rate())->getRates();
    std::vector<std::pair<double, Arrhenius>> rates(rateMap.begin(), rateMap.end());
    EXPECT_EQ(rates.size(), (size_t) 4);
    EXPECT_NEAR(rates[0].first, 0.039474 * OneAtm, 1e-6);
    EXPECT_NEAR(rates[2].first, OneAtm, 1e-6);
    EXPECT_NEAR(rates[3].first, 10 * OneAtm, 1e-6);
    EXPECT_DOUBLE_EQ(rates[0].second.preExponentialFactor(), 2.72e6);
    EXPECT_DOUBLE_EQ(rates[3].second.preExponentialFactor(), 1.68e16);
    EXPECT_DOUBLE_EQ(rates[3].second.temperatureExponent(), -0.6);
}

TEST(Reaction, ChebyshevFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "equation: 'CH4 <=> CH3 + H'\n"
        "type: Chebyshev\n"
        "temperature-range: [290, 3000]\n"
        "pressure-range: [0.0098692326671601278 atm, 98.692326671601279 atm]\n"
        "data: [[-1.44280e+01,  2.59970e-01, -2.24320e-02, -2.78700e-03],\n"
        "       [ 2.20630e+01,  4.88090e-01, -3.96430e-02, -5.48110e-03],\n"
        "       [-2.32940e-01,  4.01900e-01, -2.60730e-02, -5.04860e-03],\n"
        "       [-2.93660e-01,  2.85680e-01, -9.33730e-03, -4.01020e-03],\n"
        "       [-2.26210e-01,  1.69190e-01,  4.85810e-03, -2.38030e-03],\n"
        "       [-1.43220e-01,  7.71110e-02,  1.27080e-02, -6.41540e-04]]\n");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.size(), (size_t) 1);
    const auto rate = std::dynamic_pointer_cast<ChebyshevRate3>(R->rate());
    double logP = std::log10(2e6);
    double T = 1800;
    rate->update_C(&logP);
    EXPECT_EQ(rate->data().nRows(), (size_t) 6);
    EXPECT_EQ(rate->data().nColumns(), (size_t) 4);
    EXPECT_DOUBLE_EQ(rate->Tmax(), 3000);
    EXPECT_DOUBLE_EQ(rate->Pmin(), 1000);
    EXPECT_NEAR(rate->updateRC(std::log(T), 1.0/T), 130512.2773948636, 1e-9);
}

TEST(Reaction, BlowersMaselFromYaml)
{
    auto sol = newSolution("gri30.yaml", "", "None");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: O + H2 <=> H + OH,"
        " type: Blowers-Masel,"
        " rate-constant: [-3.87e+04 cm^3/mol/s, 2.7, 6260.0 cal/mol, 1e9 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.at("H2"), 1);
    EXPECT_EQ(R->products.at("OH"), 1);
    EXPECT_EQ(R->reaction_type, BLOWERSMASEL_RXN);

    auto ER = dynamic_cast<BlowersMaselReaction&>(*R);
    doublereal E_intrinsic = 6260 / GasConst_cal_mol_K * GasConstant; // J/kmol
    doublereal H_big = 5 * E_intrinsic;
    doublereal H_small = -5 * E_intrinsic;
    doublereal H_mid = 4 * E_intrinsic;
    doublereal w = 1e9 / GasConst_cal_mol_K * GasConstant; // J/kmol
    doublereal vp = 2 * w * ((w + E_intrinsic) / (w - E_intrinsic));
    doublereal Ea = (w + H_mid / 2) * (vp - 2 * w + H_mid) * (vp - 2 * w + H_mid)
                    / (vp * vp - 4 * w * w + H_mid * H_mid );
    EXPECT_DOUBLE_EQ(ER.rate.preExponentialFactor(), -38.7);
    EXPECT_DOUBLE_EQ(ER.rate.activationEnergy_R0(), 6260 / GasConst_cal_mol_K);
    EXPECT_DOUBLE_EQ(ER.rate.bondEnergy(), 1e9 / GasConst_cal_mol_K);
    EXPECT_DOUBLE_EQ(ER.rate.activationEnergy_R(H_big), H_big / GasConstant);
    EXPECT_DOUBLE_EQ(ER.rate.activationEnergy_R(H_small), 0);
    EXPECT_NEAR(ER.rate.activationEnergy_R(H_mid), Ea / GasConstant, 1e-7);
    EXPECT_TRUE(ER.allow_negative_pre_exponential_factor);
    EXPECT_FALSE(ER.allow_negative_orders);
}

TEST(Kinetics, GasKineticsFromYaml1)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode = infile["phases"].getMapWhere("name", "simple-kinetics");
    shared_ptr<ThermoPhase> thermo = newPhase(phaseNode, infile);
    auto kin = newKinetics({thermo.get()}, phaseNode, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 2);
    const auto& R = kin->reaction(0);
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->id, "NOx-R1");
    const auto& ER = std::dynamic_pointer_cast<ElementaryReaction3>(R);
    const auto& rate = std::dynamic_pointer_cast<ArrheniusRate>(ER->rate());
    EXPECT_DOUBLE_EQ(rate->preExponentialFactor(), 2.7e10);
}

TEST(Kinetics, GasKineticsFromYaml2)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode = infile["phases"].getMapWhere("name", "remote-kinetics");
    shared_ptr<ThermoPhase> thermo = newPhase(phaseNode, infile);
    auto kin = newKinetics({thermo.get()}, phaseNode, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 3);
}

TEST(Kinetics, EfficienciesFromYaml)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode1 = infile["phases"].getMapWhere("name", "efficiency-error");
    shared_ptr<ThermoPhase> thermo1 = newPhase(phaseNode1, infile);
    // Reaction with efficiency defined for undeclared species "AR"
    EXPECT_THROW(newKinetics({thermo1.get()}, phaseNode1, infile), CanteraError);

    auto& phaseNode2 = infile["phases"].getMapWhere("name", "efficiency-skip");
    shared_ptr<ThermoPhase> thermo2 = newPhase(phaseNode2, infile);
    auto kin = newKinetics({thermo2.get()}, phaseNode2, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 1);
}

TEST(Kinetics, InterfaceKineticsFromYaml)
{
    shared_ptr<ThermoPhase> gas(newPhase("surface-phases.yaml", "gas"));
    shared_ptr<ThermoPhase> surf_tp(newPhase("surface-phases.yaml", "Pt-surf"));
    shared_ptr<SurfPhase> surf = std::dynamic_pointer_cast<SurfPhase>(surf_tp);
    auto kin = newKinetics({surf_tp.get(), gas.get()},
                           "surface-phases.yaml", "Pt-surf");
    EXPECT_EQ(kin->nReactions(), (size_t) 3);
    EXPECT_EQ(kin->nTotalSpecies(), (size_t) 6);
    auto R1 = kin->reaction(0);
    auto IR1 = std::dynamic_pointer_cast<InterfaceReaction>(R1);
    EXPECT_DOUBLE_EQ(R1->orders["Pt(s)"], 1.0);
    EXPECT_DOUBLE_EQ(IR1->rate.preExponentialFactor(), 4.4579e7);

    auto R2 = kin->reaction(1);
    auto IR2 = std::dynamic_pointer_cast<InterfaceReaction>(R2);
    EXPECT_DOUBLE_EQ(IR2->rate.preExponentialFactor(), 3.7e20);
    EXPECT_DOUBLE_EQ(IR2->coverage_deps["H(s)"].E, -6e6 / GasConstant);
    EXPECT_FALSE(IR2->is_sticking_coefficient);

    auto R3 = kin->reaction(2);
    auto IR3 = std::dynamic_pointer_cast<InterfaceReaction>(R3);
    EXPECT_TRUE(IR3->is_sticking_coefficient);
}

TEST(Kinetics, BMInterfaceKineticsFromYaml)
{
    shared_ptr<ThermoPhase> gas(newPhase("BM_test.yaml", "gas"));
    shared_ptr<ThermoPhase> surf_tp(newPhase("BM_test.yaml", "Pt_surf"));
    shared_ptr<SurfPhase> surf = std::dynamic_pointer_cast<SurfPhase>(surf_tp);
    auto kin = newKinetics({surf_tp.get(), gas.get()}, "BM_test.yaml", "Pt_surf");
    EXPECT_EQ(kin->nReactions(), (size_t) 6);
    EXPECT_EQ(kin->nTotalSpecies(), (size_t) 14);
    auto R1 = kin->reaction(5);
    auto IR1 = std::dynamic_pointer_cast<BlowersMaselInterfaceReaction>(R1);
    EXPECT_DOUBLE_EQ(R1->orders["PT(s)"], 1.0);
    EXPECT_DOUBLE_EQ(IR1->rate.preExponentialFactor(), 4.4579e7);

    auto R2 = kin->reaction(0);
    auto IR2 = std::dynamic_pointer_cast<BlowersMaselInterfaceReaction>(R2);
    EXPECT_DOUBLE_EQ(IR2->rate.preExponentialFactor(), 3.7e20);
    EXPECT_DOUBLE_EQ(IR2->coverage_deps["H(S)"].E, -6e6 / GasConstant);
    EXPECT_FALSE(IR2->is_sticking_coefficient);

    auto R3 = kin->reaction(1);
    auto IR3 = std::dynamic_pointer_cast<BlowersMaselInterfaceReaction>(R3);
    EXPECT_TRUE(IR3->is_sticking_coefficient);
}

TEST(Kinetics, ElectrochemFromYaml)
{
    shared_ptr<ThermoPhase> graphite(newPhase("surface-phases.yaml", "graphite"));
    shared_ptr<ThermoPhase> electrolyte(newPhase("surface-phases.yaml", "electrolyte"));
    shared_ptr<ThermoPhase> anode(newPhase("surface-phases.yaml", "anode-surface"));
    auto kin = newKinetics({anode.get(), graphite.get(), electrolyte.get()},
                           "surface-phases.yaml", "anode-surface");
    graphite->setElectricPotential(0.4);
    vector_fp ropf(kin->nReactions()), ropr(kin->nReactions());
    kin->getFwdRatesOfProgress(ropf.data());
    kin->getRevRatesOfProgress(ropr.data());

    EXPECT_NEAR(ropf[0], 0.279762338, 1e-8);
    EXPECT_NEAR(ropr[0], 0.045559670, 1e-8);
}

TEST(KineticsFromYaml, NoKineticsModelOrReactionsField1)
{
    auto soln = newSolution("phase-reaction-spec1.yaml",
                            "nokinetics-noreactions");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "Kinetics");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, NoKineticsModelOrReactionsField2)
{
    auto soln = newSolution("phase-reaction-spec2.yaml",
                            "nokinetics-noreactions");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "Kinetics");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithReactionsNone1)
{
    auto soln = newSolution("phase-reaction-spec1.yaml",
                            "kinetics-reactions-none");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "Gas");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithReactionsNone2)
{
    auto soln = newSolution("phase-reaction-spec2.yaml",
                            "kinetics-reactions-none");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "Gas");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 0);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsSection1)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "kinetics-no-reaction-section1"),
                 InputFileError);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsSection2)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "kinetics-no-reaction-section2"),
                 InputFileError);
}

TEST(KineticsFromYaml, KineticsModelWithoutReactionsField)
{
    auto soln = newSolution("phase-reaction-spec2.yaml",
                            "kinetics-noreactions");
    EXPECT_EQ(soln->kinetics()->kineticsType(), "Gas");
    EXPECT_EQ(soln->kinetics()->nReactions(), (size_t) 1);
}

TEST(KineticsFromYaml, ReactionsFieldWithoutKineticsModel1)
{
    EXPECT_THROW(newSolution("phase-reaction-spec1.yaml",
                             "nokinetics-reactions"),
                 InputFileError);
}

TEST(KineticsFromYaml, ReactionsFieldWithoutKineticsModel2)
{
    EXPECT_THROW(newSolution("phase-reaction-spec2.yaml",
                             "nokinetics-reactions"),
                 InputFileError);
}

class ReactionToYaml : public testing::Test
{
public:
    void duplicateReaction(size_t i) {
        auto kin = soln->kinetics();
        iOld = i;
        // Exclude the raw input data, to make sure this test actually relies on
        // the fields being populated by the Reaction types.
        AnyMap rdata1 = kin->reaction(iOld)->parameters(false);
        AnyMap rdata2 = AnyMap::fromYamlString(rdata1.toYamlString());
        duplicate = newReaction(rdata2, *kin);
        kin->addReaction(duplicate);
        iNew = kin->nReactions() - 1;
    }

    void compareReactions() {
        auto kin = soln->kinetics();
        EXPECT_EQ(kin->reactionString(iOld), kin->reactionString(iNew));
        EXPECT_EQ(kin->isReversible(iOld), kin->isReversible(iNew));

        vector_fp kf(kin->nReactions()), kr(kin->nReactions());
        vector_fp ropf(kin->nReactions()), ropr(kin->nReactions());
        kin->getFwdRateConstants(kf.data());
        kin->getRevRateConstants(kr.data());
        kin->getFwdRatesOfProgress(ropf.data());
        kin->getRevRatesOfProgress(ropr.data());
        EXPECT_DOUBLE_EQ(kf[iOld], kf[iNew]);
        EXPECT_DOUBLE_EQ(kr[iOld], kr[iNew]);
        EXPECT_DOUBLE_EQ(ropf[iOld], ropf[iNew]);
        EXPECT_DOUBLE_EQ(ropr[iOld], ropr[iNew]);
    }

    shared_ptr<Solution> soln;
    shared_ptr<Reaction> duplicate;

    size_t iOld;
    size_t iNew;
};

TEST_F(ReactionToYaml, elementary)
{
    soln = newSolution("h2o2.yaml", "", "None");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, O2:0.5, O:1e-8, OH:3e-8");
    duplicateReaction(2);
    EXPECT_TRUE(std::dynamic_pointer_cast<ElementaryReaction3>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, threeBody)
{
    soln = newSolution("h2o2.yaml", "", "None");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, O2:0.5, O:1e-8, OH:3e-8, H:2e-7");
    duplicateReaction(1);
    EXPECT_TRUE(std::dynamic_pointer_cast<ThreeBodyReaction3>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, TroeFalloff)
{
    soln = newSolution("h2o2.yaml", "", "None");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, O2:0.5, H2O2:1e-8, OH:3e-8");
    duplicateReaction(21);
    EXPECT_TRUE(std::dynamic_pointer_cast<FalloffReaction>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, SriFalloff)
{
    soln = newSolution("sri-falloff.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "R1A: 0.1, R1B:0.2, H: 0.2, R2:0.5");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<FalloffReaction>(duplicate));
    compareReactions();
    duplicateReaction(1);
    compareReactions();
}

TEST_F(ReactionToYaml, TsangFalloff)
{
    soln = newSolution("tsang-falloff.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "NO:1.0, OH:1.0, H:1.0, CN:1.0");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<FalloffReaction>(duplicate));
    compareReactions();
    duplicateReaction(1);
    compareReactions();
}

TEST_F(ReactionToYaml, chemicallyActivated)
{
    soln = newSolution("chemically-activated-reaction.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "H2:1.0, ch2o:0.1, ch3:1e-8, oh:3e-6");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<ChemicallyActivatedReaction>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, pdepArrhenius)
{
    soln = newSolution("pdep-test.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "R2:1, H:0.1, P2A:2, P2B:0.3");
    duplicateReaction(1);
    EXPECT_TRUE(std::dynamic_pointer_cast<PlogReaction3>(duplicate));
    compareReactions();
    soln->thermo()->setState_TPY(1100, 1e3, "R2:1, H:0.2, P2A:2, P2B:0.3");
    compareReactions();
}

TEST_F(ReactionToYaml, Chebyshev)
{
    soln = newSolution("pdep-test.yaml");
    soln->thermo()->setState_TPY(1000, 2e5, "R6:1, P6A:2, P6B:0.3");
    duplicateReaction(5);
    EXPECT_TRUE(std::dynamic_pointer_cast<ChebyshevReaction3>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, surface)
{
    auto gas = newSolution("diamond.yaml", "gas");
    auto solid = newSolution("diamond.yaml", "diamond");
    soln = newSolution("diamond.yaml", "diamond_100", "None", {gas, solid});
    auto surf = std::dynamic_pointer_cast<SurfPhase>(soln->thermo());
    surf->setCoveragesByName("c6HH:0.1, c6H*:0.6, c6**:0.1");
    gas->thermo()->setMassFractionsByName("H2:0.7, CH4:0.3");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<InterfaceReaction>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, electrochemical)
{
    auto gas = newSolution("sofc.yaml", "gas");
    auto metal = newSolution("sofc.yaml", "metal");
    auto oxide_bulk = newSolution("sofc.yaml", "oxide_bulk");
    auto metal_surf = newSolution("sofc.yaml", "metal_surface", "None", {gas});
    auto oxide_surf = newSolution("sofc.yaml", "oxide_surface", "None",
                                  {gas, oxide_bulk});
    soln = newSolution("sofc.yaml", "tpb", "None",
                       {metal, metal_surf, oxide_surf});
    auto ox_surf = std::dynamic_pointer_cast<SurfPhase>(oxide_surf->thermo());
    oxide_bulk->thermo()->setElectricPotential(-3.4);
    oxide_surf->thermo()->setElectricPotential(-3.4);
    ox_surf->setCoveragesByName("O''(ox):0.2, OH'(ox):0.3, H2O(ox):0.5");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<ElectrochemicalReaction>(duplicate));
    compareReactions();
    compareReactions();
}

TEST_F(ReactionToYaml, unconvertible1)
{
    ElementaryReaction2 R({{"H2", 1}, {"OH", 1}},
                          {{"H2O", 1}, {"H", 1}},
                          Arrhenius(1e5, -1.0, 12.5));
    AnyMap params = R.parameters();
    UnitSystem U{"g", "cm", "mol"};
    params.setUnits(U);
    EXPECT_THROW(params.applyUnits(), CanteraError);
}

TEST_F(ReactionToYaml, unconvertible2)
{
    Array2D coeffs(2, 2, 1.0);
    ChebyshevReaction2 R({{"H2", 1}, {"OH", 1}},
                         {{"H2O", 1}, {"H", 1}},
                         Chebyshev(273., 3000., 1.e2, 1.e7, coeffs));
    UnitSystem U{"g", "cm", "mol"};
    AnyMap params = R.parameters();
    params.setUnits(U);
    EXPECT_THROW(params.applyUnits(), CanteraError);
}

TEST_F(ReactionToYaml, BlowersMasel)
{
    soln = newSolution("BM_test.yaml", "gas");
    soln->thermo()->setState_TPY(1100, 0.1 * OneAtm, "O:0.01, H2:0.8, O2:0.19");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<BlowersMaselReaction>(duplicate));
    compareReactions();
}

TEST_F(ReactionToYaml, BlowersMaselInterface)
{
    auto gas = newSolution("BM_test.yaml", "gas");
    soln = newSolution("BM_test.yaml", "Pt_surf", "None", {gas});
    gas->thermo()->setState_TPY(1100, 0.1 * OneAtm, "O:0.01, H2:0.8, O2:0.19");
    soln->thermo()->setState_TP(1100, 0.1 * OneAtm);
    auto surf = std::dynamic_pointer_cast<SurfPhase>(soln->thermo());
    surf->setCoveragesByName("H(S):0.1, PT(S):0.8, H2O(S):0.1");
    duplicateReaction(0);
    EXPECT_TRUE(std::dynamic_pointer_cast<BlowersMaselInterfaceReaction>(duplicate));
    compareReactions();
}
