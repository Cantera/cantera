#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/kinetics/GasKinetics.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/KineticsFactory.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace Cantera;

TEST(Reaction, ElementaryFromYaml)
{
    auto sol = newSolution("gri30.yaml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: N + NO <=> N2 + O,"
        " rate-constant: [-2.70000E+13 cm^3/mol/s, 0, 355 cal/mol],"
        " negative-A: true}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->reaction_type, ELEMENTARY_RXN);

    auto ER = dynamic_cast<ElementaryReaction&>(*R);
    EXPECT_DOUBLE_EQ(ER.rate.preExponentialFactor(), -2.7e10);
    EXPECT_DOUBLE_EQ(ER.rate.activationEnergy_R(), 355 / GasConst_cal_mol_K);
    EXPECT_TRUE(ER.allow_negative_pre_exponential_factor);
    EXPECT_FALSE(ER.allow_negative_orders);
}

TEST(Reaction, ThreeBodyFromYaml1)
{
    auto sol = newSolution("gri30.yaml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O + M = O2 + M,"
        " type: three-body,"
        " rate-constant: [1.20000E+17 cm^6/mol^2/s, -1, 0],"
        " efficiencies: {AR: 0.83, H2O: 5}}");

    auto R = newReaction(rxn, *(sol->kinetics()));
    EXPECT_EQ(R->reactants.count("M"), (size_t) 0);

    auto TBR = dynamic_cast<ThreeBodyReaction&>(*R);
    EXPECT_DOUBLE_EQ(TBR.rate.preExponentialFactor(), 1.2e11);
    EXPECT_DOUBLE_EQ(TBR.third_body.efficiencies["H2O"], 5.0);
    EXPECT_DOUBLE_EQ(TBR.third_body.default_efficiency, 1.0);
}

TEST(Reaction, ThreeBodyFromYaml2)
{
    auto sol = newSolution("gri30.yaml");
    AnyMap rxn = AnyMap::fromYamlString(
        "{equation: 2 O <=> O2," // Missing "M" on each side of the equation
        " type: three-body,"
        " rate-constant: [1.20000E+17, -1, 0]}");

    EXPECT_THROW(newReaction(rxn, *(sol->kinetics())), CanteraError);
}

TEST(Reaction, FalloffFromYaml1)
{
    auto sol = newSolution("gri30.yaml");
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
    auto sol = newSolution("gri30.yaml");
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

TEST(Reaction, ChemicallyActivatedFromYaml)
{
    auto sol = newSolution("gri30.yaml");
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
    auto sol = newSolution("gri30.yaml");
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
    auto PR = dynamic_cast<PlogReaction&>(*R);
    const auto& rates = PR.rate.rates();
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
    auto sol = newSolution("gri30.yaml");
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
    auto CR = dynamic_cast<ChebyshevReaction&>(*R);
    double logP = std::log10(2e6);
    double T = 1800;
    CR.rate.update_C(&logP);
    EXPECT_EQ(CR.rate.nTemperature(), (size_t) 6);
    EXPECT_EQ(CR.rate.nPressure(), (size_t) 4);
    EXPECT_DOUBLE_EQ(CR.rate.Tmax(), 3000);
    EXPECT_DOUBLE_EQ(CR.rate.Pmin(), 1000);
    EXPECT_NEAR(CR.rate.updateRC(std::log(T), 1.0/T), 130512.2773948636, 1e-9);
}

TEST(Kinetics, GasKineticsFromYaml1)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode = infile["phases"].getMapWhere("name", "simple-kinetics");
    shared_ptr<ThermoPhase> thermo = newPhase(phaseNode, infile);
    std::vector<ThermoPhase*> phases{thermo.get()};
    auto kin = newKinetics(phases, phaseNode, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 2);
    const auto& R = kin->reaction(0);
    EXPECT_EQ(R->reactants.at("NO"), 1);
    EXPECT_EQ(R->products.at("N2"), 1);
    EXPECT_EQ(R->id, "NOx-R1");
    const auto& ER = std::dynamic_pointer_cast<ElementaryReaction>(R);
    EXPECT_DOUBLE_EQ(ER->rate.preExponentialFactor(), 2.7e10);
}

TEST(Kinetics, GasKineticsFromYaml2)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode = infile["phases"].getMapWhere("name", "remote-kinetics");
    shared_ptr<ThermoPhase> thermo = newPhase(phaseNode, infile);
    std::vector<ThermoPhase*> phases{thermo.get()};
    auto kin = newKinetics(phases, phaseNode, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 3);
}

TEST(Kinetics, EfficienciesFromYaml)
{
    AnyMap infile = AnyMap::fromYamlFile("ideal-gas.yaml");
    auto& phaseNode1 = infile["phases"].getMapWhere("name", "efficiency-error");
    shared_ptr<ThermoPhase> thermo1 = newPhase(phaseNode1, infile);
    std::vector<ThermoPhase*> phases1{thermo1.get()};
    // Reaction with efficiency defined for undeclared species "AR"
    EXPECT_THROW(newKinetics(phases1, phaseNode1, infile), CanteraError);

    auto& phaseNode2 = infile["phases"].getMapWhere("name", "efficiency-skip");
    shared_ptr<ThermoPhase> thermo2 = newPhase(phaseNode2, infile);
    std::vector<ThermoPhase*> phases2{thermo2.get()};
    auto kin = newKinetics(phases2, phaseNode2, infile);
    EXPECT_EQ(kin->nReactions(), (size_t) 1);
}

TEST(Kinetics, InterfaceKineticsFromYaml)
{
    shared_ptr<ThermoPhase> gas(newPhase("surface-phases.yaml", "gas"));
    shared_ptr<ThermoPhase> surf_tp(newPhase("surface-phases.yaml", "Pt-surf"));
    shared_ptr<SurfPhase> surf = std::dynamic_pointer_cast<SurfPhase>(surf_tp);
    std::vector<ThermoPhase*> phases{surf_tp.get(), gas.get()};
    auto kin = newKinetics(phases, "surface-phases.yaml", "Pt-surf");
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

TEST(Kinetics, ElectrochemFromYaml)
{
    shared_ptr<ThermoPhase> graphite(newPhase("surface-phases.yaml", "graphite"));
    shared_ptr<ThermoPhase> electrolyte(newPhase("surface-phases.yaml", "electrolyte"));
    shared_ptr<ThermoPhase> anode(newPhase("surface-phases.yaml", "anode-surface"));
    std::vector<ThermoPhase*> phases{anode.get(), graphite.get(), electrolyte.get()};
    auto kin = newKinetics(phases, "surface-phases.yaml", "anode-surface");
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
