#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Elements.h"
#include "cantera/thermo/MolalityVPSSTP.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/SurfPhase.h"
#include <fstream>

using namespace Cantera;

namespace {

shared_ptr<ThermoPhase> newThermo(const std::string& fileName,
                                  const std::string& phaseName)
{
    return shared_ptr<ThermoPhase>(newPhase(fileName, phaseName));
}

} // namespace

TEST(ThermoFromYaml, simpleIdealGas)
{
    IdealGasPhase thermo("ideal-gas.yaml", "simple");
    EXPECT_EQ(thermo.nSpecies(), (size_t) 3);
    EXPECT_DOUBLE_EQ(thermo.density(), 7.0318220966379288);
    EXPECT_DOUBLE_EQ(thermo.cp_mass(), 1037.7546065787594);
}

TEST(ThermoFromYaml, failDuplicateSpecies)
{
    EXPECT_THROW(newThermo("ideal-gas.yaml", "duplicate-species"), CanteraError);
}

TEST(ThermoFromYaml, elementOverride)
{
    auto thermo = newThermo("ideal-gas.yaml", "element-override");
    EXPECT_EQ(thermo->nElements(), (size_t) 3);
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(0), getElementWeight("N"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(1), getElementWeight("O"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(2), 36);
}

TEST(ThermoFromYaml, elementFromDifferentFile)
{
    auto thermo = newThermo("ideal-gas.yaml", "element-remote");
    EXPECT_EQ(thermo->nElements(), (size_t) 3);
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(0), getElementWeight("N"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(1), getElementWeight("O"));
    EXPECT_DOUBLE_EQ(thermo->atomicWeight(2), 38);
}

TEST(ThermoFromYaml, speciesFromDifferentFile)
{
    IdealGasPhase thermo("ideal-gas.yaml", "species-remote");
    EXPECT_EQ(thermo.nElements(), (size_t) 2);
    EXPECT_EQ(thermo.nSpecies(), (size_t) 4);
    EXPECT_EQ(thermo.species(0)->composition["O"], 2);
    EXPECT_EQ(thermo.species(3)->composition["O"], 1);
    EXPECT_EQ(thermo.species(2)->name, "NO2");
    EXPECT_DOUBLE_EQ(thermo.moleFraction(3), 0.3);
}

TEST(ThermoFromYaml, speciesAll)
{
    auto thermo = newThermo("ideal-gas.yaml", "species-all");
    EXPECT_EQ(thermo->nElements(), (size_t) 3);
    EXPECT_EQ(thermo->nSpecies(), (size_t) 6);
    EXPECT_EQ(thermo->species(1)->name, "NO");
    EXPECT_EQ(thermo->species(2)->name, "N2");
}

TEST(ThermoFromYaml, StoichSubstance1)
{
    auto thermo = newThermo("thermo-models.yaml", "NaCl(s)");
    EXPECT_EQ(thermo->type(), "StoichSubstance");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    EXPECT_EQ(thermo->nElements(), (size_t) 2);
    EXPECT_DOUBLE_EQ(thermo->density(), 2165.0);
    EXPECT_DOUBLE_EQ(thermo->cp_mass(), 864.88371960557095); // Regression test based on XML
}

TEST(ThermoFromYaml, StoichSubstance2)
{
    auto thermo = newThermo("thermo-models.yaml", "KCl(s)");
    EXPECT_EQ(thermo->type(), "StoichSubstance");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    EXPECT_EQ(thermo->nElements(), (size_t) 2);
    EXPECT_NEAR(thermo->density(), 1980, 0.1);
}

TEST(ThermoFromYaml, SurfPhase)
{
    auto thermo = newThermo("surface-phases.yaml", "Pt-surf");
    EXPECT_EQ(thermo->type(), "Surf");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 3);
    auto surf = std::dynamic_pointer_cast<SurfPhase>(thermo);
    EXPECT_DOUBLE_EQ(surf->siteDensity(), 2.7063e-8);
    vector_fp cov(surf->nSpecies());
    surf->getCoverages(cov.data());
    EXPECT_DOUBLE_EQ(cov[surf->speciesIndex("Pt(s)")], 0.5);
    EXPECT_DOUBLE_EQ(cov[surf->speciesIndex("H(s)")], 0.4);
}

TEST(ThermoFromYaml, EdgePhase)
{
    auto thermo = newThermo("surface-phases.yaml", "TPB");
    EXPECT_EQ(thermo->type(), "Edge");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    auto edge = std::dynamic_pointer_cast<SurfPhase>(thermo);
    EXPECT_DOUBLE_EQ(edge->siteDensity(), 5e-18);
}

TEST(ThermoFromYaml, WaterSSTP)
{
    auto thermo = newThermo("thermo-models.yaml", "liquid-water");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    thermo->setState_TP(350, 2*OneAtm);
    // Regression tests based on XML
    EXPECT_NEAR(thermo->density(), 973.7736331, 1e-6);
    EXPECT_NEAR(thermo->enthalpy_mass(), -15649685.52296013, 1e-6);
}

//! @todo Remove after Cantera 2.5 - class FixedChemPotSSTP is deprecated
TEST(ThermoFromYaml, FixedChemPot)
{
    suppress_deprecation_warnings();
    auto thermo = newThermo("thermo-models.yaml", "Li-fixed");
    EXPECT_EQ(thermo->nSpecies(), (size_t) 1);
    double mu;
    thermo->getChemPotentials(&mu);
    EXPECT_DOUBLE_EQ(mu, -2.3e7);
    make_deprecation_warnings_fatal();
}

TEST(ThermoFromYaml, Margules)
{
    auto thermo = newThermo("thermo-models.yaml", "molten-salt-Margules");
    EXPECT_EQ(thermo->type(), "Margules");

    // Regression test based on LiKCl_liquid.xml
    EXPECT_NEAR(thermo->density(), 2041.9831422315351, 1e-9);
    EXPECT_NEAR(thermo->gibbs_mass(), -9683614.0890585743, 1e-5);
    EXPECT_NEAR(thermo->cp_mole(), 67478.48085733457, 1e-8);
}

TEST(ThermoFromYaml, IdealMolalSoln)
{
    auto thermo = newThermo("thermo-models.yaml", "ideal-molal-aqueous");
    EXPECT_EQ(thermo->type(), "IdealMolalSoln");

    EXPECT_NEAR(thermo->enthalpy_mole(), 0.013282, 1e-6);
    EXPECT_NEAR(thermo->gibbs_mole(), -3.8986e7, 1e3);
    EXPECT_NEAR(thermo->density(), 12.058, 1e-3);
}

TEST(ThermoFromYaml, DebyeHuckel_bdot_ak)
{
    auto thermo = newThermo("thermo-models.yaml", "debye-huckel-B-dot-ak");

    // Regression test based on XML input file
    EXPECT_EQ(thermo->type(), "DebyeHuckel");
    EXPECT_NEAR(thermo->density(), 60.296, 1e-2);
    EXPECT_NEAR(thermo->cp_mass(), 1.58216e5, 1e0);
    EXPECT_NEAR(thermo->entropy_mass(), 4.04233e3, 1e-2);

    vector_fp actcoeff(thermo->nSpecies());
    vector_fp mu_ss(thermo->nSpecies());
    auto& molphase = dynamic_cast<MolalityVPSSTP&>(*thermo);
    molphase.getMolalityActivityCoefficients(actcoeff.data());
    thermo->getStandardChemPotentials(mu_ss.data());
    double act_ref[] = {0.849231, 1.18392, 0.990068, 1.69245, 1.09349, 1.0};
    double mu_ss_ref[] = {-3.06816e+08, -2.57956e+08, -1.84117e+08, 0.0,
        -2.26855e+08, -4.3292e+08};
    for (size_t k = 0; k < thermo->nSpecies(); k++) {
        EXPECT_NEAR(actcoeff[k], act_ref[k], 1e-5);
        EXPECT_NEAR(mu_ss[k], mu_ss_ref[k], 1e3);
    }
}

TEST(ThermoFromYaml, DebyeHuckel_beta_ij)
{
    auto thermo = newThermo("thermo-models.yaml", "debye-huckel-beta_ij");

    // Regression test based on XML input file
    EXPECT_EQ(thermo->type(), "DebyeHuckel");
    EXPECT_NEAR(thermo->density(), 122.262, 1e-3);
    EXPECT_NEAR(thermo->cp_mass(), 81263.5, 1e-1);
    EXPECT_NEAR(thermo->entropy_mass(), 4022.35, 1e-2);

    vector_fp actcoeff(thermo->nSpecies());
    vector_fp mu_ss(thermo->nSpecies());
    auto& molphase = dynamic_cast<MolalityVPSSTP&>(*thermo);
    molphase.getMolalityActivityCoefficients(actcoeff.data());
    thermo->getStandardChemPotentials(mu_ss.data());
    double act_ref[] = {0.959912, 1.16955, 1.16955, 2.40275, 0.681552, 1.0};
    double mu_ss_ref[] = {-3.06816e+08, -2.57956e+08, -1.84117e+08, 0,
        -2.26855e+08, -4.3292e+08};
    for (size_t k = 0; k < thermo->nSpecies(); k++) {
        EXPECT_NEAR(actcoeff[k], act_ref[k], 1e-5);
        EXPECT_NEAR(mu_ss[k], mu_ss_ref[k], 1e3);
    }
}

TEST(ThermoFromYaml, IonsFromNeutral)
{
    auto thermo = newThermo("thermo-models.yaml", "ions-from-neutral-molecule");
    ASSERT_EQ((int) thermo->nSpecies(), 2);
    vector_fp mu(thermo->nSpecies());
    thermo->getChemPotentials(mu.data());

    // Values for regression testing only -- same as "fromScratch" test
    EXPECT_NEAR(thermo->density(), 1984.2507319669949, 1e-6);
    EXPECT_NEAR(thermo->enthalpy_mass(), -14738312.44316336, 1e-6);
    EXPECT_NEAR(mu[0], -4.66404010e+08, 1e1);
    EXPECT_NEAR(mu[1], -2.88157316e+06, 1e-1);
}

TEST(ThermoFromYaml, IonsFromNeutral_fromString)
{
    // A little different because we can't re-read the input file to get the
    // phase definition for the neutral phase
    std::ifstream infile("../data/thermo-models.yaml");
    std::stringstream buffer;
    buffer << infile.rdbuf();
    AnyMap input = AnyMap::fromYamlString(buffer.str());
    auto thermo = newPhase(
        input["phases"].getMapWhere("name", "ions-from-neutral-molecule"),
        input);
    ASSERT_EQ((int) thermo->nSpecies(), 2);
    vector_fp mu(thermo->nSpecies());
    thermo->getChemPotentials(mu.data());

    // Values for regression testing only -- same as "fromScratch" test
    EXPECT_NEAR(thermo->density(), 1984.2507319669949, 1e-6);
    EXPECT_NEAR(thermo->enthalpy_mass(), -14738312.44316336, 1e-6);
    EXPECT_NEAR(mu[0], -4.66404010e+08, 1e1);
    EXPECT_NEAR(mu[1], -2.88157316e+06, 1e-1);
}

//! @todo Remove after Cantera 2.5 - "gas" mode of IdealSolnGasVPSS is
//!     deprecated
TEST(ThermoFromYaml, IdealSolnGas_gas)
{
    suppress_deprecation_warnings();
    auto thermo = newThermo("thermo-models.yaml", "IdealSolnGas-gas");
    make_deprecation_warnings_fatal();
    thermo->equilibrate("HP");
    EXPECT_NEAR(thermo->temperature(), 479.929, 1e-3); // based on h2o2.cti
    EXPECT_NEAR(thermo->moleFraction("H2O"), 0.01, 1e-4);
    EXPECT_NEAR(thermo->moleFraction("H2"), 0.0, 1e-4);
}

TEST(ThermoFromYaml, IdealSolnGas_liquid)
{
    auto thermo = newThermo("thermo-models.yaml", "IdealSolnGas-liquid");
    thermo->setState_TP(300, OneAtm);
    EXPECT_NEAR(thermo->density(), 505.42393940, 2e-8);
    EXPECT_NEAR(thermo->gibbs_mole(), -7801634.1184443515, 2e-8);
    thermo->setState_TP(400, 2*OneAtm);
    EXPECT_NEAR(thermo->density(), 495.06986080, 2e-8);
    EXPECT_NEAR(thermo->molarVolume(), 0.014018223587243668, 2e-12);
    thermo->setState_TP(500, 2*OneAtm);
    EXPECT_NEAR(thermo->density(), 484.66590, 2e-8);
    EXPECT_NEAR(thermo->enthalpy_mass(), 1236701.0904197122, 2e-8);
    EXPECT_NEAR(thermo->entropy_mole(), 49848.488477407751, 2e-8);
}

TEST(ThermoFromYaml, RedlichKister)
{
    auto thermo = newThermo("thermo-models.yaml", "Redlich-Kister-LiC6");
    vector_fp chemPotentials(2);
    vector_fp dlnActCoeffdx(2);
    thermo->setState_TP(298.15, OneAtm);
    thermo->setMoleFractionsByName("Li(C6): 0.6375, V(C6): 0.3625");
    thermo->getChemPotentials(chemPotentials.data());
    thermo->getdlnActCoeffdlnX_diag(dlnActCoeffdx.data());
    EXPECT_NEAR(chemPotentials[0], -1.2618554573674981e+007, 1e-6);
    EXPECT_NEAR(dlnActCoeffdx[0], 0.200612, 1e-6);

    thermo->setMoleFractionsByName("Li(C6): 0.8625, V(C6): 0.1375");
    thermo->getChemPotentials(chemPotentials.data());
    thermo->getdlnActCoeffdlnX_diag(dlnActCoeffdx.data());
    EXPECT_NEAR(chemPotentials[0], -1.179299486233677e+07, 1e-6);
    EXPECT_NEAR(dlnActCoeffdx[0], -0.309379, 1e-6);
}

TEST(ThermoFromYaml, MaskellSolidSoln)
{
    auto thermo = newThermo("thermo-models.yaml", "MaskellSolidSoln");
    vector_fp chemPotentials(2);
    thermo->getChemPotentials(chemPotentials.data());
    EXPECT_NEAR(chemPotentials[0], -4.989677789060059e6, 1e-6);
    EXPECT_NEAR(chemPotentials[1], 4.989677789060059e6 + 1000, 1e-6);
}

TEST(ThermoFromYaml, HMWSoln)
{
    auto thermo = newThermo("thermo-models.yaml", "HMW-NaCl-electrolyte");
    size_t N = thermo->nSpecies();
    auto HMW = dynamic_cast<MolalityVPSSTP*>(thermo.get());
    vector_fp acMol(N), mf(N), activities(N), moll(N), mu0(N);
    thermo->getMoleFractions(mf.data());
    HMW->getMolalities(moll.data());
    HMW->getMolalityActivityCoefficients(acMol.data());
    thermo->getActivities(activities.data());
    thermo->getStandardChemPotentials(mu0.data());

    double acMolRef[] = {0.9341, 1.0191, 3.9637, 1.0191, 0.4660};
    double mfRef[] = {0.8198, 0.0901, 0.0000, 0.0901, 0.0000};
    double activitiesRef[] = {0.7658, 6.2164, 0.0000, 6.2164, 0.0000};
    double mollRef[] = {55.5093, 6.0997, 0.0000, 6.0997, 0.0000};
    double mu0Ref[] = {-317.175792, -186.014569, 0.0017225, -441.615456, -322.000432}; // kJ/gmol

    for (size_t k = 0 ; k < N; k++) {
        EXPECT_NEAR(acMol[k], acMolRef[k], 2e-4);
        EXPECT_NEAR(mf[k], mfRef[k], 2e-4);
        EXPECT_NEAR(activities[k], activitiesRef[k], 2e-4);
        EXPECT_NEAR(moll[k], mollRef[k], 2e-4);
        EXPECT_NEAR(mu0[k]/1e6, mu0Ref[k], 2e-6);
    }
    EXPECT_EQ("liquid", HMW->phaseOfMatter());
}

TEST(ThermoFromYaml, HMWSoln_HKFT)
{
    auto thermo = newThermo("thermo-models.yaml", "HMW-NaCl-HKFT");
    double mvRef[] = {0.01815224, 0.00157182, 0.01954605, 0.00173137, -0.0020266};
    double hRef[] = {-2.84097587e+08, -2.38159643e+08, -1.68846908e+08,
                     3.59728865e+06, -2.29291570e+08};
    double acoeffRef[] = {0.922403480, 1.21859875, 1.21859855, 5.08171133,
                          0.5983205};

    // Regression test based on HMWSoln.fromScratch_HKFT
    size_t N = thermo->nSpecies();
    vector_fp mv(N), h(N), acoeff(N);
    thermo->getPartialMolarVolumes(mv.data());
    thermo->getPartialMolarEnthalpies(h.data());
    thermo->getActivityCoefficients(acoeff.data());
    for (size_t k = 0; k < N; k++) {
        EXPECT_NEAR(mv[k], mvRef[k], 2e-8);
        EXPECT_NEAR(h[k], hRef[k], 2e0);
        EXPECT_NEAR(acoeff[k], acoeffRef[k], 2e-8);
    }
}

TEST(ThermoFromYaml, RedlichKwong_CO2)
{
    auto thermo = newThermo("thermo-models.yaml", "CO2-RK");
    EXPECT_NEAR(thermo->density(), 892.404657616, 1e-8);
    EXPECT_NEAR(thermo->enthalpy_mass(), -9199911.5290408, 1e-6);
    EXPECT_NEAR(thermo->cp_mass(), 2219.940330064, 1e-8);

    thermo->setState_TPX(350, 180*OneAtm, "CO2:0.6, H2O:0.02, H2:0.38");
    EXPECT_NEAR(thermo->density(), 181.564971902, 1e-8);
    EXPECT_NEAR(thermo->enthalpy_mass(), -8873033.2793978, 1e-6);
    EXPECT_NEAR(thermo->cp_mass(), 3358.492543261, 1e-8);
}

TEST(ThermoFromYaml, PureFluid_nitrogen)
{
    auto thermo = newThermo("thermo-models.yaml", "nitrogen");
    thermo->setState_TP(70, 2*OneAtm);
    EXPECT_NEAR(thermo->density(), 841.0420151, 1e-6);
    EXPECT_NEAR(thermo->gibbs_mole(), -17654454.0912211, 1e-6);
}

TEST(ThermoFromYaml, PureFluid_CO2)
{
    auto thermo = newThermo("thermo-models.yaml", "CO2-purefluid");
    EXPECT_NEAR(thermo->vaporFraction(), 0.1, 1e-6);
    EXPECT_NEAR(thermo->density(), 513.27928388, 1e-6);
}

TEST(ThermoFromYaml, PureFluid_Unknown)
{
    AnyMap root = AnyMap::fromYamlString(
        "phases:\n"
        "- name: unknown-purefluid\n"
        "  species: [N2]\n"
        "  thermo: pure-fluid\n"
        "  pure-fluid-name: unknown-purefluid\n"
    );
    AnyMap& phase = root["phases"].getMapWhere("name", "unknown-purefluid");
    EXPECT_THROW(newPhase(phase, root), CanteraError);
}

TEST(ThermoFromYaml, ConstDensityThermo)
{
    suppress_deprecation_warnings();
    auto thermo = newThermo("thermo-models.yaml", "const-density");
    EXPECT_DOUBLE_EQ(thermo->density(), 700.0);
    make_deprecation_warnings_fatal();
}

TEST(ThermoFromYaml, IdealSolidSolnPhase)
{
    auto thermo = newThermo("thermo-models.yaml", "IdealSolidSolnPhase");

    // Regression test following IdealSolidSolnPhase.fromScratch
    EXPECT_NEAR(thermo->density(), 10.1787080, 1e-6);
    EXPECT_NEAR(thermo->enthalpy_mass(), -15642788.8547624, 1e-4);
    EXPECT_NEAR(thermo->gibbs_mole(), -313642312.7114608, 1e-4);
}

TEST(ThermoFromYaml, Lattice)
{
    auto thermo = newThermo("thermo-models.yaml", "Li7Si3_and_interstitials");

    // Regression test based on modified version of Li7Si3_ls.xml
    EXPECT_NEAR(thermo->enthalpy_mass(), -2077955.0584538165, 1e-6);
    double mu_ref[] = {-4.62717474e+08, -4.64248485e+07, 1.16370186e+05};
    double vol_ref[] = {0.095564748201438871, 0.2, 0.09557086};
    vector_fp mu(thermo->nSpecies());
    vector_fp vol(thermo->nSpecies());
    thermo->getChemPotentials(mu.data());
    thermo->getPartialMolarVolumes(vol.data());

    for (size_t k = 0; k < thermo->nSpecies(); k++) {
        EXPECT_NEAR(mu[k], mu_ref[k], 1e-7*fabs(mu_ref[k]));
        EXPECT_NEAR(vol[k], vol_ref[k], 1e-7);
    }
}

TEST(ThermoFromYaml, Lattice_fromString)
{
    // A little different because we can't re-read the input file to get the
    // phase definition for the neutral phase
    std::ifstream infile("../data/thermo-models.yaml");
    std::stringstream buffer;
    buffer << infile.rdbuf();
    AnyMap input = AnyMap::fromYamlString(buffer.str());
    auto thermo = newPhase(
        input["phases"].getMapWhere("name", "Li7Si3_and_interstitials"),
        input);

    // Regression test based on modified version of Li7Si3_ls.xml
    EXPECT_NEAR(thermo->enthalpy_mass(), -2077955.0584538165, 1e-6);
    double mu_ref[] = {-4.62717474e+08, -4.64248485e+07, 1.16370186e+05};
    double vol_ref[] = {0.095564748201438871, 0.2, 0.09557086};
    vector_fp mu(thermo->nSpecies());
    vector_fp vol(thermo->nSpecies());
    thermo->getChemPotentials(mu.data());
    thermo->getPartialMolarVolumes(vol.data());

    for (size_t k = 0; k < thermo->nSpecies(); k++) {
        EXPECT_NEAR(mu[k], mu_ref[k], 1e-7*fabs(mu_ref[k]));
        EXPECT_NEAR(vol[k], vol_ref[k], 1e-7);
    }
}

TEST(ThermoFromYaml, Metal)
{
    auto thermo = newThermo("thermo-models.yaml", "Metal");
    EXPECT_DOUBLE_EQ(thermo->density(), 9.0);
    EXPECT_DOUBLE_EQ(thermo->gibbs_mass(), 0.0);
}

TEST(ThermoFromYaml, BinarySolutionTabulatedThermo)
{
    auto thermo = newThermo("thermo-models.yaml", "graphite-anode");
    EXPECT_NEAR(thermo->density(), 5031.7, 1e-5);
    EXPECT_NEAR(thermo->enthalpy_mass(), -32501.245047302145, 1e-9);
    EXPECT_NEAR(thermo->entropy_mass(), 90.443481807823474, 1e-12);
    thermo->setMoleFractionsByName("Li[anode]: 0.55, V[anode]: 0.45");
    EXPECT_NEAR(thermo->gibbs_mass(), -87066.246182649265, 1e-9);
}

TEST(ThermoFromYaml, DeprecatedPhase)
{
    // The deprecation warning in this file is turned into an
    // error by make_deprecation_warnings_fatal() called in main()
    // for the test suite.
    EXPECT_THROW(newThermo("gri30.yaml", "gri30_mix"), CanteraError);
}
