#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/YamlWriter.h"
#include "cantera/thermo/Species.h"
#include "cantera/thermo/PlasmaPhase.h"

using namespace Cantera;
typedef vector<string> strvec;

class ThermoToYaml : public testing::Test
{
public:
    void setup(const string& fileName, const string& phaseName="") {
        thermo = newThermo(fileName, phaseName);
        // Because ThermoPhase::input may already contain the data we are trying
        // to check for here, clear it so that the only parameters are those
        // added by the overrides of getParameters.
        thermo->input().clear();
        data = thermo->parameters();
        data.applyUnits();

        speciesData.resize(thermo->nSpecies());
        eosData.resize(thermo->nSpecies());
        for (size_t k = 0; k < thermo->nSpecies(); k++) {
            thermo->getSpeciesParameters(thermo->speciesName(k), speciesData[k]);
            speciesData[k].applyUnits();
            if (speciesData[k].hasKey("equation-of-state")) {
                // Get the first EOS node, for convenience
                eosData[k] = speciesData[k]["equation-of-state"].asVector<AnyMap>()[0];
            }
        }
    }

    shared_ptr<ThermoPhase> thermo;
    AnyMap data;
    vector<AnyMap> speciesData;
    vector<AnyMap> eosData;
};

TEST_F(ThermoToYaml, simpleIdealGas)
{
    setup("ideal-gas.yaml", "simple");
    thermo->setState_TP(1010, 2e5);
    double rho = thermo->density();
    data = thermo->parameters();
    data.applyUnits();

    ASSERT_EQ(data["thermo"], "ideal-gas");
    ASSERT_EQ(data["state"]["T"], 1010);
    ASSERT_EQ(data["state"]["density"], rho);
}

TEST_F(ThermoToYaml, IdealSolidSoln)
{
    setup("thermo-models.yaml", "IdealSolidSolnPhase2");
    EXPECT_EQ(data["name"], "IdealSolidSolnPhase2");
    EXPECT_EQ(data["species"].asVector<string>().size(), thermo->nSpecies());
    EXPECT_EQ(data["standard-concentration-basis"], "solvent-molar-volume");

    EXPECT_DOUBLE_EQ(eosData[0]["molar-volume"].asDouble(), 1.5);
    EXPECT_DOUBLE_EQ(eosData[2]["molar-volume"].asDouble(), 0.1);
}

TEST_F(ThermoToYaml, BinarySolutionTabulated)
{
    setup("thermo-models.yaml", "graphite-anode");
    EXPECT_EQ(data["tabulated-species"], "Li[anode]");
    auto& tabThermo = data["tabulated-thermo"].as<AnyMap>();
    auto& X = tabThermo["mole-fractions"].asVector<double>();
    auto& h = tabThermo["enthalpy"].asVector<double>();
    auto& s = tabThermo["entropy"].asVector<double>();
    EXPECT_DOUBLE_EQ(X[0], 5.75e-3);
    EXPECT_DOUBLE_EQ(h[1], -9.69664e6);
    EXPECT_DOUBLE_EQ(s[2], 1.27000e4);
}

TEST_F(ThermoToYaml, StoichSubstance1)
{
    setup("thermo-models.yaml", "NaCl(s)");
    EXPECT_EQ(eosData[0]["model"], "constant-volume");
    EXPECT_DOUBLE_EQ(eosData[0]["density"].asDouble(), 2165.0);
}

TEST_F(ThermoToYaml, StoichSubstance2)
{
    setup("thermo-models.yaml", "KCl(s)");
    EXPECT_EQ(eosData[0]["model"], "constant-volume");
    EXPECT_DOUBLE_EQ(eosData[0]["molar-volume"].asDouble(), 0.0376521717);
}

TEST_F(ThermoToYaml, Metal)
{
    setup("thermo-models.yaml", "Metal");
    EXPECT_EQ(data["thermo"], "electron-cloud");
    EXPECT_DOUBLE_EQ(data["density"].asDouble(), 9.0);
}

TEST_F(ThermoToYaml, PureFluid)
{
    setup("thermo-models.yaml", "nitrogen");
    EXPECT_EQ(data["thermo"], "pure-fluid");
    EXPECT_EQ(data["pure-fluid-name"], "nitrogen");
}

TEST_F(ThermoToYaml, RedlichKwong)
{
    setup("thermo-models.yaml", "CO2-RK");
    auto a = eosData[0]["a"].asVector<double>();
    EXPECT_DOUBLE_EQ(a[0], 7.54e6);
    EXPECT_DOUBLE_EQ(a[1], -4.13e3);
    EXPECT_DOUBLE_EQ(eosData[0]["b"].asDouble(), 27.80e-3);
}

TEST_F(ThermoToYaml, Surface)
{
    setup("surface-phases.yaml", "Pt-surf");
    EXPECT_EQ(data["thermo"], "ideal-surface");
    EXPECT_DOUBLE_EQ(data["site-density"].asDouble(), 2.7063e-8);
}

TEST_F(ThermoToYaml, Edge)
{
    setup("surface-phases.yaml", "TPB");
    EXPECT_EQ(data["thermo"], "edge");
    EXPECT_DOUBLE_EQ(data["site-density"].asDouble(), 5e-18);
}

TEST_F(ThermoToYaml, CoverageDependentSurface)
{
    setup("copt_covdepsurf_example.yaml", "covdep");
    EXPECT_EQ(data["thermo"], "coverage-dependent-surface");
    EXPECT_DOUBLE_EQ(data["site-density"].asDouble(), 2.72e-8);
    EXPECT_DOUBLE_EQ(data["reference-state-coverage"].asDouble(), 0.22);
    UnitSystem us;
    EXPECT_EQ(speciesData[1]["coverage-dependencies"]["OC_Pt"]["model"], "linear");
    EXPECT_DOUBLE_EQ(
        speciesData[1]["coverage-dependencies"]["OC_Pt"]["enthalpy"]
        .asDouble(), us.convertFrom(0.48, "eV/molec"));
    EXPECT_EQ(
        speciesData[2]["coverage-dependencies"]["OC_Pt"]["model"], "piecewise-linear");
    EXPECT_DOUBLE_EQ(
        speciesData[2]["coverage-dependencies"]["OC_Pt"]["entropy-high"]
        .asDouble(), us.convertFrom(-0.1, "eV/molec"));
    EXPECT_EQ(
        speciesData[1]["coverage-dependencies"]["O_Pt"]["model"], "interpolative");
    EXPECT_DOUBLE_EQ(
        speciesData[1]["coverage-dependencies"]["O_Pt"]["enthalpies"]
        .asVector<double>()[3], us.convertFrom(2.7, "kcal/mol"));
    EXPECT_EQ(
        speciesData[4]["coverage-dependencies"]["O_Pt"]["model"], "polynomial");
    EXPECT_DOUBLE_EQ(
        speciesData[4]["coverage-dependencies"]["O_Pt"]["enthalpy-coefficients"]
        .asVector<double>()[3], us.convertFrom(2.11, "eV/molec"));
    EXPECT_DOUBLE_EQ(
        speciesData[3]["coverage-dependencies"]["C_Pt"]["heat-capacity-a"]
        .asDouble(), us.convertFrom(0.07e-3, "eV/molec"));
}

TEST_F(ThermoToYaml, Margules)
{
    setup("thermo-models.yaml", "molten-salt-Margules");
    auto& interactions = data["interactions"].asVector<AnyMap>();
    EXPECT_EQ(interactions.size(), (size_t) 1);
    EXPECT_EQ(interactions[0]["species"].asVector<string>()[0], "KCl(l)");
    EXPECT_EQ(interactions[0]["excess-enthalpy"].asVector<double>()[1], -377e3);
}

TEST_F(ThermoToYaml, RedlichKister)
{
    setup("thermo-models.yaml", "Redlich-Kister-LiC6");
    auto& interactions = data["interactions"].asVector<AnyMap>();
    EXPECT_EQ(interactions.size(), (size_t) 1);
    auto& I = interactions[0];
    EXPECT_EQ(I["excess-enthalpy"].asVector<double>().size(), (size_t) 15);
    EXPECT_EQ(I["excess-entropy"].asVector<double>().size(), (size_t) 1);
}

TEST_F(ThermoToYaml, DebyeHuckel_B_dot_ak)
{
    setup("thermo-models.yaml", "debye-huckel-B-dot-ak");
    auto& ac = data["activity-data"];
    EXPECT_EQ(ac["model"], "B-dot-with-variable-a");
    EXPECT_DOUBLE_EQ(ac["B-dot"].asDouble(), 0.0410);
    EXPECT_DOUBLE_EQ(ac["max-ionic-strength"].asDouble(), 50.0);
    EXPECT_DOUBLE_EQ(ac["default-ionic-radius"].asDouble(), 4e-10);
    EXPECT_FALSE(ac.as<AnyMap>().hasKey("A_Debye"));
    EXPECT_FALSE(ac.as<AnyMap>().hasKey("B_Debye"));

    EXPECT_EQ(eosData[0]["model"], "liquid-water-IAPWS95");
    EXPECT_EQ(eosData[1]["model"], "constant-volume");
    EXPECT_DOUBLE_EQ(eosData[1]["molar-volume"].asDouble(), 1.3);

    EXPECT_FALSE(speciesData[0].hasKey("Debye-Huckel"));
    EXPECT_FALSE(speciesData[1].hasKey("Debye-Huckel")); // defaults are ok
    EXPECT_DOUBLE_EQ(speciesData[2]["Debye-Huckel"]["ionic-radius"].asDouble(), 3e-10);
    EXPECT_DOUBLE_EQ(speciesData[5]["Debye-Huckel"]["weak-acid-charge"].asDouble(), -1);
}

TEST_F(ThermoToYaml, DebyeHuckel_beta_ij)
{
    setup("thermo-models.yaml", "debye-huckel-beta_ij");
    EXPECT_EQ(data["activity-data"]["model"], "beta_ij");
    EXPECT_TRUE(data["activity-data"]["use-Helgeson-fixed-form"].asBool());
    auto& beta = data["activity-data"]["beta"].asVector<AnyMap>();
    ASSERT_EQ(beta.size(), (size_t) 3);
    for (size_t i = 0; i < 3; i++) {
        auto species = beta[i]["species"].asVector<string>();
        std::sort(species.begin(), species.end());
        if (species[0] == "Cl-" && species[1] == "H+") {
            EXPECT_DOUBLE_EQ(beta[i]["beta"].asDouble(), 0.27);
        } else if (species[0] == "Cl-" && species[1] == "Na+") {
            EXPECT_DOUBLE_EQ(beta[i]["beta"].asDouble(), 0.15);
        } else {
            EXPECT_EQ(species[0], "Na+");
            EXPECT_EQ(species[1], "OH-");
            EXPECT_DOUBLE_EQ(beta[i]["beta"].asDouble(), 0.06);
        }
    }
}

TEST_F(ThermoToYaml, HMWSoln1)
{
    setup("thermo-models.yaml", "HMW-NaCl-electrolyte");
    EXPECT_EQ(data["activity-data"]["temperature-model"], "complex");
    auto& interactions = data["activity-data"]["interactions"].asVector<AnyMap>();
    EXPECT_EQ(interactions.size(), (size_t) 7);
    for (auto& item : interactions) {
        auto species = item["species"].asVector<string>();
        std::sort(species.begin(), species.end());
        if (species == strvec{"Cl-", "Na+"}) {
            auto& beta0 = item["beta0"].asVector<double>();
            EXPECT_EQ(beta0.size(), (size_t) 5);
            EXPECT_DOUBLE_EQ(beta0[1], 0.008946);
        } else if (species == strvec{"Cl-", "H+"}) {
            EXPECT_TRUE(item.hasKey("beta2"));
            EXPECT_TRUE(item.hasKey("Cphi"));
        } else if (species == strvec{"Na+", "OH-"}) {
            EXPECT_DOUBLE_EQ(item["beta2"].asDouble(), 0.0);
        } else if (species == strvec{"Cl-", "OH-"}) {
            EXPECT_DOUBLE_EQ(item["theta"].asDouble(), -0.05);
        } else if (species == strvec{"Cl-", "Na+", "OH-"}) {
            EXPECT_DOUBLE_EQ(item["psi"].asDouble(), -0.006);
        } else if (species == strvec{"H+", "Na+"}) {
            EXPECT_DOUBLE_EQ(item["theta"].asDouble(), 0.036);
        } else if (species == strvec{"Cl-", "H+", "Na+"}) {
            EXPECT_DOUBLE_EQ(item["psi"].asDouble(), -0.004);
        } else {
            FAIL(); // unexpected set of species
        }
    }
    EXPECT_EQ(eosData[0]["model"], "liquid-water-IAPWS95");
    EXPECT_EQ(eosData[1]["model"], "constant-volume");
    EXPECT_DOUBLE_EQ(eosData[2]["molar-volume"].asDouble(), 1.3);
}

TEST_F(ThermoToYaml, HMWSoln2)
{
    setup("thermo-models.yaml", "HMW-bogus");
    EXPECT_EQ(data["activity-data"]["temperature-model"], "linear");
    auto& interactions = data["activity-data"]["interactions"].asVector<AnyMap>();
    EXPECT_EQ(interactions.size(), (size_t) 4);
    for (auto& item : interactions) {
        auto species = item["species"].asVector<string>();
        std::sort(species.begin(), species.end());
        if (species == strvec{"Cl-", "NaCl(aq)"}) {
            EXPECT_DOUBLE_EQ(item["lambda"].asVector<double>()[0], 0.3);
        } else if (species == strvec{"Na+", "NaCl(aq)"}) {
            EXPECT_DOUBLE_EQ(item["lambda"].asVector<double>()[1], 0.02);
        } else if (species == strvec{"Na+", "NaCl(aq)", "OH-"}) {
            EXPECT_DOUBLE_EQ(item["zeta"].asVector<double>()[0], 0.5);
        } else if (species == strvec{"NaCl(aq)"}) {
            EXPECT_DOUBLE_EQ(item["mu"].asVector<double>()[1], 0.3);
        } else {
            FAIL(); // unexpected set of species
        }
    }
    auto& crop = data["activity-data"]["cropping-coefficients"];
    EXPECT_DOUBLE_EQ(crop["ln_gamma_k_min"].asDouble(), -8.0);
    EXPECT_DOUBLE_EQ(crop["ln_gamma_k_max"].asDouble(), 20);
}

TEST_F(ThermoToYaml, HMWSoln_HKFT)
{
    setup("thermo-models.yaml", "HMW-NaCl-HKFT");
    EXPECT_DOUBLE_EQ(eosData[1]["h0"].asDouble(), -57433 * 4184);
    EXPECT_DOUBLE_EQ(eosData[1]["s0"].asDouble(), 13.96 * 4184);
    EXPECT_DOUBLE_EQ(eosData[2]["a"].asVector<double>()[2], 5.563 * 4184 / 1e5);
    EXPECT_DOUBLE_EQ(eosData[4]["c"].asVector<double>()[1], -103460 * 4184);
    EXPECT_DOUBLE_EQ(eosData[4]["omega"].asDouble(), 172460 * 4184);
}

TEST_F(ThermoToYaml, IdealMolalSolution)
{
    setup("thermo-models.yaml", "ideal-molal-aqueous");
    auto& cutoff = data["cutoff"];
    EXPECT_EQ(cutoff["model"], "polyexp");
    EXPECT_EQ(cutoff.as<AnyMap>().size(), (size_t) 2); // other values are defaults
    EXPECT_DOUBLE_EQ(cutoff["gamma_o"].asDouble(), 0.0001);

    EXPECT_EQ(eosData[2]["model"], "constant-volume");
    EXPECT_DOUBLE_EQ(eosData[2]["molar-density"].asDouble(), 10);

    EXPECT_EQ(eosData[3]["model"], "constant-volume");
    EXPECT_DOUBLE_EQ(eosData[3]["density"].asDouble(), 160.43);
}

TEST_F(ThermoToYaml, IsotropicElectronEnergyPlasma)
{
    setup("oxygen-plasma.yaml", "isotropic-electron-energy-plasma");
    auto& electronEnergyDist = data["electron-energy-distribution"].as<AnyMap>();
    EXPECT_EQ(electronEnergyDist["type"], "isotropic");
    EXPECT_DOUBLE_EQ(electronEnergyDist["shape-factor"].asDouble(), 2.0);
}

TEST_F(ThermoToYaml, DiscretizedElectronEnergyPlasma)
{
    setup("oxygen-plasma.yaml", "discretized-electron-energy-plasma");
    auto& electronEnergyDist = data["electron-energy-distribution"].as<AnyMap>();
    vector<double> levels = electronEnergyDist["energy-levels"].asVector<double>();
    vector<double> dist = electronEnergyDist["distribution"].asVector<double>();
    EXPECT_EQ(electronEnergyDist["type"], "discretized");
    EXPECT_DOUBLE_EQ(levels[3], 10.0);
    EXPECT_DOUBLE_EQ(dist[3], 0.01);
}


class ThermoYamlRoundTrip : public testing::Test
{
public:
    void roundtrip(const string& fileName, const string& phaseName="",
        const vector<string> extraPhases={}) {
        original = newThermo(fileName, phaseName);
        YamlWriter writer;
        writer.addPhase(original);
        for (const auto& name : extraPhases) {
            shared_ptr<ThermoPhase> p(newThermo(fileName, name));
            writer.addPhase(p);
        }
        writer.skipUserDefined();
        AnyMap input1 = AnyMap::fromYamlString(writer.toYamlString());
        duplicate = newThermo(input1["phases"].getMapWhere("name", phaseName), input1);
        skip_cp = false;
        skip_activities = false;
        skip_entropy = false;
        rtol = 1e-14;
    }

    void compareThermo(double T, double P, const string& X="") {
        size_t kk = original->nSpecies();
        ASSERT_EQ(original->nSpecies(), duplicate->nSpecies());

        if (X.empty()) {
            original->setState_TP(T, P);
            duplicate->setState_TP(T, P);
        } else {
            original->setState_TPX(T, P, X);
            duplicate->setState_TPX(T, P, X);
        }

        double rhoOrig = original->density();
        double rhoDup = duplicate->density();
        if (rhoOrig != rhoDup) {
            EXPECT_NEAR(original->density(), duplicate->density(),
                        rtol * original->density());
        }
        if (!skip_cp) {
            EXPECT_NEAR(original->cp_mass(), duplicate->cp_mass(),
                        rtol * original->cp_mass());
        }
        if (!skip_entropy) {
            EXPECT_NEAR(original->entropy_mass(), duplicate->entropy_mass(),
                        rtol * fabs(original->entropy_mass()));
        }
        EXPECT_NEAR(original->enthalpy_mole(), duplicate->enthalpy_mole(),
                    rtol * fabs(original->enthalpy_mole()));

        vector<double> Y1(kk), Y2(kk), h1(kk), h2(kk), s1(kk), s2(kk);
        vector<double> mu1(kk), mu2(kk), v1(kk), v2(kk), a1(kk), a2(kk);
        original->getMassFractions(Y1);
        duplicate->getMassFractions(Y2);
        original->getPartialMolarEnthalpies(h1);
        duplicate->getPartialMolarEnthalpies(h2);
        original->getPartialMolarEntropies(s1);
        duplicate->getPartialMolarEntropies(s2);
        original->getChemPotentials(mu1);
        duplicate->getChemPotentials(mu2);
        original->getPartialMolarVolumes(v1);
        duplicate->getPartialMolarVolumes(v2);
        if (!skip_activities) {
            original->getActivityCoefficients(a1);
            duplicate->getActivityCoefficients(a2);
        }

        for (size_t k = 0; k < kk; k++) {
            EXPECT_NEAR(Y1[k], Y2[k], 1e-20 + rtol*fabs(Y1[k])) << k;
            EXPECT_NEAR(h1[k], h2[k], 1e-20 + rtol*fabs(h1[k])) << k;
            EXPECT_NEAR(s1[k], s2[k], 1e-20 + rtol*fabs(s1[k])) << k;
            EXPECT_NEAR(mu1[k], mu2[k], 1e-20 + rtol*fabs(mu1[k])) << k;
            EXPECT_NEAR(v1[k], v2[k], 1e-20 + rtol*fabs(v1[k])) << k;
            EXPECT_NEAR(a1[k], a2[k], 1e-20 + rtol*fabs(a1[k])) << k;
        }
    }

    shared_ptr<ThermoPhase> original;
    shared_ptr<ThermoPhase> duplicate;
    bool skip_cp;
    bool skip_activities;
    bool skip_entropy;
    double rtol;
};

TEST_F(ThermoYamlRoundTrip, RedlichKwong)
{
    roundtrip("nDodecane_Reitz.yaml", "nDodecane_RK");
    compareThermo(500, 6e5, "c12h26: 0.2, o2: 0.1, co2: 0.4, c2h2: 0.3");
}

TEST_F(ThermoYamlRoundTrip, RedlichKwong_crit_props)
{
    roundtrip("thermo-models.yaml", "CO2-RK-params");
    compareThermo(400, 1e6, "CO2:0.8, H2O:0.1, H2:0.1");
    auto params = duplicate->species("CO2")->parameters(duplicate.get());
    params.applyUnits();
    double Tc = params["critical-parameters"]["critical-temperature"].asDouble();
    EXPECT_NEAR(Tc, 304.128, 1e-3);
}

TEST_F(ThermoYamlRoundTrip, PengRobinson)
{
    roundtrip("co2_PR_example.yaml");
    compareThermo(400, 20e5, "CO2:0.9, H2O:0.07, CH4:0.03");
}

TEST_F(ThermoYamlRoundTrip, PengRobinson_crit_props)
{
    roundtrip("thermo-models.yaml", "CO2-PR-params");
    // rtol = 1e-13;
    compareThermo(400, 1e6, "CO2:0.8, H2O:0.1, H2:0.1");
    auto params = duplicate->species("CO2")->parameters(duplicate.get());
    params.applyUnits();
    double Tc = params["critical-parameters"]["critical-temperature"].asDouble();
    EXPECT_NEAR(Tc, 304.128, 1e-3);
}

TEST_F(ThermoYamlRoundTrip, BinarySolutionTabulated)
{
    roundtrip("BinarySolutionTabulatedThermo.yaml", "anode");
    compareThermo(310, 2e5, "Li[anode]:0.4, V[anode]:0.6");
}

TEST_F(ThermoYamlRoundTrip, Margules)
{
    roundtrip("LiKCl_liquid.yaml", "MoltenSalt_electrolyte");
    compareThermo(920, 3e5, "KCl(L):0.35, LiCl(L):0.65");
}

TEST_F(ThermoYamlRoundTrip, DebyeHuckel)
{
    roundtrip("thermo-models.yaml", "debye-huckel-B-dot-ak");
    compareThermo(305, 2e5);
}

TEST_F(ThermoYamlRoundTrip, IdealMolalSolution)
{
    roundtrip("thermo-models.yaml", "ideal-molal-aqueous");
    compareThermo(308, 1.1e5, "H2O(l): 0.95, H2S(aq): 0.01, CO2(aq): 0.04");
}

TEST_F(ThermoYamlRoundTrip, IdealSolutionVpss)
{
    roundtrip("thermo-models.yaml", "IdealSolnGas-liquid");
    compareThermo(320, 1.5e5, "Li(l):1.0");
}

TEST_F(ThermoYamlRoundTrip, HMWSoln)
{
    roundtrip("thermo-models.yaml", "HMW-NaCl-electrolyte");
    rtol = 1e-10; // @todo: Determine why more stringent tolerances can't be met
    compareThermo(350.15, 101325,
                  "H2O(L): 0.8198, Na+:0.09, Cl-:0.09, H+:4.4e-6, OH-:4.4e-6");
}

TEST_F(ThermoYamlRoundTrip, PureFluid_Nitrogen)
{
    roundtrip("thermo-models.yaml", "nitrogen");
    compareThermo(90, 19e5);
}

TEST_F(ThermoYamlRoundTrip, RedlichKister)
{
    roundtrip("thermo-models.yaml", "Redlich-Kister-LiC6");
    compareThermo(310, 2e5);
}

TEST_F(ThermoYamlRoundTrip, Surface)
{
    roundtrip("surface-phases.yaml", "Pt-surf");
    skip_activities = true;
    compareThermo(800, 2*OneAtm);
    auto origSurf = std::dynamic_pointer_cast<SurfPhase>(original);
    auto duplSurf = std::dynamic_pointer_cast<SurfPhase>(duplicate);
    EXPECT_DOUBLE_EQ(origSurf->siteDensity(), duplSurf->siteDensity());
}

TEST_F(ThermoYamlRoundTrip, CoverageDependentSurface)
{
    roundtrip("copt_covdepsurf_example.yaml", "covdep");
    skip_activities = true;
    compareThermo(800, 2*OneAtm,
        "Pt: 0.2, OC_Pt: 0.2, CO2_Pt: 0.2, C_Pt: 0.2, O_Pt: 0.2");
}

TEST_F(ThermoYamlRoundTrip, IsotropicElectronEnergyPlasma)
{
    roundtrip("oxygen-plasma.yaml", "isotropic-electron-energy-plasma");
    skip_cp = true; // Not implemented for PlasmaPhase
    skip_entropy = true; // Not implemented for PlasmaPhase
    compareThermo(800, 2*OneAtm);
    auto origPlasma = std::dynamic_pointer_cast<PlasmaPhase>(original);
    auto duplPlasma = std::dynamic_pointer_cast<PlasmaPhase>(duplicate);
    vector<double> origDist(origPlasma->nElectronEnergyLevels());
    vector<double> duplDist(duplPlasma->nElectronEnergyLevels());
    origPlasma->getElectronEnergyLevels(origDist);
    duplPlasma->getElectronEnergyLevels(duplDist);
    EXPECT_DOUBLE_EQ(origDist[2], duplDist[2]);
}

TEST_F(ThermoYamlRoundTrip, DiscretizedElectronEnergyPlasma)
{
    roundtrip("oxygen-plasma.yaml", "discretized-electron-energy-plasma");
    skip_cp = true; // Not implemented for PlasmaPhase
    skip_entropy = true; // Not implemented for PlasmaPhase
    compareThermo(800, 2*OneAtm);
    EXPECT_DOUBLE_EQ(original->electronTemperature(), duplicate->electronTemperature());
}
