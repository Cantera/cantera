#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace Cantera;
typedef std::vector<std::string> strvec;

class ThermoToYaml : public testing::Test
{
public:
    void setup(const std::string& fileName, const std::string& phaseName="") {
        thermo.reset(newPhase(fileName, phaseName));
        // Because ThermoPhase::input may already contain the data we are trying
        // to check for here, clear it so that the only parameters are those
        // added by the overrides of getParameters.
        thermo->input().clear();
        thermo->getParameters(data);

        speciesData.resize(thermo->nSpecies());
        eosData.resize(thermo->nSpecies());
        for (size_t k = 0; k < thermo->nSpecies(); k++) {
            thermo->getSpeciesParameters(thermo->speciesName(k), speciesData[k]);
            if (speciesData[k].hasKey("equation-of-state")) {
                // Get the first EOS node, for convenience
                eosData[k] = speciesData[k]["equation-of-state"].asVector<AnyMap>()[0];
            }
        }
    }

    shared_ptr<ThermoPhase> thermo;
    AnyMap data;
    std::vector<AnyMap> speciesData;
    std::vector<AnyMap> eosData;
};

TEST_F(ThermoToYaml, simpleIdealGas)
{
    setup("ideal-gas.yaml", "simple");
    thermo->setState_TP(1010, 2e5);
    double rho = thermo->density();
    thermo->getParameters(data);

    ASSERT_EQ(data["thermo"], "ideal-gas");
    ASSERT_EQ(data["state"]["T"], 1010);
    ASSERT_EQ(data["state"]["density"], rho);
}

TEST_F(ThermoToYaml, IdealSolidSoln)
{
    setup("thermo-models.yaml", "IdealSolidSolnPhase2");
    EXPECT_EQ(data["name"], "IdealSolidSolnPhase2");
    EXPECT_EQ(data["species"].asVector<std::string>().size(), thermo->nSpecies());
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

TEST_F(ThermoToYaml, Lattice)
{
    setup("thermo-models.yaml", "Li7Si3-interstitial");
    EXPECT_DOUBLE_EQ(data["site-density"].asDouble(), 1.046344e+01);
    EXPECT_DOUBLE_EQ(eosData[0]["molar-volume"].asDouble(), 0.2);
    EXPECT_EQ(eosData[1].size(), (size_t) 0);
}

TEST_F(ThermoToYaml, LatticeSolid)
{
    setup("thermo-models.yaml", "Li7Si3_and_interstitials");
    EXPECT_DOUBLE_EQ(data["composition"]["Li7Si3(s)"].asDouble(), 1.0);
    EXPECT_DOUBLE_EQ(data["composition"]["Li7Si3-interstitial"].asDouble(), 1.0);
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

TEST_F(ThermoToYaml, IonsFromNeutral)
{
    setup("thermo-models.yaml", "ions-from-neutral-molecule");
    EXPECT_EQ(data["neutral-phase"], "KCl-neutral");
}

TEST_F(ThermoToYaml, Margules)
{
    setup("thermo-models.yaml", "molten-salt-Margules");
    auto& interactions = data["interactions"].asVector<AnyMap>();
    EXPECT_EQ(interactions.size(), (size_t) 1);
    EXPECT_EQ(interactions[0]["species"].asVector<std::string>()[0], "KCl(l)");
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

TEST_F(ThermoToYaml, MaskellSolidSolution)
{
    setup("thermo-models.yaml", "MaskellSolidSoln");
    EXPECT_EQ(data["product-species"], "H(s)");
    EXPECT_DOUBLE_EQ(data["excess-enthalpy"].asDouble(), 5e3);
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
}

TEST_F(ThermoToYaml, DebyeHuckel_beta_ij)
{
    setup("thermo-models.yaml", "debye-huckel-beta_ij");
    EXPECT_EQ(data["activity-data"]["model"], "beta_ij");
    EXPECT_TRUE(data["activity-data"]["use-Helgeson-fixed-form"].asBool());
    auto& beta = data["activity-data"]["beta"].asVector<AnyMap>();
    ASSERT_EQ(beta.size(), (size_t) 3);
    for (size_t i = 0; i < 3; i++) {
        auto species = beta[i]["species"].asVector<std::string>();
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
        auto species = item["species"].asVector<std::string>();
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
}

TEST_F(ThermoToYaml, HMWSoln2)
{
    setup("thermo-models.yaml", "HMW-bogus");
    EXPECT_EQ(data["activity-data"]["temperature-model"], "linear");
    auto& interactions = data["activity-data"]["interactions"].asVector<AnyMap>();
    EXPECT_EQ(interactions.size(), (size_t) 4);
    for (auto& item : interactions) {
        auto species = item["species"].asVector<std::string>();
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

TEST_F(ThermoToYaml, IdealMolalSolution)
{
    setup("thermo-models.yaml", "ideal-molal-aqueous");
    auto& cutoff = data["cutoff"];
    EXPECT_EQ(cutoff["model"], "polyexp");
    EXPECT_EQ(cutoff.as<AnyMap>().size(), (size_t) 2); // other values are defaults
    EXPECT_DOUBLE_EQ(cutoff["gamma_o"].asDouble(), 0.0001);
}
