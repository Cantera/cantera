#include "gtest/gtest.h"
#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"


namespace Cantera
{

class PengRobinson_Test : public testing::Test
{
public:
    PengRobinson_Test() {
        test_phase.reset(newPhase("../data/thermo-models.yaml", "CO2-PR"));
    }

    //vary the composition of a co2-h2 mixture:
    void set_r(const double r) {
        vector_fp moleFracs(7);
        moleFracs[0] = r;
        moleFracs[2] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(PengRobinson_Test, construct_from_yaml)
{
    PengRobinson* peng_robinson_phase = dynamic_cast<PengRobinson*>(test_phase.get());
    EXPECT_TRUE(peng_robinson_phase != NULL);
}

TEST_F(PengRobinson_Test, chem_potentials)
{
    test_phase->setState_TP(298.15, 101325.);
    /* Chemical potential should increase with increasing co2 mole fraction:
    *      mu = mu_0 + RT ln(gamma_k*X_k).
    *  where gamma_k is the activity coefficient. Run regression test against values
    *  calculated using the model.
    */
    const double expected_result[9] = {
        -457338129.70445037,
        -457327078.87912911,
        -457317214.31077951,
        -457308354.65227401,
        -457300353.74028891,
        -457293092.45485628,
        -457286472.73969948,
        -457280413.14238912,
        -457274845.44186872
    };

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(7);
    for(int i=0; i < numSteps; ++i)
    {
        set_r(xmin + i*dx);
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
    }
}

TEST_F(PengRobinson_Test, chemPotentials_RT)
{
    double T = 410.0;
    test_phase->setState_TP(T, 130 * OneAtm);

    // Test that chemPotentials_RT*RT = chemPotentials
    const double RT = GasConstant * T;
    vector_fp mu(7);
    vector_fp mu_RT(7);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);

    for(int i=0; i < numSteps; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getChemPotentials(&mu[0]);
        test_phase->getChemPotentials_RT(&mu_RT[0]);
        EXPECT_NEAR(mu[0], mu_RT[0]*RT, 1.e-6);
        EXPECT_NEAR(mu[2], mu_RT[2]*RT, 1.e-6);
    }
}

TEST_F(PengRobinson_Test, activityCoeffs)
{
    double T = 330.0;
    test_phase->setState_TP(T, 120 * OneAtm);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * T;
    vector_fp mu0(7);
    vector_fp activityCoeffs(7);
    vector_fp chemPotentials(7);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);

    for(int i=0; i < numSteps; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getChemPotentials(&chemPotentials[0]);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        test_phase->getStandardChemPotentials(&mu0[0]);
        EXPECT_NEAR(chemPotentials[0], mu0[0] + RT*std::log(activityCoeffs[0] * r), 1.e-6);
        EXPECT_NEAR(chemPotentials[2], mu0[2] + RT*std::log(activityCoeffs[2] * (1-r)), 1.e-6);
    }
}

TEST_F(PengRobinson_Test, standardConcentrations)
{
    EXPECT_DOUBLE_EQ(test_phase->pressure()/(test_phase->temperature()*GasConstant),
                     test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(test_phase->pressure()/(test_phase->temperature()*GasConstant),
                     test_phase->standardConcentration(1));
}

TEST_F(PengRobinson_Test, activityConcentrations)
{
    // Check to make sure activityConcentration_i == standardConcentration_i * gamma_i * X_i
    vector_fp standardConcs(7);
    vector_fp activityCoeffs(7);
    vector_fp activityConcentrations(7);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    test_phase->setState_TP(350, 100 * OneAtm);

    for(int i=0; i < numSteps; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        standardConcs[0] = test_phase->standardConcentration(0);
        standardConcs[2] = test_phase->standardConcentration(2);
        test_phase->getActivityConcentrations(&activityConcentrations[0]);

        EXPECT_NEAR(standardConcs[0] * r * activityCoeffs[0], activityConcentrations[0], 1.e-6);
        EXPECT_NEAR(standardConcs[2] * (1-r) * activityCoeffs[2], activityConcentrations[2], 1.e-6);
    }
}

TEST_F(PengRobinson_Test, setTP)
{
    // Check to make sure that the phase diagram is accurately reproduced for a few select isobars

    // All sub-cooled liquid:
    const double rho1[6] = {
        6.6603507723749249e+002,
        1.6824762614489907e+002,
        1.6248354709581241e+002,
        1.5746729362032696e+002,
        1.5302217175386241e+002,
        1.4902908974486667e+002
    };
    // Phase change between temperatures 4 & 5:
    const double rho2[6] = {
        7.5732259810273172e+002,
        7.2766981078381912e+002,
        6.935475475396446e+002,
        6.5227027102964917e+002,
        5.9657442842753153e+002,
        3.9966973266966875e+002
    };
    // Supercritical; no discontinuity in rho values:
    const double rho3[6] = {
        8.0601205067780199e+002,
        7.8427655940884574e+002,
        7.6105347579146576e+002,
        7.3605202492828505e+002,
        7.0887891410210011e+002,
        6.7898591969734434e+002
    };

    for(int i=0; i<6; ++i)
    {
        const double temp = 294 + i*2;
        set_r(0.999);
        test_phase->setState_TP(temp, 5542027.5);
        EXPECT_NEAR(test_phase->density(), rho1[i], 1.e-8);

        test_phase->setState_TP(temp, 7388370.);
        EXPECT_NEAR(test_phase->density(), rho2[i], 1.e-7);

        test_phase->setState_TP(temp, 9236712.5);
        EXPECT_NEAR(test_phase->density(), rho3[i], 1.e-8);
    }
}

TEST_F(PengRobinson_Test, getPressure)
{
    // Check to make sure that the P-R equation is accurately reproduced for a few selected values

    /* This test uses CO2 as the only species.
    *  Values of a_coeff, b_coeff are calculated based on the the critical temperature
    *  and pressure values of CO2 as follows:
    *       a_coeff = 0.457235(RT_crit)^2/p_crit
    *       b_coeff = 0.077796(RT_crit)/p_crit
    *  The temperature dependent parameter in P-R EoS is calculated as
    *       \alpha = [1 + \kappa(1 - sqrt{T/T_crit}]^2
    *  kappa is a function calculated based on the acentric factor.
    */

    double a_coeff = 3.958095109E+5;
    double b_coeff = 26.62616317/1000;
    double acc_factor = 0.228;
    double pres_theoretical, kappa, alpha, mv;
    const double rho = 1.0737;

    //Calculate kappa value
    kappa = 0.37464 + 1.54226*acc_factor - 0.26992*acc_factor*acc_factor;

    for (int i = 0; i < 15; i++)
    {
        const double temp = 296 + i * 50;
        set_r(1.0);
        test_phase->setState_TR(temp, rho);
        const double Tcrit = test_phase->critTemperature();
        mv = 1 / rho * test_phase->meanMolecularWeight();
        //Calculate pressure using Peng-Robinson EoS
        alpha = pow(1 + kappa*(1 - sqrt(temp / Tcrit)), 2.0);
        pres_theoretical = GasConstant*temp / (mv - b_coeff)
                          - a_coeff*alpha / (mv*mv + 2*b_coeff*mv - b_coeff*b_coeff);
        EXPECT_NEAR(test_phase->pressure(), pres_theoretical, 3);
    }
}

TEST_F(PengRobinson_Test, gibbsEnergy)
{
    // Test that g == h - T*s
    const double T = 360.;
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax - xmin) / (numSteps - 1);
    double gibbs_theoretical;

    for (int i = 0; i < numSteps; ++i)
    {
        const double r = xmin + i * dx;
        test_phase->setState_TP(T, 150e5);
        set_r(r);
        gibbs_theoretical = test_phase->enthalpy_mole() - T * (test_phase->entropy_mole());
        EXPECT_NEAR(test_phase->gibbs_mole(), gibbs_theoretical, 1.e-6);
    }
}

TEST_F(PengRobinson_Test, totalEnthalpy)
{
    // Test that hbar = \sum (h_k*x_k)
    double hbar, sum = 0.0;
    vector_fp partialEnthalpies(7);
    vector_fp moleFractions(7);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax - xmin) / (numSteps - 1);

    for (int i = 0; i < numSteps; ++i)
    {
        sum = 0.0;
        const double r = xmin + i * dx;
        test_phase->setState_TP(430., 120e5);
        set_r(r);
        hbar = test_phase->enthalpy_mole();
        test_phase->getMoleFractions(&moleFractions[0]);
        test_phase->getPartialMolarEnthalpies(&partialEnthalpies[0]);
        for (int k = 0; k < 7; k++)
        {
            sum += moleFractions[k] * partialEnthalpies[k];
        }
        EXPECT_NEAR(hbar, sum, 1.e-6);
    }
}

TEST_F(PengRobinson_Test, cpValidate)
{
    // Test that cp = dH/dT at constant pressure using finite difference method

    double p = 200e5;
    double Tmin = 298;
    int numSteps = 20;
    double dT = 1e-4;
    test_phase->setMoleFractionsByName("CO2: 0.7, H2O: 0.1, H2: 0.2");

    for (int i = 0; i < numSteps; ++i) {
        const double T = Tmin + 10 * i;
        test_phase->setState_TP(T - dT, p);
        double h1 = test_phase->enthalpy_mole();  // J/kmol
        test_phase->setState_TP(T, p);
        double cp = test_phase->cp_mole();        // unit is J/kmol/K
        test_phase->setState_TP(T + dT, p);
        double h2 = test_phase->enthalpy_mole();

        double dh_dT = (h2 - h1) / (2 * dT);
        EXPECT_NEAR(cp, dh_dT, 1e-6 * cp);
    }
}

TEST_F(PengRobinson_Test, CoolPropValidate)
{
    // Validate the P-R EoS in Cantera with P-R EoS from CoolProp

    const double rhoCoolProp[10] = {
        9.067928191884574,
        8.318322900591179,
        7.6883521740498155,
        7.150504298001246,
        6.685330199667018,
        6.278630757480957,
        5.919763108091383,
        5.600572727499541,
        5.314694056926007,
        5.057077678380463
    };

    double p = 5e5;

    // Calculate density using Peng-Robinson EoS from Cantera
    for(int i=0; i<10; i++)
    {
        const double temp = 300 + i*25;
        set_r(1.0);
        test_phase->setState_TP(temp, p);
        EXPECT_NEAR(test_phase->density(),rhoCoolProp[i],1.e-5);
    }
}

TEST_F(PengRobinson_Test, partialMolarPropertyIdentities)
{
    // unique_ptr<ThermoPhase> phase(newPhase("co2_PR_example.yaml"));
    vector_fp hk(test_phase->nSpecies());
    vector_fp uk(test_phase->nSpecies());
    vector_fp sk(test_phase->nSpecies());
    vector_fp gk(test_phase->nSpecies());
    vector_fp vk(test_phase->nSpecies());
    vector_fp X(test_phase->nSpecies());

    test_phase->setMoleFractionsByName("CO2: 0.7, H2O: 0.1, H2: 0.2");
    test_phase->getMoleFractions(X.data());
    double P = 100 * OneAtm;

    for (int i = 0; i < 5; i++) {
        double T = 300 + i*60;
        test_phase->setState_TP(T, P);
        test_phase->getPartialMolarEnthalpies(hk.data());
        test_phase->getPartialMolarIntEnergies(uk.data());
        test_phase->getPartialMolarEntropies(sk.data());
        test_phase->getChemPotentials(gk.data());
        test_phase->getPartialMolarVolumes(vk.data());

        double h_mix = test_phase->enthalpy_mole();
        double u_mix = test_phase->intEnergy_mole();
        double s_mix = test_phase->entropy_mole();
        double g_mix = test_phase->gibbs_mole();
        double v_mix = test_phase->molarVolume();

        double h = dot(X.begin(), X.end(), hk.begin());
        EXPECT_NEAR(h, h_mix, 1e-11 * std::abs(h_mix));
        double u = dot(X.begin(), X.end(), uk.begin());
        EXPECT_NEAR(u, u_mix, 1e-11 * std::abs(u_mix));
        double s = dot(X.begin(), X.end(), sk.begin());
        EXPECT_NEAR(s, s_mix, 1e-11 * std::abs(s_mix));
        double g = dot(X.begin(), X.end(), gk.begin());
        EXPECT_NEAR(g, g_mix, 1e-11 * std::abs(g_mix));
        double v = dot(X.begin(), X.end(), vk.begin());
        EXPECT_NEAR(v, v_mix, 1e-11 * std::abs(v_mix));

        for (size_t k = 0; k < test_phase->nSpecies(); k++) {
            EXPECT_NEAR(uk[k] + P * vk[k], hk[k], 1e-11 * std::abs(h_mix));
            EXPECT_NEAR(hk[k] - T * sk[k], gk[k], 1e-11 * std::abs(g_mix));
        }
    }
}

TEST(PengRobinson, lookupSpeciesProperties)
{
    AnyMap phase_def = AnyMap::fromYamlString(
        "{name: test, species: [{gri30.yaml/species: [CO2, CH4, N2]}],"
        " thermo: Peng-Robinson}"
    );
    unique_ptr<ThermoPhase> test(newPhase(phase_def));

    // Check for correspondence to values in critical properties "database"
    test->setState_TPX(330, 100 * OneAtm, "CH4: 1.0");
    EXPECT_NEAR(test->critPressure(), 4.63e+06, 1e-2);
    EXPECT_NEAR(test->critTemperature(), 190.7, 1e-3);

    test->setState_TPX(330, 180 * OneAtm, "N2: 1.0");
    EXPECT_NEAR(test->critPressure(), 3.39e+06, 1e-2);
    EXPECT_NEAR(test->critTemperature(), 126.2, 1e-3);
}

TEST(PengRobinson, lookupSpeciesPropertiesMissing)
{
    AnyMap phase_def = AnyMap::fromYamlString(
        "{name: test, species: [{gri30.yaml/species: [CO2, CH3, N2]}],"
        " thermo: Peng-Robinson}"
    );

    // CH3 is not in the critical properties database, so this should be
    // detected as an error
    EXPECT_THROW(newPhase(phase_def), CanteraError);
}

TEST(PengRobinson, localCritProperties)
{
    // Test calculation based on critical properties stored in the YAML species
    // definition, in the "critical-parameters" field
    unique_ptr<ThermoPhase> test(newPhase("thermo-models.yaml", "CO2-PR-params"));
    test->setState_TPX(400, 1.2e6, "CO2: 1.0");
    EXPECT_NEAR(test->critTemperature(), 304.128, 1e-5);
    EXPECT_NEAR(test->critPressure(), 7.3773e6, 1e-4);

    test->setState_TPX(400, 1.2e6, "H2O: 1.0");
    EXPECT_NEAR(test->critTemperature(), 647.096, 1e-5);
    EXPECT_NEAR(test->critPressure(), 22.064e6, 1e-4);
}

};
