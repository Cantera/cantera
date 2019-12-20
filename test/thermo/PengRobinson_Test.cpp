#include "gtest/gtest.h"
#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"


namespace Cantera
{

class PengRobinson_Test : public testing::Test
{
public:
    PengRobinson_Test() {
<<<<<<< HEAD
        test_phase.reset(newPhase("../data/co2_PR_example.yaml"));
=======
        test_phase.reset(newPhase("../data/co2_PR_example.cti"));
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
TEST_F(PengRobinson_Test, construct_from_yaml)
=======
TEST_F(PengRobinson_Test, construct_from_cti)
>>>>>>> Removing MFTP suffix from PengRobinson class
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
<<<<<<< HEAD
        -457361607.71983075,
        -457350560.54839599,
        -457340699.25698096,
        -457331842.5539279,
        -457323844.32100844,
        -457316585.4752928,
        -457309967.99120748,
        -457303910.44199038,
        -457298344.62820804
=======
        -4.5736182681761962e+008,
        -4.5733771904416579e+008,
        -4.5732943831449223e+008,
        -4.5732206687414169e+008,
        -4.5731546826955432e+008,
        -4.5730953161186475e+008,
        -4.5730416590547645e+008,
        -4.5729929581635743e+008,
        -4.5729485847173005e+008
>>>>>>> Removing MFTP suffix from PengRobinson class
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

<<<<<<< HEAD
TEST_F(PengRobinson_Test, chemPotentials_RT)
{
    test_phase->setState_TP(298., 1.);

    // Test that chemPotentials_RT*RT = chemPotentials
    const double RT = GasConstant * 298.;
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

=======
>>>>>>> Removing MFTP suffix from PengRobinson class
TEST_F(PengRobinson_Test, activityCoeffs)
{
    test_phase->setState_TP(298., 1.);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * 298.;
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
    EXPECT_DOUBLE_EQ(test_phase->pressure()/(test_phase->temperature()*GasConstant), test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(test_phase->pressure()/(test_phase->temperature()*GasConstant), test_phase->standardConcentration(1));
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
    const double p1[6] = {
<<<<<<< HEAD
        1.7084253549322079e+002,
        1.6543121742659784e+002,
        1.6066148681014121e+002,
        1.5639259178086871e+002,
        1.525268502259365e+002,
        1.4899341020317422e+002
    };
    // Phase change between temperatures 4 & 5:
    const double p2[6] = {
        7.3025772038910179e+002,
        7.0307777665949902e+002,
        6.7179381832541878e+002,
        6.3389192023192868e+002,
        5.8250166044528487e+002,
        3.8226318921022073e+002
    };
    // Supercritical; no discontinuity in rho values:
    const double p3[6] = {
        7.8626284773748239e+002,
        7.6743023689871097e+002,
        7.4747463917603955e+002,
        7.2620055080831412e+002,
        7.0335270498118734e+002,
        6.7859003092723128e+002
=======
        1.7474528924963985e+002,
        1.6800540828415956e+002,
        1.62278413743154e+002,
        1.5728963799103039e+002,
        1.5286573762819748e+002,
        1.4888956030449546e+002
    };
    // Phase change between temperatures 4 & 5:
    const double p2[6] = {
        7.5565889855724288e+002,
        7.2577747673480337e+002,
        6.913183942651284e+002,
        6.494661249672663e+002,
        5.9240469307757724e+002,
        3.645826047440932e+002
    };
    // Supercritical; no discontinuity in rho values:
    const double p3[6] = {
        8.047430802847415e+002,
        7.8291565113595595e+002,
        7.5958477920749681e+002,
        7.3445460137134626e+002,
        7.0712433093853724e+002,
        6.77034438769492e+002
>>>>>>> Removing MFTP suffix from PengRobinson class
    };

    for(int i=0; i<6; ++i)
    {
        const double temp = 294 + i*2;
        set_r(0.999);
        test_phase->setState_TP(temp, 5542027.5);
        EXPECT_NEAR(test_phase->density(),p1[i],1.e-8);

<<<<<<< HEAD
        test_phase->setState_TP(temp, 7388370.);
=======
        test_phase->setState_TP(temp, 7389370.);
>>>>>>> Removing MFTP suffix from PengRobinson class
        EXPECT_NEAR(test_phase->density(),p2[i],1.e-8);

        test_phase->setState_TP(temp, 9236712.5);
        EXPECT_NEAR(test_phase->density(),p3[i],1.e-8);
    }
}

TEST_F(PengRobinson_Test, getPressure)
{
    // Check to make sure that the P-R equation is accurately reproduced for a few selected values

<<<<<<< HEAD
    /* This test uses CO2 as the only species (mole fraction 99.9%, balance H2).
=======
    /* This test uses CO2 as the only species.
>>>>>>> Removing MFTP suffix from PengRobinson class
    *  Values of a_coeff, b_coeff are calculated based on the the critical temperature and pressure values of CO2 as follows:
    *       a_coeff = 0.457235(RT_crit)^2/p_crit
    *       b_coeff = 0.077796(RT_crit)/p_crit
    *  The temperature dependent parameter in P-R EoS is calculated as
    *       \alpha = [1 + \kappa(1 - sqrt{T/T_crit}]^2
    *  kappa is a function calulated based on the accentric factor.
    */

<<<<<<< HEAD
    double a_coeff = 3.958095109E+5;
    double b_coeff = 26.62616317/1000;
    double acc_factor = 0.228;
    double pres_theoretical, kappa, alpha, mv;
    const double rho = 1.0737; 
=======
    double a_coeff = 3.958134E+5;
    double b_coeff = 26.6275/1000;
    double acc_factor = 0.228;
    double pres_theoretical, kappa, alpha, mv;
    const double rho = 1.0737; 
    const double Tcrit = test_phase->critTemperature();
>>>>>>> Removing MFTP suffix from PengRobinson class

    //Calculate kappa value 
    kappa = 0.37464 + 1.54226*acc_factor - 0.26992*acc_factor*acc_factor;

    for (int i = 0; i<10; i++)
    {
        const double temp = 296 + i * 2;
<<<<<<< HEAD
        set_r(1.0);        
        test_phase->setState_TR(temp, rho);
        const double Tcrit = test_phase->critTemperature();
        mv = 1 / rho * test_phase->meanMolecularWeight();
        //Calculate pressure using Peng-Robinson EoS
        alpha = pow(1 + kappa*(1 - sqrt(temp / Tcrit)), 2.0);
        pres_theoretical = GasConstant*temp / (mv - b_coeff) - a_coeff*alpha / (mv*mv + 2*b_coeff*mv - b_coeff*b_coeff);
        EXPECT_NEAR(test_phase->pressure(), pres_theoretical, 3);
=======
        set_r(0.999);        
        test_phase->setState_TR(temp, 1.0737);
        mv = 1 / rho * test_phase->meanMolecularWeight();
        //Calculate pressure using Peng-Robinson EoS
        alpha = pow(1 + kappa*(1 - sqrt(temp / Tcrit)), 2);
        pres_theoretical = GasConstant*temp / (mv - b_coeff) - a_coeff*alpha / (mv*mv + 2*b_coeff*mv - b_coeff*b_coeff);
        EXPECT_NEAR(test_phase->pressure(), pres_theoretical, 2);
>>>>>>> Removing MFTP suffix from PengRobinson class
    }
}
};
