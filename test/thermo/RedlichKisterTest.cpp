#include "gtest/gtest.h"
#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/PDSS_IdealGas.h"

namespace Cantera
{

static const double expected_chempot[9] = {
    -1.2791500499152161e+007,
    -1.2618554573674981e+007,
    -1.2445418333486753e+007,
    -1.2282611731533309e+007,
    -1.2134110797552738e+007,
    -1.1999465433876401e+007,
    -1.1882669440244285e+007,
    -1.1792994862336770e+007,
    -1.1730896003312804e+007
};

class RedlichKister_Test : public testing::Test
{
public:
    RedlichKister_Test() {}

    void initXML() {
        test_phase.reset(newPhase("../data/RedlichKisterVPSSTP_valid.xml"));
    }

    void set_r(const double r) {
        vector_fp moleFracs(2);
        moleFracs[0] = r;
        moleFracs[1] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(RedlichKister_Test, construct_from_xml)
{
    initXML();
    RedlichKisterVPSSTP* redlich_kister_phase = dynamic_cast<RedlichKisterVPSSTP*>(test_phase.get());
    ASSERT_TRUE(redlich_kister_phase != NULL);
}

TEST_F(RedlichKister_Test, chem_potentials)
{
    initXML();
    test_phase->setState_TP(298.15, 101325.);

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        set_r(xmin + i*dx);
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_chempot[i], chemPotentials[0], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, dlnActivities)
{
    initXML();
    test_phase->setState_TP(298.15, 101325.);

    const double expected_result[9] = {
        0.0907127,
        0.200612,
        0.229316,
        0.193278,
        0.142257,
        0.0766133,
        -0.0712113,
        -0.309379,
        -0.492206
    };

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp dlnActCoeffdx(2);
    for(int i=0; i < 9; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getdlnActCoeffdlnX_diag(&dlnActCoeffdx[0]);
        EXPECT_NEAR(expected_result[i], dlnActCoeffdx[0], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, activityCoeffs)
{
    initXML();
    test_phase->setState_TP(298., 1.);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * 298.;
    vector_fp mu0(2);
    vector_fp activityCoeffs(2);
    vector_fp chemPotentials(2);
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
        EXPECT_NEAR(chemPotentials[1], mu0[1] + RT*std::log(activityCoeffs[1] * (1-r)), 1.e-6);
    }
}

TEST_F(RedlichKister_Test, standardConcentrations)
{
    initXML();
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(1));
}

TEST_F(RedlichKister_Test, activityConcentrations)
{
    initXML();
    // Check to make sure activityConcentration_i == standardConcentration_i * gamma_i * X_i
    vector_fp standardConcs(2);
    vector_fp activityCoeffs(2);
    vector_fp activityConcentrations(2);
    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);

    for(int i=0; i < 9; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        standardConcs[0] = test_phase->standardConcentration(0);
        standardConcs[1] = test_phase->standardConcentration(1);
        test_phase->getActivityConcentrations(&activityConcentrations[0]);

        EXPECT_NEAR(standardConcs[0] * r * activityCoeffs[0], activityConcentrations[0], 1.e-6);
        EXPECT_NEAR(standardConcs[1] * (1-r) * activityCoeffs[1], activityConcentrations[1], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, fromScratch)
{
    test_phase.reset(new RedlichKisterVPSSTP());
    RedlichKisterVPSSTP& rk = dynamic_cast<RedlichKisterVPSSTP&>(*test_phase);

    auto sLiC6 = make_shared<Species>("Li(C6)", parseCompString("C:6 Li:1"));
    double coeffs1[] = {298.15, -11.65e6, 0.0, 0.0};
    sLiC6->thermo.reset(new ConstCpPoly(100, 5000, 101325, coeffs1));

    auto sVC6 = make_shared<Species>("V(C6)", parseCompString("C:6"));
    double coeffs2[] = {298.15, 0.0, 0.0, 0.0};
    sVC6->thermo.reset(new ConstCpPoly(250, 800, 101325, coeffs2));

    rk.addSpecies(sLiC6);
    rk.addSpecies(sVC6);

    std::unique_ptr<PDSS> ssLiC6(new PDSS_IdealGas());
    rk.installPDSS(0, std::move(ssLiC6));

    std::unique_ptr<PDSS> ssVC6(new PDSS_IdealGas());
    rk.installPDSS(1, std::move(ssVC6));

    double hcoeffs[] = {-3.268E6, 3.955E6, -4.573E6, 6.147E6, -3.339E6, 1.117E7,
                        2.997E5, -4.866E7, 1.362E5, 1.373E8, -2.129E7, -1.722E8,
                        3.956E7, 9.302E7, -3.280E7};
    double scoeffs[] = {0.0};

    rk.addBinaryInteraction("Li(C6)", "V(C6)", hcoeffs, 15, scoeffs, 1);
    rk.initThermo();
    rk.setState_TPX(298.15, 101325, "Li(C6):0.6,V(C6):0.4");

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    vector_fp chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        set_r(xmin + i*(xmax-xmin)/(numSteps-1));
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_chempot[i], chemPotentials[0], 1.e-6);
    }
}

};
