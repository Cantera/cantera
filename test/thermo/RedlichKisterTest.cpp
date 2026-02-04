#include "gtest/gtest.h"
#include "cantera/thermo/RedlichKisterVPSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/PDSS_ConstVol.h"

#include <array>

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

    void setup() {
        test_phase = newThermo("thermo-models.yaml", "Redlich-Kister-LiC6");
    }

    void set_r(const double r) {
        vector<double> moleFracs(2);
        moleFracs[0] = r;
        moleFracs[1] = 1-r;
        test_phase->setMoleFractions(moleFracs);
    }

    shared_ptr<ThermoPhase> test_phase;
};

TEST_F(RedlichKister_Test, chem_potentials)
{
    setup();
    test_phase->setState_TP(298.15, 101325.);

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector<double> chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        set_r(xmin + i*dx);
        test_phase->getChemPotentials(chemPotentials);
        EXPECT_NEAR(expected_chempot[i], chemPotentials[0], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, dlnActivities)
{
    setup();
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
    vector<double> dlnActCoeffdx(2);
    for(int i=0; i < 9; ++i)
    {
        const double r = xmin + i*dx;
        set_r(r);
        test_phase->getdlnActCoeffdlnX_diag(dlnActCoeffdx);
        EXPECT_NEAR(expected_result[i], dlnActCoeffdx[0], 1.e-6);
    }
}

TEST_F(RedlichKister_Test, standardConcentrations)
{
    setup();
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(1));
}

TEST_F(RedlichKister_Test, fromScratch)
{
    test_phase = make_shared<RedlichKisterVPSSTP>();
    RedlichKisterVPSSTP& rk = dynamic_cast<RedlichKisterVPSSTP&>(*test_phase);

    auto sLiC6 = make_shared<Species>("Li(C6)", parseCompString("C:6 Li:1"));
    double coeffs1[] = {298.15, -11.65e6, 0.0, 0.0};
    sLiC6->thermo = make_shared<ConstCpPoly>(100, 5000, 101325, coeffs1);

    auto sVC6 = make_shared<Species>("V(C6)", parseCompString("C:6"));
    double coeffs2[] = {298.15, 0.0, 0.0, 0.0};
    sVC6->thermo = make_shared<ConstCpPoly>(250, 800, 101325, coeffs2);

    rk.addSpecies(sLiC6);
    rk.addSpecies(sVC6);

    auto ssLiC6 = make_unique<PDSS_ConstVol>();
    rk.installPDSS(0, std::move(ssLiC6));

    auto ssVC6 = make_unique<PDSS_ConstVol>();
    rk.installPDSS(1, std::move(ssVC6));

    std::array hcoeffs {-3.268E6, 3.955E6, -4.573E6, 6.147E6, -3.339E6, 1.117E7,
                        2.997E5, -4.866E7, 1.362E5, 1.373E8, -2.129E7, -1.722E8,
                        3.956E7, 9.302E7, -3.280E7};
    std::array scoeffs {0.0};

    rk.addBinaryInteraction("Li(C6)", "V(C6)", hcoeffs, scoeffs);
    rk.initThermo();
    rk.setState_TPX(298.15, 101325, "Li(C6):0.6,V(C6):0.4");

    double xmin = 0.6;
    double xmax = 0.9;
    int numSteps = 9;
    vector<double> chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        set_r(xmin + i*(xmax-xmin)/(numSteps-1));
        test_phase->getChemPotentials(chemPotentials);
        EXPECT_NEAR(expected_chempot[i], chemPotentials[0], 1.e-6);
    }
}

};
