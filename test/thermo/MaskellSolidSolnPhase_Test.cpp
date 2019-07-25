#include "gtest/gtest.h"
#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/base/stringUtils.h"
#include <iostream>

namespace Cantera
{

static const double expected_result_0[9] = {1.2338461937651645e7, 8.0117750485066967e6, 4.9909899513507355e6, 2.4159732793453853e6, 0., -2.415973279345389e6, -4.9909899513507383e6, -8.0117750485066995e6, -1.2338461937651645e7};
static const double expected_result_5000[9] = { 1.2336254543579847e7, 8.0099571647402933e6, 4.9896777890600581e6, 2.4152804151720102e6, 0., -2.4152804151720135e6, -4.9896777890600609e6, -8.0099571647402961e6, -1.2336254543579847e7 };
static const double expected_result_minus_5000[9] = { 1.2340671804814456e7, 8.0135951995088588e6, 4.9923039182151733e6, 2.4166671660301457e6, 0., -2.4166671660301457e6, -4.9923039182151733e6, -8.0135951995088588e6, -1.2340671804814456e7};

class MaskellSolidSolnPhase_Test : public testing::Test
{
public:
    void initializeTestPhaseWithXML(const std::string & filename)
    {
        test_phase.reset(newPhase(filename));
    }

    void set_r(const double r) {
        vector_fp moleFracs(2);
        moleFracs[0] = r;
        moleFracs[1] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    void check_chemPotentials(const double expected_result[9])
    {
        vector_fp chemPotentials(2);
        for(int i=0; i < 9; ++i)
        {
            const double r = 0.1 * (i+1);
            set_r(r);
            test_phase->getChemPotentials(&chemPotentials[0]);
            EXPECT_NEAR(-expected_result[i], chemPotentials[0], 1.e-6);
            EXPECT_NEAR(1000.+expected_result[i], chemPotentials[1], 1.e-6);
        }
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(MaskellSolidSolnPhase_Test, construct_from_xml)
{
    const std::string invalid_file("../data/MaskellSolidSolnPhase_nohmix.xml");
    EXPECT_THROW(initializeTestPhaseWithXML(invalid_file), CanteraError);

    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    MaskellSolidSolnPhase* maskell_phase = dynamic_cast<MaskellSolidSolnPhase*>(test_phase.get());
    EXPECT_TRUE(maskell_phase != NULL);
}

TEST_F(MaskellSolidSolnPhase_Test, chem_potentials)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    test_phase->setState_TP(298., 1.);
    set_r(0.5);

    MaskellSolidSolnPhase* maskell_phase = dynamic_cast<MaskellSolidSolnPhase*>(test_phase.get());

    maskell_phase->set_h_mix(0.);
    check_chemPotentials(expected_result_0);

    maskell_phase->set_h_mix(5000.);
    check_chemPotentials(expected_result_5000);

    maskell_phase->set_h_mix(-5000.);
    check_chemPotentials(expected_result_minus_5000);
}

TEST_F(MaskellSolidSolnPhase_Test, partialMolarVolumes)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);

    vector_fp pmv(2);
    test_phase->getPartialMolarVolumes(&pmv[0]);
    EXPECT_EQ(0.005, pmv[0]);
    EXPECT_EQ(0.01, pmv[1]);
}

TEST_F(MaskellSolidSolnPhase_Test, activityCoeffs)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    test_phase->setState_TP(298., 1.);
    set_r(0.5);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * 298.;
    vector_fp mu0(2);
    vector_fp activityCoeffs(2);
    vector_fp chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        const double r = 0.1 * (i+1);
        set_r(r);
        test_phase->getChemPotentials(&chemPotentials[0]);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        test_phase->getStandardChemPotentials(&mu0[0]);
        EXPECT_NEAR(chemPotentials[0], mu0[0] + RT*std::log(activityCoeffs[0] * r), 1.e-6);
        EXPECT_NEAR(chemPotentials[1], mu0[1] + RT*std::log(activityCoeffs[1] * (1-r)), 1.e-6);
    }
}

TEST_F(MaskellSolidSolnPhase_Test, standardConcentrations)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);

    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(1));
}

TEST_F(MaskellSolidSolnPhase_Test, activityConcentrations)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);

    // Check to make sure activityConcentration_i == standardConcentration_i * gamma_i * X_i
    vector_fp standardConcs(2);
    vector_fp activityCoeffs(2);
    vector_fp activityConcentrations(2);
    for(int i=0; i < 9; ++i)
    {
        const double r = 0.1 * (i+1);
        set_r(r);
        test_phase->getActivityCoefficients(&activityCoeffs[0]);
        standardConcs[0] = test_phase->standardConcentration(0);
        standardConcs[1] = test_phase->standardConcentration(1);
        test_phase->getActivityConcentrations(&activityConcentrations[0]);

        EXPECT_NEAR(standardConcs[0] * r * activityCoeffs[0], activityConcentrations[0], 1.e-6);
        EXPECT_NEAR(standardConcs[1] * (1-r) * activityCoeffs[1], activityConcentrations[1], 1.e-6);
    }
}

TEST_F(MaskellSolidSolnPhase_Test, fromScratch) {
    auto sH = make_shared<Species>("H(s)", parseCompString("H:1 He:2"));
    double coeffs1[] = {1.0, 0.0, 0.0, 0.0};
    sH->thermo.reset(new ConstCpPoly(250, 800, 1e5, coeffs1));

    auto sHe = make_shared<Species>("He(s)", parseCompString("He:1"));
    double coeffs2[] = {1.0, 1000.0, 0.0, 0.0};
    sHe->thermo.reset(new ConstCpPoly(250, 800, 1e5, coeffs2));

    MaskellSolidSolnPhase* p = new MaskellSolidSolnPhase();
    test_phase.reset(p);

    p->addSpecies(sH);
    p->addSpecies(sHe);

    std::unique_ptr<PDSS_ConstVol> ssH(new PDSS_ConstVol());
    ssH->setMolarVolume(0.005);
    p->installPDSS(0, std::move(ssH));

    std::unique_ptr<PDSS_ConstVol> ssHe(new PDSS_ConstVol());
    ssHe->setMolarVolume(0.01);
    p->installPDSS(1, std::move(ssHe));

    p->set_h_mix(5000);
    p->setProductSpecies("H(s)");

    p->initThermo();
    p->setState_TPX(298, 1, "H(s):0.90   He(s):0.10");

    vector_fp pmv(2);
    p->getPartialMolarVolumes(&pmv[0]);
    EXPECT_EQ(0.005, pmv[0]);
    EXPECT_EQ(0.01, pmv[1]);

    // Compare with XML chem_potentials test
    p->set_h_mix(5000.);
    check_chemPotentials(expected_result_5000);

    p->set_h_mix(-5000.);
    check_chemPotentials(expected_result_minus_5000);
}

};
