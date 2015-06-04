#include "gtest/gtest.h"
#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/thermo/SimpleThermo.h"
#include "cantera/thermo/VPSSMgr_General.h"
#include "cantera/thermo/ThermoFactory.h"
#include <iostream>

namespace Cantera
{

class MaskellSolidSolnPhase_Test : public testing::Test
{
protected:
    ThermoPhase *test_phase;
public:
    MaskellSolidSolnPhase_Test() : test_phase(NULL) {}

    ~MaskellSolidSolnPhase_Test() { delete test_phase; }

    void initializeTestPhaseWithXML(const std::string & filename)
    {
        test_phase = newPhase(filename.c_str());
    }

    void set_r(const double r) {
        std::vector<double> moleFracs(2);
        moleFracs[0] = r;
        moleFracs[1] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    void check_chemPotentials(const double expected_result[9])
    {
        std::vector<double> chemPotentials(2);
        for(int i=0; i < 9; ++i)
        {
            const double r = 0.1 * (i+1);
            set_r(r);
            test_phase->getChemPotentials(&chemPotentials[0]);
            EXPECT_NEAR(-expected_result[i], chemPotentials[0], 1.e-6);
            EXPECT_NEAR(1000.+expected_result[i], chemPotentials[1], 1.e-6);
        }
    }
};

TEST_F(MaskellSolidSolnPhase_Test, construct_from_xml)
{
    const std::string invalid_file("../data/MaskellSolidSolnPhase_nohmix.xml");
    EXPECT_THROW(initializeTestPhaseWithXML(invalid_file), CanteraError);
    delete test_phase;

    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    MaskellSolidSolnPhase * maskell_phase = dynamic_cast<MaskellSolidSolnPhase *>(test_phase);
    EXPECT_TRUE(maskell_phase != NULL);
}

TEST_F(MaskellSolidSolnPhase_Test, chem_potentials)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    test_phase->setState_TP(298., 1.);
    set_r(0.5);

    MaskellSolidSolnPhase * maskell_phase = dynamic_cast<MaskellSolidSolnPhase *>(test_phase);
    ASSERT_TRUE(maskell_phase != NULL);

    maskell_phase->set_h_mix(0.);
    const double expected_result_0[9] = {1.2338461168724738e7, 8.011774549216799e6, 4.990989640314685e6, 2.415973128783114e6, 0., -2.415973128783114e6, -4.99098964031469e6, -8.0117745492168e6, -1.2338461168724738e7};
    check_chemPotentials(expected_result_0);

    maskell_phase->set_h_mix(5000.);
    const double expected_result_5000[9] = { 1.233625377465302e7, 8.00995666545047e6, 4.989677478024063e6, 2.41528026460977e6, 0., -2.415280264609771e6, -4.989677478024068e6, -8.00995666545047e6, -1.233625377465302e7 };
    check_chemPotentials(expected_result_5000);

    maskell_phase->set_h_mix(-5000.);
    const double expected_result_minus_5000[9] = { 1.2340671035887627e7, 8.013594700219031e6, 4.992303607179179e6, 2.4166670154679064e6, 0., -2.4166670154679064e6, -4.9923036071791835e6, -8.013594700219034e6, -1.2340671035887627e7};
    check_chemPotentials(expected_result_minus_5000);
}

TEST_F(MaskellSolidSolnPhase_Test, partialMolarVolumes)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    ASSERT_TRUE(dynamic_cast<MaskellSolidSolnPhase *>(test_phase) != NULL);

    std::vector<double> pmv(2);
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

    MaskellSolidSolnPhase * maskell_phase = dynamic_cast<MaskellSolidSolnPhase *>(test_phase);
    ASSERT_TRUE(maskell_phase != NULL);

    // Test that mu0 + RT log(activityCoeff * MoleFrac) == mu
    const double RT = GasConstant * 298.;
    std::vector<double> mu0(2);
    std::vector<double> activityCoeffs(2);
    std::vector<double> chemPotentials(2);
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
    ASSERT_TRUE(dynamic_cast<MaskellSolidSolnPhase *>(test_phase) != NULL);

    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(0));
    EXPECT_DOUBLE_EQ(1.0, test_phase->standardConcentration(1));
}

TEST_F(MaskellSolidSolnPhase_Test, activityConcentrations)
{
    const std::string valid_file("../data/MaskellSolidSolnPhase_valid.xml");
    initializeTestPhaseWithXML(valid_file);
    ASSERT_TRUE(dynamic_cast<MaskellSolidSolnPhase *>(test_phase) != NULL);

    // Check to make sure activityConcentration_i == standardConcentration_i * gamma_i * X_i
    std::vector<double> standardConcs(2);
    std::vector<double> activityCoeffs(2);
    std::vector<double> activityConcentrations(2);
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

};
