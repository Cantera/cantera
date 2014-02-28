#include "gtest/gtest.h"
#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/thermo/SimpleThermo.h"
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

    void initializeTestPhaseWithSimpleThermo()
    {
        test_phase = new MaskellSolidSolnPhase();
        test_phase->addElement("A", 1.);
        test_phase->addElement("B", 2.);
        std::vector<double> comp(2);
        comp[0] = 1.;
        comp[1] = 0.;
        test_phase->addSpecies("A", &comp[0], 0., 1.);
        comp[0] = 0.;
        comp[1] = 1.;
        test_phase->addSpecies("B", &comp[0], 0., 1.);

        // Setup simple thermo  so that the standard state enthalpy and
        // gibbs free energies are always 0 so that we can just test the
        // additional contribution from the Maskell model
        SimpleThermo * spec_thermo = new SimpleThermo();
        std::vector<double> coeffs(4);
        coeffs[0] = 1;
        coeffs[1] = 0;
        coeffs[2] = 0;
        coeffs[3] = 0;
        spec_thermo->install("A", 0, 0, &coeffs[0], 0., 1000., 1.);
        coeffs[1] = 1000;
        spec_thermo->install("B", 1, 0, &coeffs[0], 0., 1000., 1.);
        test_phase->setSpeciesThermo(spec_thermo);

        test_phase->setState_TP(298., 1.);
        set_r(0.5);
    }

    void initializeTestPhaseWithXML(const std::string & filename)
    {
        test_phase = newPhase(filename.c_str(), "");
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
            EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
            EXPECT_NEAR(1000.-expected_result[i], chemPotentials[1], 1.e-6);
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
}

TEST_F(MaskellSolidSolnPhase_Test, chem_potentials)
{
    initializeTestPhaseWithSimpleThermo();

    MaskellSolidSolnPhase * maskell_phase = dynamic_cast<MaskellSolidSolnPhase *>(test_phase);

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

};
