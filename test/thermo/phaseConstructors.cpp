#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ShomatePoly.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include <fstream>
#include "thermo_data.h"

namespace Cantera
{

class FixedChemPotSstpConstructorTest : public testing::Test
{
};

TEST_F(FixedChemPotSstpConstructorTest, fromXML)
{
    std::unique_ptr<ThermoPhase> p(newPhase("../data/LiFixed.xml"));
    ASSERT_EQ((int) p->nSpecies(), 1);
    double mu;
    p->getChemPotentials(&mu);
    ASSERT_DOUBLE_EQ(-2.3e7, mu);
}

TEST_F(FixedChemPotSstpConstructorTest, SimpleConstructor)
{
    FixedChemPotSSTP p("Li", -2.3e7);
    ASSERT_EQ((int) p.nSpecies(), 1);
    double mu;
    p.getChemPotentials(&mu);
    ASSERT_DOUBLE_EQ(-2.3e7, mu);
}

#ifndef HAS_NO_PYTHON // skip these tests if the Python converter is unavailable
class CtiConversionTest : public testing::Test
{
public:
    CtiConversionTest() {
        appdelete();
    }

    std::unique_ptr<ThermoPhase> p1;
    std::unique_ptr<ThermoPhase> p2;
    void compare()
    {
        ASSERT_EQ(p1->nSpecies(), p2->nSpecies());
        for (size_t i = 0; i < p1->nSpecies(); i++) {
            ASSERT_EQ(p1->speciesName(i), p2->speciesName(i));
            ASSERT_EQ(p1->molecularWeight(i), p2->molecularWeight(i));
        }
    }
};

TEST_F(CtiConversionTest, ExplicitConversion) {
    p1.reset(newPhase("../data/air-no-reactions.xml"));
    ct2ctml("../data/air-no-reactions.cti");
    p2.reset(newPhase("air-no-reactions.xml", ""));
    compare();
}

TEST_F(CtiConversionTest, ImplicitConversion) {
    p1.reset(newPhase("../data/air-no-reactions.xml"));
    p2.reset(newPhase("../data/air-no-reactions.cti"));
    compare();
}

class ChemkinConversionTest : public testing::Test {
public:
    void copyInputFile(const std::string& name) {
        std::string in_name = "../data/" + name;
        std::ifstream source(in_name, std::ios::binary);
        std::ofstream dest(name, std::ios::binary);
        dest << source.rdbuf();
    }
};

TEST_F(ChemkinConversionTest, ValidConversion) {
    copyInputFile("pdep-test.inp");
    ck2cti("pdep-test.inp");
    std::unique_ptr<ThermoPhase> p(newPhase("pdep-test.cti"));
    ASSERT_GT(p->temperature(), 0.0);
}

TEST_F(ChemkinConversionTest, MissingInputFile) {
    ASSERT_THROW(ck2cti("nonexistent-file.inp"),
                 CanteraError);
}

TEST_F(ChemkinConversionTest, FailedConversion) {
    copyInputFile("h2o2_missingThermo.inp");
    ASSERT_THROW(ck2cti("h2o2_missingThermo.inp"),
                 CanteraError);
}
#endif

class ConstructFromScratch : public testing::Test
{
public:
    ConstructFromScratch()
        : sH2O(new Species("H2O", parseCompString("H:2 O:1")))
        , sH2(new Species("H2", parseCompString("H:2")))
        , sO2(new Species("O2", parseCompString("O:2")))
        , sOH(new Species("OH", parseCompString("H:1 O:1")))
        , sCO(new Species("CO", parseCompString("C:1 O:1")))
        , sCO2(new Species("CO2", parseCompString("C:1 O:2")))
    {
        sH2O->thermo.reset(new NasaPoly2(200, 3500, 101325, h2o_nasa_coeffs));
        sH2->thermo.reset(new NasaPoly2(200, 3500, 101325, h2_nasa_coeffs));
        sO2->thermo.reset(new NasaPoly2(200, 3500, 101325, o2_nasa_coeffs));
        sOH->thermo.reset(new NasaPoly2(200, 3500, 101325, oh_nasa_coeffs));
        sCO->thermo.reset(new NasaPoly2(200, 3500, 101325, o2_nasa_coeffs));
        sCO2->thermo.reset(new ShomatePoly2(200, 3500, 101325, co2_shomate_coeffs));
    }

    shared_ptr<Species> sH2O, sH2, sO2, sOH, sCO, sCO2;
};

TEST_F(ConstructFromScratch, AddElements)
{
    IdealGasPhase p;
    p.addElement("H");
    p.addElement("O");
    ASSERT_EQ((size_t) 2, p.nElements());
    ASSERT_EQ("H", p.elementName(0));
    ASSERT_EQ((size_t) 1, p.elementIndex("O"));
}

TEST_F(ConstructFromScratch, AddSpeciesDefaultBehavior)
{
    IdealGasPhase p;
    p.addElement("H");
    p.addElement("O");
    p.addSpecies(sH2O);
    p.addSpecies(sH2);

    ASSERT_EQ((size_t) 2, p.nSpecies());

    p.addSpecies(sO2);
    p.addSpecies(sOH);

    ASSERT_EQ((size_t) 4, p.nSpecies());
    ASSERT_EQ("H2", p.speciesName(1));
    ASSERT_EQ(2, p.nAtoms(2, 1)); // O in O2
    ASSERT_EQ(2, p.nAtoms(0, 0)); // H in H2O
    ASSERT_THROW(p.addSpecies(sCO), CanteraError);
}

TEST_F(ConstructFromScratch, ignoreUndefinedElements)
{
    IdealGasPhase p;
    p.addElement("H");
    p.addElement("O");
    p.ignoreUndefinedElements();

    p.addSpecies(sO2);
    p.addSpecies(sOH);
    ASSERT_EQ((size_t) 2, p.nSpecies());

    p.addSpecies(sCO);
    p.addSpecies(sCO2);
    ASSERT_EQ((size_t) 2, p.nSpecies());
    ASSERT_EQ((size_t) 2, p.nElements());
    ASSERT_EQ(npos, p.speciesIndex("CO2"));
}

TEST_F(ConstructFromScratch, addUndefinedElements)
{
    IdealGasPhase p;
    p.addElement("H");
    p.addElement("O");
    p.addUndefinedElements();

    p.addSpecies(sH2);
    p.addSpecies(sOH);
    ASSERT_EQ((size_t) 2, p.nSpecies());
    ASSERT_EQ((size_t) 2, p.nElements());

    p.addSpecies(sCO);
    p.addSpecies(sCO2);
    ASSERT_EQ((size_t) 4, p.nSpecies());
    ASSERT_EQ((size_t) 3, p.nElements());
    ASSERT_EQ((size_t) 1, p.nAtoms(p.speciesIndex("CO2"), p.elementIndex("C")));
    ASSERT_EQ((size_t) 2, p.nAtoms(p.speciesIndex("co2"), p.elementIndex("O")));
    p.setMassFractionsByName("H2:0.5, CO2:0.5");
    ASSERT_DOUBLE_EQ(0.5, p.massFraction("CO2"));
}

TEST_F(ConstructFromScratch, RedlichKwongMFTP)
{
    RedlichKwongMFTP p;
    p.addUndefinedElements();
    p.addSpecies(sCO2);
    p.addSpecies(sH2O);
    p.addSpecies(sH2);
    double fa = toSI("bar-cm6/mol2");
    double fb = toSI("cm3/mol");
    p.setBinaryCoeffs("H2", "H2O", 4 * fa, 40 * fa);
    p.setSpeciesCoeffs("CO2", 7.54e7 * fa, -4.13e4 * fa, 27.80 * fb);
    p.setBinaryCoeffs("CO2", "H2O", 7.897e7 * fa, 0.0);
    p.setSpeciesCoeffs("H2O", 1.7458e8 * fa, -8e4 * fa, 18.18 * fb);
    p.setSpeciesCoeffs("H2", 30e7 * fa, -330e4 * fa, 31 * fb);
    p.initThermo();
    p.setMoleFractionsByName("CO2:0.9998, H2O:0.0002");
    p.setState_TP(300, 200 * OneAtm);
    EXPECT_NEAR(p.pressure(), 200 * OneAtm, 1e-5);
    // Arbitrary regression test values
    EXPECT_NEAR(p.density(), 892.421, 2e-3);
    EXPECT_NEAR(p.enthalpy_mole(), -404848642.3797, 1e-3);
}

TEST(PureFluidFromScratch, CarbonDioxide)
{
    PureFluidPhase p;
    auto sCO2 = make_shared<Species>("CO2", parseCompString("C:1 O:2"));
    sCO2->thermo.reset(new ShomatePoly2(200, 6000, 101325, co2_shomate_coeffs));
    p.addUndefinedElements();
    p.addSpecies(sCO2);
    p.setSubstance("carbondioxide");
    p.initThermo();
    p.setState_Tsat(280, 0.5);
    EXPECT_NEAR(p.pressure(), 4160236.987, 1e-2);
}

TEST(WaterSSTP, fromScratch)
{
    WaterSSTP water;
    auto sH2O = make_shared<Species>("H2O", parseCompString("H:2 O:1"));
    sH2O->thermo.reset(new NasaPoly2(200, 3500, 101325, h2o_nasa_coeffs)); // unused
    water.addUndefinedElements();
    water.addSpecies(sH2O);
    water.initThermo();
    water.setState_TP(298.15, 1e5);
    EXPECT_NEAR(water.enthalpy_mole() / 1e6, -285.83, 2e-2);
}

} // namespace Cantera
