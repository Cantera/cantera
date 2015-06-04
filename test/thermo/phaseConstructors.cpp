#include "gtest/gtest.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/NasaPoly2.h"
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
    ThermoPhase* p = newPhase("../data/LiFixed.xml");
    ASSERT_EQ((int) p->nSpecies(), 1);
    double mu;
    p->getChemPotentials(&mu);
    ASSERT_DOUBLE_EQ(-2.3e7, mu);
    delete p;
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
    CtiConversionTest() : p1(0), p2(0) {
        appdelete();
    }
    ~CtiConversionTest() {
        delete p1;
        delete p2;
    }

    ThermoPhase* p1;
    ThermoPhase* p2;
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
    p1 = newPhase("../data/air-no-reactions.xml");
    ct2ctml("../data/air-no-reactions.cti");
    p2 = newPhase("air-no-reactions.xml", "");
    compare();
}

TEST_F(CtiConversionTest, ImplicitConversion) {
    p1 = newPhase("../data/air-no-reactions.xml");
    p2 = newPhase("../data/air-no-reactions.cti");
    compare();
}

class ChemkinConversionTest : public testing::Test {
public:
    void copyInputFile(const std::string& name) {
        std::string in_name = "../data/" + name;
        std::ifstream source(in_name.c_str(), std::ios::binary);
        std::ofstream dest(name.c_str(), std::ios::binary);
        dest << source.rdbuf();
    }
};

TEST_F(ChemkinConversionTest, ValidConversion) {
    copyInputFile("pdep-test.inp");
    ck2cti("pdep-test.inp");
    ThermoPhase* p = newPhase("pdep-test.cti");
    ASSERT_GT(p->temperature(), 0.0);
    delete p;
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
        sCO2->thermo.reset(new NasaPoly2(200, 3500, 101325, h2o_nasa_coeffs));
    }

    IdealGasPhase p;
    shared_ptr<Species> sH2O, sH2, sO2, sOH, sCO, sCO2;
};

TEST_F(ConstructFromScratch, AddElements)
{
    p.addElement("H");
    p.addElement("O");
    ASSERT_EQ((size_t) 2, p.nElements());
    ASSERT_EQ("H", p.elementName(0));
    ASSERT_EQ((size_t) 1, p.elementIndex("O"));
}

TEST_F(ConstructFromScratch, AddSpeciesDefaultBehavior)
{
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
    ASSERT_EQ((size_t) 2, p.nAtoms(p.speciesIndex("CO2"), p.elementIndex("O")));
    p.setMassFractionsByName("H2:0.5, CO2:0.5");
    ASSERT_DOUBLE_EQ(0.5, p.massFraction("CO2"));
}

} // namespace Cantera
