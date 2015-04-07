#include "gtest/gtest.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/ctml.h"
#include <fstream>

namespace Cantera
{

class FixedChemPotSstpConstructorTest : public testing::Test
{
};

TEST_F(FixedChemPotSstpConstructorTest, fromXML)
{
    ThermoPhase* p = newPhase("../data/LiFixed.xml", "");
    ASSERT_EQ((int) p->nSpecies(), 1);
    double mu;
    p->getChemPotentials(&mu);
    ASSERT_FLOAT_EQ(-2.3e7, mu);
    delete p;
}

TEST_F(FixedChemPotSstpConstructorTest, SimpleConstructor)
{
    FixedChemPotSSTP p("Li", -2.3e7);
    ASSERT_EQ((int) p.nSpecies(), 1);
    double mu;
    p.getChemPotentials(&mu);
    ASSERT_FLOAT_EQ(-2.3e7, mu);
}

#ifndef HAS_NO_PYTHON // skip these tests if the Python converter is unavailable
class CtiConversionTest : public testing::Test
{
public:
    CtiConversionTest() {
        appdelete();
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
    p1 = newPhase("../data/air-no-reactions.xml", "");
    ctml::ct2ctml("../data/air-no-reactions.cti");
    p2 = newPhase("air-no-reactions.xml", "");
    compare();
}

TEST_F(CtiConversionTest, ImplicitConversion) {
    p1 = newPhase("../data/air-no-reactions.xml", "");
    p2 = newPhase("../data/air-no-reactions.cti", "");
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
    ctml::ck2cti("pdep-test.inp");
    ThermoPhase* p = newPhase("pdep-test.cti", "");
    ASSERT_GT(p->temperature(), 0.0);
}

TEST_F(ChemkinConversionTest, MissingInputFile) {
    ASSERT_THROW(ctml::ck2cti("nonexistent-file.inp"),
                 CanteraError);
}

TEST_F(ChemkinConversionTest, FailedConversion) {
    copyInputFile("h2o2_missingThermo.inp");
    ASSERT_THROW(ctml::ck2cti("h2o2_missingThermo.inp"),
                 CanteraError);
}


#endif

} // namespace Cantera
