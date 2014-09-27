#include "gtest/gtest.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/base/ctml.h"
#include <fstream>

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
    p1 = newPhase("../data/air-no-reactions.xml");
    ctml::ct2ctml("../data/air-no-reactions.cti");
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
    ctml::ck2cti("pdep-test.inp");
    ThermoPhase* p = newPhase("pdep-test.cti");
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

// 2-region NASA coefficients; Order is significantly different from the
// standard NASA format.
double h2o_coeffs[] = {
    1000.0, -3.029372670E+04, -8.490322080E-01, 4.198640560E+00,
    -2.036434100E-03, 6.520402110E-06, -5.487970620E-09, 1.771978170E-12,
    -3.000429710E+04, 4.966770100E+00,  3.033992490E+00, 2.176918040E-03,
    -1.640725180E-07, -9.704198700E-11, 1.682009920E-14};
double h2o_comp[] = {2.0, 1.0, 0.0};

double h2_coeffs[] = {
    1000.0, -9.17935173E+02, 6.83010238E-01, 2.34433112E+00,
    7.98052075E-03, -1.94781510E-05, 2.01572094E-08, -7.37611761E-12,
    -9.50158922E+02, -3.20502331E+00, 3.33727920E+00, -4.94024731E-05,
    4.99456778E-07, -1.79566394E-10, 2.00255376E-14};
double h2_comp[] = {2.0, 0.0, 0.0};

double o2_coeffs[] = {
    1000.0, -1.063943560E+03, 3.657675730E+00, 3.782456360E+00,
    -2.996734160E-03, 9.847302010E-06, -9.681295090E-09, 3.243728370E-12,
    -1.088457720E+03, 5.453231290E+00, 3.282537840E+00, 1.483087540E-03,
    -7.579666690E-07, 2.094705550E-10, -2.167177940E-14};
double o2_comp[] = {0.0, 2.0, 0.0};

double oh_coeffs[] = {
    1000.0, 3.615080560E+03, -1.039254580E-01, 3.992015430E+00,
    -2.401317520E-03, 4.617938410E-06, -3.881133330E-09, 1.364114700E-12,
    3.858657000E+03, 4.476696100E+00, 3.092887670E+00, 5.484297160E-04,
    1.265052280E-07, -8.794615560E-11, 1.174123760E-14};
double oh_comp[] = {1.0, 1.0, 0.0};

// 2-region Shomate coefficients
double co2_coeffs[] = {
    1200.0, 24.99735, 55.18696, -33.69137, 7.948387, -0.136638, -403.6075, 228.2431,
    58.16639, 2.720074, -0.492289, 0.038844, -6.447293, -425.9186, 263.6125};
double co2_comp[] = {0.0, 2.0, 1.0};

class ConstructFromScratch : public testing::Test
{
public:
    ConstructFromScratch() {
        p.addElement("H");
        p.addElement("O");
        p.addElement("C");
    }
    IdealGasPhase p;
};

TEST_F(ConstructFromScratch, AddElements)
{
    ASSERT_EQ(3, p.nElements());
    ASSERT_EQ("H", p.elementName(0));
    ASSERT_EQ(1, p.elementIndex("O"));
}

TEST_F(ConstructFromScratch, AddSpeciesNasa)
{
    p.setSpeciesThermo(newSpeciesThermoMgr(NASA));
    SpeciesThermo& sp = p.speciesThermo();

    p.addUniqueSpecies("H2O", h2o_comp);
    sp.install("H2O", 0, NASA, h2o_coeffs, 200.0, 3500.0, 101325.0);
    p.addUniqueSpecies("H2", h2_comp);
    sp.install("H2", 1, NASA, h2_coeffs, 200.0, 3500.0, 101325.0);

    ASSERT_EQ(2, p.nSpecies());

    p.addUniqueSpecies("O2", o2_comp);
    sp.install("O2", 2, NASA, o2_coeffs, 200.0, 3500.0, 101325.0);
    p.addUniqueSpecies("OH", oh_comp);
    sp.install("OH", 3, NASA, oh_coeffs, 200.0, 3500.0, 101325.0);

    ASSERT_EQ(4, p.nSpecies());
    ASSERT_EQ("H2", p.speciesName(1));
    ASSERT_EQ(2, p.nAtoms(2, 1)); // O in O2
    ASSERT_EQ(2, p.nAtoms(0, 0)); // H in H2O
}

} // namespace Cantera
