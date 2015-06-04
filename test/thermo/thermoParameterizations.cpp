#include "gtest/gtest.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SimpleThermo.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ShomatePoly.h"
#include "thermo_data.h"

using namespace Cantera;

class SpeciesThermoInterpTypeTest : public testing::Test
{
public:
    SpeciesThermoInterpTypeTest() {
        p.addElement("H");
        p.addElement("O");
        p.addElement("C");
    }

    IdealGasPhase p;
};

// {T0, h0, s0, cp0} (in J/kmol)
double c_o2[] = {298.15, 0.0, 2.05152e5, 2.939e4};
double c_h2[] = {298.15, 0.0, 1.3068e5, 2.885e4};
double c_h2o[] = {298.15, -2.41826e8, 1.8884e5, 3.522e4};
double c_co2[] = {298.15, -3.9351e8, 2.13785e5, 3.712e4};

TEST_F(SpeciesThermoInterpTypeTest, install_const_cp)
{
    // Compare against instantiation from CTI file
    IdealGasPhase p2("../data/simplephases.cti", "simple1");
    shared_ptr<Species> sO2(new Species("O2", parseCompString("O:2")));
    shared_ptr<Species> sH2(new Species("H2", parseCompString("H:2")));
    shared_ptr<Species> sH2O(new Species("H2O", parseCompString("H:2 O:1")));
    sO2->thermo.reset(new ConstCpPoly(200, 5000, 101325, c_o2));
    sH2->thermo.reset(new ConstCpPoly(200, 5000, 101325, c_h2));
    sH2O->thermo.reset(new ConstCpPoly(200, 5000, 101325, c_h2o));
    p.addSpecies(sO2);
    p.addSpecies(sH2);
    p.addSpecies(sH2O);
    p.initThermo();
    p2.setState_TPX(298.15, 101325, "H2:0.2, O2:0.7, H2O:0.1");
    p.setState_TPX(298.15, 101325, "H2:0.2, O2:0.7, H2O:0.1");
    EXPECT_DOUBLE_EQ(p2.meanMolecularWeight(), p.meanMolecularWeight());
    EXPECT_DOUBLE_EQ(p2.enthalpy_mass(), p.enthalpy_mass());
    EXPECT_DOUBLE_EQ(p2.entropy_mass(), p.entropy_mass());
    EXPECT_DOUBLE_EQ(p2.cp_mass(), p.cp_mass());
}

TEST_F(SpeciesThermoInterpTypeTest, DISABLED_install_bad_pref)
{
    // Currently broken because GeneralSpeciesThermo does not enforce reference
    // pressure consistency.
    shared_ptr<Species> sO2(new Species("O2", parseCompString("O:2")));
    shared_ptr<Species> sH2(new Species("H2", parseCompString("H:2")));
    sO2->thermo.reset(new ConstCpPoly(200, 5000, 101325, c_o2));
    sH2->thermo.reset(new ConstCpPoly(200, 5000, 100000, c_h2));
    p.addSpecies(sO2);
    // Pref does not match
    ASSERT_THROW(p.addSpecies(sH2), CanteraError);
}

TEST_F(SpeciesThermoInterpTypeTest, install_nasa)
{
    // Compare against instantiation from CTI file
    IdealGasPhase p2("../data/simplephases.cti", "nasa1");
    shared_ptr<Species> sO2(new Species("O2", parseCompString("O:2")));
    shared_ptr<Species> sH2(new Species("H2", parseCompString("H:2")));
    shared_ptr<Species> sH2O(new Species("H2O", parseCompString("H:2 O:1")));
    sO2->thermo.reset(new NasaPoly2(200, 3500, 101325, o2_nasa_coeffs));
    sH2->thermo.reset(new NasaPoly2(200, 3500, 101325, h2_nasa_coeffs));
    sH2O->thermo.reset(new NasaPoly2(200, 3500, 101325, h2o_nasa_coeffs));
    p.addSpecies(sO2);
    p.addSpecies(sH2);
    p.addSpecies(sH2O);
    p.initThermo();
    p2.setState_TPX(900, 101325, "H2:0.2, O2:0.7, H2O:0.1");
    p.setState_TPX(900, 101325, "H2:0.2, O2:0.7, H2O:0.1");
    EXPECT_DOUBLE_EQ(p2.meanMolecularWeight(), p.meanMolecularWeight());
    EXPECT_DOUBLE_EQ(p2.enthalpy_mass(), p.enthalpy_mass());
    EXPECT_DOUBLE_EQ(p2.entropy_mass(), p.entropy_mass());
    EXPECT_DOUBLE_EQ(p2.cp_mass(), p.cp_mass());
}

TEST_F(SpeciesThermoInterpTypeTest, install_shomate)
{
    // Compare against instantiation from CTI file
    IdealGasPhase p2("../data/simplephases.cti", "shomate1");
    shared_ptr<Species> sCO(new Species("CO", parseCompString("C:1 O:1")));
    shared_ptr<Species> sCO2(new Species("CO2", parseCompString("C:1 O:2")));
    sCO->thermo.reset(new ShomatePoly2(200, 6000, 101325, co_shomate_coeffs));
    sCO2->thermo.reset(new ShomatePoly2(200, 6000, 101325, co2_shomate_coeffs));
    p.addSpecies(sCO);
    p.addSpecies(sCO2);
    p.initThermo();
    p2.setState_TPX(900, 101325, "CO:0.2, CO2:0.8");
    p.setState_TPX(900, 101325, "CO:0.2, CO2:0.8");
    EXPECT_DOUBLE_EQ(p2.meanMolecularWeight(), p.meanMolecularWeight());
    EXPECT_DOUBLE_EQ(p2.enthalpy_mass(), p.enthalpy_mass());
    EXPECT_DOUBLE_EQ(p2.entropy_mass(), p.entropy_mass());
    EXPECT_DOUBLE_EQ(p2.cp_mass(), p.cp_mass());
}
