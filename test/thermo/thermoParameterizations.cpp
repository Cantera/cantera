#include "gtest/gtest.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ShomatePoly.h"
#include "cantera/thermo/PDSS_HKFT.h"
#include "cantera/base/stringUtils.h"
#include "thermo_data.h"
#include <sstream>

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
    auto sO2 = make_shared<Species>("O2", parseCompString("O:2"));
    auto sH2 = make_shared<Species>("H2", parseCompString("H:2"));
    auto sH2O = make_shared<Species>("H2O", parseCompString("H:2 O:1"));
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
    // Currently broken because MultiSpeciesThermo does not enforce reference
    // pressure consistency.
    auto sO2 = make_shared<Species>("O2", parseCompString("O:2"));
    auto sH2 = make_shared<Species>("H2", parseCompString("H:2"));
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
    auto sO2 = make_shared<Species>("O2", parseCompString("O:2"));
    auto sH2 = make_shared<Species>("H2", parseCompString("H:2"));
    auto sH2O = make_shared<Species>("H2O", parseCompString("H:2 O:1"));
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
    auto sCO = make_shared<Species>("CO", parseCompString("C:1 O:1"));
    auto sCO2 = make_shared<Species>("CO2", parseCompString("C:1 O:2"));
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

TEST(Shomate, modifyOneHf298)
{
    ShomatePoly2 S(200, 6000, 101325, co2_shomate_coeffs);

    double hf = S.reportHf298();
    EXPECT_NEAR(-393.5224e6, hf, 1e4);
    double Htest = -400e6;
    S.modifyOneHf298(npos, Htest);
    double cp, h, s;
    S.updatePropertiesTemp(298.15, &cp, &h, &s);
    EXPECT_DOUBLE_EQ(Htest, h * 298.15 * GasConstant);
    EXPECT_DOUBLE_EQ(Htest, S.reportHf298());
    S.resetHf298();
    S.updatePropertiesTemp(298.15, &cp, &h, &s);
    EXPECT_DOUBLE_EQ(hf, h * 298.15 * GasConstant);
}
