#include "gtest/gtest.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
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

TEST(SpeciesThermo, NasaPoly2FromYaml1) {
    AnyMap data = AnyMap::fromYamlString(
        "model: NASA7\n"
        "reference-pressure: 1 atm\n"
        "temperature-ranges: [200, 1000, 6000]\n"
        "data:\n"
        "- [3.944031200E+00, -1.585429000E-03, 1.665781200E-05, -2.047542600E-08,\n"
        "   7.835056400E-12, 2.896617900E+03, 6.311991700E+00]\n"
        "- [4.884754200E+00, 2.172395600E-03, -8.280690600E-07, 1.574751000E-10,\n"
        "   -1.051089500E-14, 2.316498300E+03, -1.174169500E-01]\n");
    double cp_R, h_RT, s_R;
    auto st = newSpeciesThermo(data);
    st->validate("NO2");
    st->updatePropertiesTemp(300, &cp_R, &h_RT, &s_R);
    EXPECT_DOUBLE_EQ(st->refPressure(), OneAtm);
    EXPECT_DOUBLE_EQ(cp_R, 4.47823303484);
    EXPECT_DOUBLE_EQ(h_RT, 13.735827875868003);
    EXPECT_DOUBLE_EQ(s_R, 28.913447733267262);
}

TEST(SpeciesThermo, NasaPoly2FromYaml2) {
    AnyMap data = AnyMap::fromYamlString(
        "model: NASA7\n"
        "units: {pressure: atm}\n"
        "reference-pressure: 1\n"
        "temperature-ranges: [200 K, 1000 K]\n"
        "data:\n"
        "- [3.944031200E+00, -1.585429000E-03, 1.665781200E-05, -2.047542600E-08,\n"
        "   7.835056400E-12, 2.896617900E+03, 6.311991700E+00]\n");
    double cp_R, h_RT, s_R;
    auto st = newSpeciesThermo(data);
    st->validate("NO2");
    EXPECT_DOUBLE_EQ(st->refPressure(), OneAtm);
    st->updatePropertiesTemp(300, &cp_R, &h_RT, &s_R);
    EXPECT_DOUBLE_EQ(st->maxTemp(), 1000);
    EXPECT_DOUBLE_EQ(cp_R, 4.47823303484);
    EXPECT_DOUBLE_EQ(h_RT, 13.735827875868003);
    EXPECT_DOUBLE_EQ(s_R, 28.913447733267262);
}

TEST(SpeciesThermo, Shomate2FromYaml1) {
    AnyMap data = AnyMap::fromYamlString(
        "model: Shomate\n"
        "temperature-ranges: [298, 1300, 6000]\n"
        "data:\n"
        "- [25.56759, 6.096130, 4.054656, -2.671301, 0.131021, -118.0089, 227.3665]\n"
        "- [35.15070, 1.300095, -0.205921, 0.013550, -3.282780, -127.8375, 231.7120]\n");
    double cp_R, h_RT, s_R;
    auto st = newSpeciesThermo(data);
    st->validate("CO");
    st->updatePropertiesTemp(1500, &cp_R, &h_RT, &s_R);
    EXPECT_DOUBLE_EQ(st->refPressure(), OneAtm);
    EXPECT_DOUBLE_EQ(cp_R, 4.2365020788908732);
    EXPECT_DOUBLE_EQ(h_RT, -5.747000804338211);
    EXPECT_DOUBLE_EQ(s_R, 29.878974213540165);
}

TEST(SpeciesThermo, Nasa9PolyFromYaml) {
    AnyMap data = AnyMap::fromYamlString(
        "model: NASA9\n"
        "temperature-ranges: [200.00, 1000.00, 6000.0, 20000]\n"
        "reference-pressure: 1 bar\n"
        "data:\n"
        "- [2.210371497E+04, -3.818461820E+02, 6.082738360E+00, -8.530914410E-03,\n"
        "   1.384646189E-05, -9.625793620E-09, 2.519705809E-12, 7.108460860E+02,\n"
        "   -1.076003744E+01]\n"
        "- [5.877124060E+05, -2.239249073E+03, 6.066949220E+00, -6.139685500E-04,\n"
        "   1.491806679E-07,  -1.923105485E-11, 1.061954386E-15, 1.283210415E+04,\n"
        "   -1.586640027E+01]\n"
        "- [8.310139160E+08, -6.420733540E+05, 2.020264635E+02, -3.065092046E-02,\n"
        "   2.486903333E-06, -9.705954110E-11, 1.437538881E-15, 4.938707040E+06,\n"
        "   -1.672099740E+03]");
    double cp_R, h_RT, s_R;
    auto st = newSpeciesThermo(data);
    EXPECT_DOUBLE_EQ(st->refPressure(), 1e5);
    st->updatePropertiesTemp(2000, &cp_R, &h_RT, &s_R);
    EXPECT_DOUBLE_EQ(cp_R, 4.326181187976);
    EXPECT_DOUBLE_EQ(h_RT, 3.3757914517886856);
    EXPECT_DOUBLE_EQ(s_R, 30.31743870437559);
}

TEST(SpeciesThermo, ConstCpPolyFromYaml) {
    AnyMap data = AnyMap::fromYamlString(
        "model: constant-cp # was 'const_cp'\n"
        "T0: 1000 K\n"
        "h0: 9.22 kcal/mol\n"
        "s0: -3.02 cal/mol/K\n"
        "cp0: 5.95 cal/mol/K\n");
    double cp_R, h_RT, s_R;
    auto st = newSpeciesThermo(data);
    st->updatePropertiesTemp(1100, &cp_R, &h_RT, &s_R);
    EXPECT_DOUBLE_EQ(cp_R * GasConst_cal_mol_K, 5.95);
    EXPECT_DOUBLE_EQ(h_RT * GasConst_cal_mol_K * 1100, 9.22e3 + 100 * 5.95);
    EXPECT_DOUBLE_EQ(s_R * GasConst_cal_mol_K, -3.02 + 5.95 * log(1100.0/1000.0));
}

TEST(SpeciesThermo, Mu0PolyFromYaml) {
    AnyMap data = AnyMap::fromYamlString(
        "{model: piecewise-Gibbs,"
        " h0: -890 kJ/mol,"
        " dimensionless: true,"
        " data: {298.15: -363.2104, 323.15: -300}}");
    auto st = newSpeciesThermo(data);
    double cp_R, h_RT, s_R;
    st->updatePropertiesTemp(310, &cp_R, &h_RT, &s_R);
    EXPECT_DOUBLE_EQ(cp_R, -11226.315743743922);
    EXPECT_DOUBLE_EQ(h_RT, -774.43302435932878);
    EXPECT_DOUBLE_EQ(s_R, -433.36374417010006);
}
