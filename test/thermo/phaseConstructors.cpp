#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/PDSSFactory.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/PureFluidPhase.h"
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/thermo/RedlichKwongMFTP.h"
#include "cantera/thermo/IonsFromNeutralVPSSTP.h"
#include "cantera/thermo/PDSS_IonsFromNeutral.h"
#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/MargulesVPSSTP.h"
#include "cantera/thermo/LatticePhase.h"
#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/LatticeSolidPhase.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/HMWSoln.h"
#include "cantera/thermo/PDSS_HKFT.h"

#include "cantera/thermo/NasaPoly2.h"
#include "cantera/thermo/ConstCpPoly.h"
#include "cantera/thermo/ShomatePoly.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/thermo/Mu0Poly.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include <fstream>
#include "thermo_data.h"

namespace Cantera
{

shared_ptr<Species> make_species(const std::string& name,
     const std::string& composition, const double* nasa_coeffs)
{
    auto species = make_shared<Species>(name, parseCompString(composition));
    species->thermo.reset(new NasaPoly2(200, 3500, 101325, nasa_coeffs));
    return species;
}

shared_ptr<Species> make_shomate_species(const std::string& name,
     const std::string& composition, const double* shomate_coeffs)
{
    auto species = make_shared<Species>(name, parseCompString(composition));
    species->thermo.reset(new ShomatePoly(200, 3500, 101325, shomate_coeffs));
    return species;
}

shared_ptr<Species> make_shomate2_species(const std::string& name,
     const std::string& composition, const double* shomate_coeffs)
{
    auto species = make_shared<Species>(name, parseCompString(composition));
    species->thermo.reset(new ShomatePoly2(200, 3500, 101325, shomate_coeffs));
    return species;
}

shared_ptr<Species> make_species(const std::string& name,
    const std::string& composition, double h298,
    double T1, double mu1, double T2, double mu2, double pref=101325)
{
    auto species = make_shared<Species>(name, parseCompString(composition));
    double coeffs[] = {2, h298, T1, mu1*GasConstant*T1, T2, mu2*GasConstant*T2};
    species->thermo.reset(new Mu0Poly(200, 3500, pref, coeffs));
    return species;
}

shared_ptr<Species> make_const_cp_species(const std::string& name,
    const std::string& composition, double T0, double h0, double s0, double cp)
{
    auto species = make_shared<Species>(name, parseCompString(composition));
    double coeffs[] = {T0, h0, s0, cp};
    species->thermo.reset(new ConstCpPoly(200, 3500, 101325, coeffs));
    return species;
}

//! @todo Remove after Cantera 2.5 - class FixedChemPotSSTP is deprecated
class FixedChemPotSstpConstructorTest : public testing::Test
{
};

TEST_F(FixedChemPotSstpConstructorTest, fromXML)
{
    suppress_deprecation_warnings();
    std::unique_ptr<ThermoPhase> p(newPhase("../data/LiFixed.xml"));
    ASSERT_EQ((int) p->nSpecies(), 1);
    double mu;
    p->getChemPotentials(&mu);
    ASSERT_DOUBLE_EQ(-2.3e7, mu);
    make_deprecation_warnings_fatal();
}

TEST_F(FixedChemPotSstpConstructorTest, SimpleConstructor)
{
    suppress_deprecation_warnings();
    FixedChemPotSSTP p("Li", -2.3e7);
    ASSERT_EQ((int) p.nSpecies(), 1);
    double mu;
    p.getChemPotentials(&mu);
    ASSERT_DOUBLE_EQ(-2.3e7, mu);
    make_deprecation_warnings_fatal();
}

TEST(IonsFromNeutralConstructor, fromXML)
{
    std::unique_ptr<ThermoPhase> p(newPhase("../data/mock_ion.xml",
                                            "mock_ion_phase"));
    ASSERT_EQ((int) p->nSpecies(), 2);
    p->setState_TPX(500, 2e5, "K+:0.1, Cl-:0.1");
    vector_fp mu(p->nSpecies());
    p->getChemPotentials(mu.data());

    // Values for regression testing only -- no reference values known for comparison
    EXPECT_NEAR(p->density(), 1984.2507319669949, 1e-6);
    EXPECT_NEAR(p->enthalpy_mass(), -14738312.44316336, 1e-6);
    EXPECT_NEAR(mu[0], -4.66404010e+08, 1e1);
    EXPECT_NEAR(mu[1], -2.88157316e+06, 1e-1);
}

TEST(IonsFromNeutralConstructor, fromScratch)
{
    auto neutral = make_shared<MargulesVPSSTP>();
    auto sKCl = make_shomate_species("KCl(L)", "K:1 Cl:1", kcl_shomate_coeffs);
    neutral->addSpecies(sKCl);
    std::unique_ptr<PDSS_ConstVol> ssKCl(new PDSS_ConstVol());
    ssKCl->setMolarVolume(0.03757);
    neutral->installPDSS(0, std::move(ssKCl));
    neutral->initThermo();

    IonsFromNeutralVPSSTP p;
    p.setNeutralMoleculePhase(neutral);

    auto sKp = make_shared<Species>("K+", parseCompString("K:1"), 1);
    auto sClm = make_shared<Species>("Cl-", parseCompString("Cl:1"), -1);
    sClm->input["equation-of-state"]["special-species"] = true;
    sClm->input["equation-of-state"]["model"] = "ions-from-neutral-molecule";
    p.addSpecies(sKp);
    p.addSpecies(sClm);
    std::unique_ptr<PDSS_IonsFromNeutral> ssKp(new PDSS_IonsFromNeutral());
    std::unique_ptr<PDSS_IonsFromNeutral> ssClm(new PDSS_IonsFromNeutral());
    ssKp->setNeutralSpeciesMultiplier("KCl(L)", 1.2);
    ssClm->setNeutralSpeciesMultiplier("KCl(L)", 1.5);
    ssClm->setSpecialSpecies();
    p.installPDSS(0, std::move(ssKp));
    p.installPDSS(1, std::move(ssClm));
    p.initThermo();

    ASSERT_EQ((int) p.nSpecies(), 2);
    p.setState_TPX(500, 2e5, "K+:0.1, Cl-:0.1");
    vector_fp mu(p.nSpecies());
    p.getChemPotentials(mu.data());

    // Values for regression testing only -- same as XML test
    EXPECT_NEAR(p.density(), 1984.2507319669949, 1e-6);
    EXPECT_NEAR(p.enthalpy_mass(), -14738312.44316336, 1e-6);
    EXPECT_NEAR(mu[0], -4.66404010e+08, 1e1);
    EXPECT_NEAR(mu[1], -2.88157316e+06, 1e-1);
}

class CtiConversionTest : public testing::Test
{
public:
    CtiConversionTest() {
        close_XML_File("all");
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

class ConstructFromScratch : public testing::Test
{
public:
    ConstructFromScratch()
        : sH2O(make_species("H2O", "H:2 O:1", h2o_nasa_coeffs))
        , sH2(make_species("H2", "H:2", h2_nasa_coeffs))
        , sO2(make_species("O2", "O:2", o2_nasa_coeffs))
        , sOH(make_species("OH", "H:1 O:1", oh_nasa_coeffs))
        , sCO(make_species("CO", "C:1 O:1", o2_nasa_coeffs))
        , sCO2(new Species("CO2", parseCompString("C:1 O:2")))
    {
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

TEST_F(ConstructFromScratch, throwUndefindElements)
{
    IdealGasPhase p;
    p.throwUndefinedElements();
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
    p.addUndefinedElements(); // default behavior

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
    EXPECT_NEAR(p.density(), 892.405, 2e-3);
    EXPECT_NEAR(p.enthalpy_mole(), -404848641.9409, 1e-3);

    p.setMoleFractionsByName("CO2:.6, H2O:0.02, H2:0.38");
    p.setState_TP(350, 180*OneAtm);
    EXPECT_NEAR(p.density(), 181.565, 2e-3);
    EXPECT_NEAR(p.gibbs_mass(), -1.0607e7, 2e3);
}

TEST_F(ConstructFromScratch, RedlichKwongMFTP_missing_coeffs)
{
    RedlichKwongMFTP p;
    p.addSpecies(sH2O);
    p.addSpecies(sCO2);
    p.addSpecies(sH2);
    double fa = toSI("bar-cm6/mol2");
    double fb = toSI("cm3/mol");
    p.setSpeciesCoeffs("H2O", 1.7458e8 * fa, -8e4 * fa, 18.18 * fb);
    p.setSpeciesCoeffs("H2", 30e7 * fa, -330e4 * fa, 31 * fb);
    EXPECT_THROW(p.setState_TP(300, 200e5), CanteraError);
}

//! @todo Remove after Cantera 2.5 - "gas" mode of IdealSolnGasVPSS is
//!     deprecated
TEST_F(ConstructFromScratch, IdealSolnGasVPSS_gas)
{
    suppress_deprecation_warnings();
    IdealSolnGasVPSS p;
    p.addSpecies(sH2O);
    p.addSpecies(sH2);
    p.addSpecies(sO2);
    std::unique_ptr<PDSS> pH2O(newPDSS("ideal-gas"));
    std::unique_ptr<PDSS> pH2(newPDSS("ideal-gas"));
    std::unique_ptr<PDSS> pO2(newPDSS("ideal-gas"));
    p.installPDSS(0, std::move(pH2O));
    p.installPDSS(1, std::move(pH2));
    p.installPDSS(2, std::move(pO2));

    p.setGasMode();
    EXPECT_THROW(p.setStandardConcentrationModel("unity"), CanteraError);
    p.initThermo();

    p.setState_TPX(400, 5*OneAtm, "H2:0.01, O2:0.99");
    p.equilibrate("HP");

    EXPECT_NEAR(p.temperature(), 479.929, 1e-3); // based on h2o2.cti
    EXPECT_NEAR(p.moleFraction("H2O"), 0.01, 1e-4);
    EXPECT_NEAR(p.moleFraction("H2"), 0.0, 1e-4);
    make_deprecation_warnings_fatal();
}

TEST(PureFluidFromScratch, CarbonDioxide)
{
    PureFluidPhase p;
    auto sCO2 = make_shared<Species>("CO2", parseCompString("C:1 O:2"));
    sCO2->thermo.reset(new ShomatePoly2(200, 6000, 101325, co2_shomate_coeffs));
    p.addSpecies(sCO2);
    p.setSubstance("carbon-dioxide");
    p.initThermo();
    p.setState_Tsat(280, 0.5);
    EXPECT_NEAR(p.pressure(), 4160236.987, 1e-2);
}

TEST(WaterSSTP, fromScratch)
{
    WaterSSTP water;
    water.addSpecies(make_species("H2O", "H:2, O:1", h2o_nasa_coeffs));
    water.initThermo();
    water.setState_TP(298.15, 1e5);
    EXPECT_NEAR(water.enthalpy_mole() / 1e6, -285.83, 2e-2);
}

TEST(IdealMolalSoln, fromScratch)
{
    IdealMolalSoln p;
    p.addSpecies(make_species("H2O(l)", "H:2, O:1", h2_nasa_coeffs));
    p.addSpecies(make_species("CO2(aq)", "C:1, O:2", h2_nasa_coeffs));
    p.addSpecies(make_species("H2S(aq)", "H:2, S:1", h2_nasa_coeffs));
    p.addSpecies(make_species("CH4(aq)", "C:1, H:4", h2_nasa_coeffs));
    size_t k = 0;
    for (double v : {1.5, 1.3, 0.1, 0.1}) {
        std::unique_ptr<PDSS_ConstVol> ss(new PDSS_ConstVol());
        ss->setMolarVolume(v);
        p.installPDSS(k++, std::move(ss));
    }
    p.setStandardConcentrationModel("solvent_volume");
    p.setCutoffModel("polyexp");
    // These propreties probably shouldn't be public
    p.IMS_X_o_cutoff_ = 0.20;
    p.IMS_gamma_o_min_ = 0.00001;
    p.IMS_gamma_k_min_ = 10.0;
    p.IMS_slopefCut_ = 0.6;
    p.IMS_slopegCut_ = 0.0;
    p.IMS_cCut_ = .05;
    p.initThermo();
    p.setState_TPM(298.15, OneAtm, "CH4(aq):0.01, H2S(aq):0.03, CO2(aq):0.1");

    EXPECT_NEAR(p.enthalpy_mole(), 0.013282, 1e-6);
    EXPECT_NEAR(p.gibbs_mole(), -3.8986e7, 1e3);
    EXPECT_NEAR(p.density(), 12.058, 1e-3);
}

TEST(DebyeHuckel, fromScratch)
{
    DebyeHuckel p;
    auto sH2O = make_species("H2O(l)", "H:2, O:1", h2oliq_nasa_coeffs);
    auto sNa = make_species("Na+", "Na:1, E:-1", -240.34e6,
                              298.15, -103.98186, 333.15, -103.98186);
    sNa->charge = 1;
    sNa->input["Debye-Huckel"]["ionic-radius"] = 4.0e-10;
    auto sCl = make_species("Cl-", "Cl:1, E:1", -167.08e6,
                              298.15, -74.20664, 333.15, -74.20664);
    sCl->charge = -1;
    sCl->input["Debye-Huckel"]["ionic-radius"] = 3.0e-10;
    auto sH = make_species("H+", "H:1, E:-1", 0.0, 298.15, 0.0, 333.15, 0.0);
    sH->charge = 1;
    sH->input["Debye-Huckel"]["ionic-radius"] = 9.0e-10;
    auto sOH = make_species("OH-", "O:1, H:1, E:1", -230.015e6,
                              298.15, -91.50963, 333.15, -85);
    sOH->charge = -1;
    sOH->input["Debye-Huckel"]["ionic-radius"] = 3.5e-10;
    auto sNaCl = make_species("NaCl(aq)", "Na:1, Cl:1", -96.03e6*4.184,
                              298.15, -174.5057463, 333.15, -174.5057463);
    sNaCl->input["Debye-Huckel"]["weak-acid-charge"] = -1.0;
    sNaCl->input["Debye-Huckel"]["electrolyte-species-type"] = "weakAcidAssociated";
    for (auto& s : {sH2O, sNa, sCl, sH, sOH, sNaCl}) {
        p.addSpecies(s);
    }
    std::unique_ptr<PDSS_Water> ss(new PDSS_Water());
    p.installPDSS(0, std::move(ss));
    size_t k = 1;
    for (double v : {1.3, 1.3, 0.0, 1.3, 1.3}) {
        std::unique_ptr<PDSS_ConstVol> ss(new PDSS_ConstVol());
        ss->setMolarVolume(v);
        p.installPDSS(k++, std::move(ss));
    }
    p.setDebyeHuckelModel("bdot_with_variable_a");
    p.setA_Debye(1.172576);
    p.setB_Debye(3.2864e9);
    p.setDefaultIonicRadius(3.5e-10);
    p.setMaxIonicStrength(3.0);
    p.useHelgesonFixedForm();
    p.initThermo();
    p.setState_TPM(300, 101325, "Na+:9.3549, Cl-:9.3549, H+:1.0499E-8,"
        "OH-:1.3765E-6,NaCl(aq):0.98492");

    // Regression test based on XML input file
    EXPECT_NEAR(p.density(), 60.296, 1e-2);
    EXPECT_NEAR(p.cp_mass(), 1.58216e5, 2e0);
    EXPECT_NEAR(p.entropy_mass(), 4.01279e3, 2e-2);
    vector_fp actcoeff(p.nSpecies());
    vector_fp mu_ss(p.nSpecies());
    p.getMolalityActivityCoefficients(actcoeff.data());
    p.getStandardChemPotentials(mu_ss.data());
    double act_ref[] = {1.21762, 0.538061, 0.472329, 0.717707, 0.507258, 1.0};
    double mu_ss_ref[] = {-3.06816e+08, -2.57956e+08, -1.84117e+08, 0.0,
        -2.26855e+08, -4.3292e+08};
    for (size_t k = 0; k < p.nSpecies(); k++) {
        EXPECT_NEAR(actcoeff[k], act_ref[k], 1e-5);
        EXPECT_NEAR(mu_ss[k], mu_ss_ref[k], 1e3);
    }
}

TEST(MargulesVPSSTP, fromScratch)
{
    MargulesVPSSTP p;
    auto sKCl = make_shomate_species("KCl(L)", "K:1 Cl:1", kcl_shomate_coeffs);
    auto sLiCl = make_shomate_species("LiCl(L)", "Li:1 Cl:1", licl_shomate_coeffs);
    p.addSpecies(sKCl);
    p.addSpecies(sLiCl);
    size_t k = 0;
    for (double v : {0.03757, 0.020304}) {
        std::unique_ptr<PDSS_ConstVol> ss(new PDSS_ConstVol());
        ss->setMolarVolume(v);
        p.installPDSS(k++, std::move(ss));
    }
    p.initThermo();
    p.setState_TPX(900, 101325, "KCl(L):0.3, LiCl(L):0.7");
    p.addBinaryInteraction("KCl(L)", "LiCl(L)",
        -1.757e7, -3.77e5, -7.627e3, 4.958e3, 0.0, 0.0, 0.0, 0.0);

    // Regression test based on LiKCl_liquid.xml
    EXPECT_NEAR(p.density(), 2041.9831422315356, 1e-9);
    EXPECT_NEAR(p.gibbs_mass(), -9683614.0890585743, 1e-5);
    EXPECT_NEAR(p.cp_mole(), 67478.48085733457, 1e-9);
}

TEST(LatticeSolidPhase, fromScratch)
{
    auto base = make_shared<StoichSubstance>();
    base->setName("Li7Si3(S)");
    double rho = 1390.0;
    base->setParameters(1, &rho);
    auto sLi7Si3 = make_shomate2_species("Li7Si3(S)", "Li:7 Si:3", li7si3_shomate_coeffs);
    base->addSpecies(sLi7Si3);
    base->initThermo();

    auto interstital = make_shared<LatticePhase>();
    interstital->setName("Li7Si3_Interstitial");
    auto sLii = make_const_cp_species("Li(i)", "Li:1", 298.15, 0, 2e4, 2e4);
    auto sVac = make_const_cp_species("V(i)", "", 298.15, 8.98e4, 0, 0);
    sLii->extra["molar_volume"] = 0.2;
    interstital->setSiteDensity(10.46344);
    interstital->addSpecies(sLii);
    interstital->addSpecies(sVac);
    interstital->initThermo();
    interstital->setMoleFractionsByName("Li(i):0.01 V(i):0.99");

    LatticeSolidPhase p;
    p.addLattice(base);
    p.addLattice(interstital);
    p.setLatticeStoichiometry(parseCompString("Li7Si3(S):1.0 Li7Si3_Interstitial:1.0"));
    p.initThermo();
    p.setState_TP(725, 10 * OneAtm);

    // Regression test based on modified version of Li7Si3_ls.xml
    EXPECT_NEAR(p.enthalpy_mass(), -2077955.0584538165, 1e-6);
    double mu_ref[] = {-4.62717474e+08, -4.64248485e+07, 1.16370186e+05};
    double vol_ref[] = {0.095564748201438857, 0.2, 0.095570863884152812};
    vector_fp mu(p.nSpecies());
    vector_fp vol(p.nSpecies());
    p.getChemPotentials(mu.data());
    p.getPartialMolarVolumes(vol.data());

    for (size_t k = 0; k < p.nSpecies(); k++) {
        EXPECT_NEAR(mu[k], mu_ref[k], 1e-7*fabs(mu_ref[k]));
        EXPECT_NEAR(vol[k], vol_ref[k], 1e-7);
    }
}

TEST(IdealSolidSolnPhase, fromScratch)
{
    // Regression test based fictitious XML input file
    IdealSolidSolnPhase p;
    auto sp1 = make_species("sp1", "C:2, H:2", o2_nasa_coeffs);
    sp1->extra["molar_volume"] = 1.5;
    auto sp2 = make_species("sp2", "C:1", h2o_nasa_coeffs);
    sp2->extra["molar_volume"] = 1.3;
    auto sp3 = make_species("sp3", "H:2", h2_nasa_coeffs);
    sp3->extra["molar_volume"] = 0.1;
    for (auto& s : {sp1, sp2, sp3}) {
        p.addSpecies(s);
    }
    p.setState_TPX(500, 2e5, "sp1:0.1, sp2:0.89, sp3:0.01");
    EXPECT_NEAR(p.density(), 10.1787080, 1e-6);
    EXPECT_NEAR(p.enthalpy_mass(), -15642788.8547624, 1e-4);
    EXPECT_NEAR(p.gibbs_mole(), -313642312.7114608, 1e-4);
}

static void set_hmw_interactions(HMWSoln& p) {
    double beta0_nacl[] = {0.0765, 0.008946, -3.3158E-6, -777.03, -4.4706};
    double beta1_nacl[] = {0.2664, 6.1608E-5, 1.0715E-6, 0.0, 0.0};
    double beta2_nacl[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double cphi_nacl[] = {0.00127, -4.655E-5, 0.0, 33.317, 0.09421};
    p.setBinarySalt("Na+", "Cl-", 5, beta0_nacl, beta1_nacl, beta2_nacl,
        cphi_nacl, 2.0, 0.0);

    double beta0_hcl[] = {0.1775, 0.0, 0.0, 0.0, 0.0};
    double beta1_hcl[] = {0.2945, 0.0, 0.0, 0.0, 0.0};
    double beta2_hcl[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double cphi_hcl[] = {0.0008, 0.0, 0.0, 0.0, 0.0};
    p.setBinarySalt("H+", "Cl-", 5, beta0_hcl, beta1_hcl, beta2_hcl,
        cphi_hcl, 2.0, 0.0);

    double beta0_naoh[] = {0.0864, 0.0, 0.0, 0.0, 0.0};
    double beta1_naoh[] = {0.253, 0.0, 0.0, 0.0, 0.0};
    double beta2_naoh[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    double cphi_naoh[] = {0.0044, 0.0, 0.0, 0.0, 0.0};
    p.setBinarySalt("Na+", "OH-", 5, beta0_naoh, beta1_naoh, beta2_naoh,
        cphi_naoh, 2.0, 0.0);

    double theta_cloh[] = {-0.05, 0.0, 0.0, 0.0, 0.0};
    double psi_nacloh[] = {-0.006, 0.0, 0.0, 0.0, 0.0};
    double theta_nah[] = {0.036, 0.0, 0.0, 0.0, 0.0};
    double psi_clnah[] = {-0.004, 0.0, 0.0, 0.0, 0.0};
    p.setTheta("Cl-", "OH-", 5, theta_cloh);
    p.setPsi("Na+", "Cl-", "OH-", 5, psi_nacloh);
    p.setTheta("Na+", "H+", 5, theta_nah);
    p.setPsi("Cl-", "Na+", "H+", 5, psi_clnah);
}

TEST(HMWSoln, fromScratch)
{
    // Regression test based on HMW_test_3
    HMWSoln p;
    auto sH2O = make_species("H2O(l)", "H:2, O:1", h2oliq_nasa_coeffs);
    auto sCl = make_species("Cl-", "Cl:1, E:1", 0.0,
                            298.15, -52.8716, 333.15, -52.8716, 1e5);
    sCl->charge = -1;
    auto sH = make_species("H+", "H:1, E:-1", 0.0, 298.15, 0.0, 333.15, 0.0, 1e5);
    sH->charge = 1;
    auto sNa = make_species("Na+", "Na:1, E:-1", 0.0,
                            298.15, -125.5213, 333.15, -125.5213, 1e5);
    sNa->charge = 1;
    auto sOH = make_species("OH-", "O:1, H:1, E:1", 0.0,
                            298.15, -91.523, 333.15, -91.523, 1e5);
    sOH->charge = -1;
    for (auto& s : {sH2O, sCl, sH, sNa, sOH}) {
        p.addSpecies(s);
    }
    std::unique_ptr<PDSS_Water> ss(new PDSS_Water());
    p.installPDSS(0, std::move(ss));
    size_t k = 1;
    for (double v : {1.3, 1.3, 1.3, 1.3}) {
        std::unique_ptr<PDSS_ConstVol> ss(new PDSS_ConstVol());
        ss->setMolarVolume(v);
        p.installPDSS(k++, std::move(ss));
    }
    p.setPitzerTempModel("complex");
    p.setA_Debye(1.175930);
    p.initThermo();

    set_hmw_interactions(p);
    p.setMolalitiesByName("Na+:6.0997 Cl-:6.0996986044628 H+:2.1628E-9 OH-:1.3977E-6");
    p.setState_TP(150 + 273.15, 101325);

    size_t N = p.nSpecies();
    vector_fp acMol(N), mf(N), activities(N), moll(N), mu0(N);
    p.getMolalityActivityCoefficients(acMol.data());
    p.getMoleFractions(mf.data());
    p.getActivities(activities.data());
    p.getMolalities(moll.data());
    p.getStandardChemPotentials(mu0.data());

    double acMolRef[] = {0.9341, 1.0191, 3.9637, 1.0191, 0.4660};
    double mfRef[] = {0.8198, 0.0901, 0.0000, 0.0901, 0.0000};
    double activitiesRef[] = {0.7658, 6.2164, 0.0000, 6.2164, 0.0000};
    double mollRef[] = {55.5093, 6.0997, 0.0000, 6.0997, 0.0000};
    double mu0Ref[] = {-317.175791, -186.014570, 0.0017225, -441.615456, -322.000432}; // kJ/gmol

    for (size_t k = 0 ; k < N; k++) {
        EXPECT_NEAR(acMol[k], acMolRef[k], 2e-4);
        EXPECT_NEAR(mf[k], mfRef[k], 2e-4);
        EXPECT_NEAR(activities[k], activitiesRef[k], 2e-4);
        EXPECT_NEAR(moll[k], mollRef[k], 2e-4);
        EXPECT_NEAR(mu0[k]/1e6, mu0Ref[k], 2e-6);
    }
}

TEST(HMWSoln, fromScratch_HKFT)
{
    HMWSoln p;
    auto sH2O = make_species("H2O(l)", "H:2, O:1", h2oliq_nasa_coeffs);
    auto sNa = make_species("Na+", "Na:1, E:-1", 0.0,
                            298.15, -125.5213, 333.15, -125.5213, 1e5);
    sNa->charge = 1;
    auto sCl = make_species("Cl-", "Cl:1, E:1", 0.0,
                            298.15, -52.8716, 333.15, -52.8716, 1e5);
    sCl->charge = -1;
    auto sH = make_species("H+", "H:1, E:-1", 0.0, 298.15, 0.0, 333.15, 0.0, 1e5);
    sH->charge = 1;
    auto sOH = make_species("OH-", "O:1, H:1, E:1", 0.0,
                            298.15, -91.523, 333.15, -91.523, 1e5);
    sOH->charge = -1;
    for (auto& s : {sH2O, sNa, sCl, sH, sOH}) {
        p.addSpecies(s);
    }
    double h0[] = {-57433, Undef, 0.0, -54977};
    double g0[] = {Undef, -31379, 0.0, -37595};
    double s0[] = {13.96, 13.56, Undef, -2.56};
    double a[][4] = {{0.1839, -228.5, 3.256, -27260},
                     {0.4032, 480.1, 5.563, -28470},
                     {0.0, 0.0, 0.0, 0.0},
                     {0.12527, 7.38, 1.8423, -27821}};
    double c[][2] = {{18.18, -29810}, {-4.4, -57140}, {0.0, 0.0}, {4.15, -103460}};
    double omega[] = {33060, 145600, 0.0, 172460};

    std::unique_ptr<PDSS_Water> ss(new PDSS_Water());
    p.installPDSS(0, std::move(ss));
    for (size_t k = 0; k < 4; k++) {
        std::unique_ptr<PDSS_HKFT> ss(new PDSS_HKFT());
        if (h0[k] != Undef) {
            ss->setDeltaH0(h0[k] * toSI("cal/gmol"));
        }
        if (g0[k] != Undef) {
            ss->setDeltaG0(g0[k] * toSI("cal/gmol"));
        }
        if (s0[k] != Undef) {
            ss->setS0(s0[k] * toSI("cal/gmol/K"));
        }
        a[k][0] *= toSI("cal/gmol/bar");
        a[k][1] *= toSI("cal/gmol");
        a[k][2] *= toSI("cal-K/gmol/bar");
        a[k][3] *= toSI("cal-K/gmol");
        c[k][0] *= toSI("cal/gmol/K");
        c[k][1] *= toSI("cal-K/gmol");
        ss->set_a(a[k]);
        ss->set_c(c[k]);
        ss->setOmega(omega[k] * toSI("cal/gmol"));
        p.installPDSS(k+1, std::move(ss));
    }
    p.setPitzerTempModel("complex");
    p.setA_Debye(-1);
    p.initThermo();

    set_hmw_interactions(p);
    p.setMolalitiesByName("Na+:6.0954 Cl-:6.0954 H+:2.1628E-9 OH-:1.3977E-6");
    p.setState_TP(50 + 273.15, 101325);

    size_t N = p.nSpecies();
    vector_fp mv(N), h(N), mu(N), ac(N), acoeff(N);
    p.getPartialMolarVolumes(mv.data());
    p.getPartialMolarEnthalpies(h.data());
    p.getChemPotentials(mu.data());
    p.getActivities(ac.data());
    p.getActivityCoefficients(acoeff.data());

    double mvRef[] = {0.01815224, 0.00157182, 0.01954605, 0.00173137, -0.0020266};

    for (size_t k = 0; k < N; k++) {
        EXPECT_NEAR(mv[k], mvRef[k], 2e-8);
    }
}

TEST(PDSS_SSVol, fromScratch)
{
    // Regression test based on comparison with using XML input file
    IdealSolnGasVPSS p;
    double coeffs[] = {700.0, 26.3072, 30.4657, -69.1692, 44.1951, 0.0776,
        -6.0337, 59.8106, 22.6832, 10.476, -6.5428, 1.3255, 0.8783, -2.0426,
        62.8859};
    auto sLi = make_shomate2_species("Li(L)", "Li:1", coeffs);
    p.addSpecies(sLi);
    p.setSolnMode();
    p.setStandardConcentrationModel("unity");
    std::unique_ptr<PDSS_SSVol> ss(new PDSS_SSVol());
    double rho_coeffs[] = {536.504, -1.04279e-1, 3.84825e-6, -5.2853e-9};
    ss->setDensityPolynomial(rho_coeffs);
    p.installPDSS(0, std::move(ss));
    p.initThermo();
    p.setState_TP(300, OneAtm);
    EXPECT_NEAR(p.density(), 505.42393940, 2e-8);
    EXPECT_NEAR(p.gibbs_mole(), -7801634.1184443515, 2e-8);
    p.setState_TP(400, 2*OneAtm);
    EXPECT_NEAR(p.density(), 495.06986080, 2e-8);
    EXPECT_NEAR(p.molarVolume(), 0.014018223587243668, 2e-12);
    p.setState_TP(500, 2*OneAtm);
    EXPECT_NEAR(p.density(), 484.66590, 2e-8);
    EXPECT_NEAR(p.enthalpy_mass(), 1236701.0904197122, 2e-8);
    EXPECT_NEAR(p.entropy_mole(), 49848.488477407751, 2e-8);
}

TEST(Species, fromYaml)
{
    AnyMap spec = AnyMap::fromYamlString(
        "name: NO2\n"
        "composition: {N: 1, O: 2}\n"
        "units: {length: cm, quantity: mol}\n"
        "molar-volume: 0.536\n"
        "thermo:\n"
        "  model: NASA7\n"
        "  temperature-ranges: [200, 1000, 6000]\n"
        "  data:\n"
        "  - [3.944031200E+00, -1.585429000E-03, 1.665781200E-05, -2.047542600E-08,\n"
        "     7.835056400E-12, 2.896617900E+03, 6.311991700E+00]\n"
        "  - [4.884754200E+00, 2.172395600E-03, -8.280690600E-07, 1.574751000E-10,\n"
        "     -1.051089500E-14, 2.316498300E+03, -1.174169500E-01]\n");

    auto S = newSpecies(spec);
    EXPECT_DOUBLE_EQ(S->thermo->minTemp(), 200);
    EXPECT_EQ(S->composition.at("N"), 1);
    // Check that units directive gets propagated to `input`
    EXPECT_DOUBLE_EQ(S->input.convert("molar-volume", "m^3/kmol"), 0.000536);
}

} // namespace Cantera
