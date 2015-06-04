#include "gtest/gtest.h"
#include "cantera/kinetics.h"
#include "cantera/thermo/IdealGasPhase.h"

namespace Cantera
{

#ifndef HAS_NO_PYTHON
TEST(FracCoeff, ConvertFracCoeff)
{
    IdealGasPhase thermo1("../data/frac.cti", "gas");
    std::vector<ThermoPhase*> phases1;
    phases1.push_back(&thermo1);
    GasKinetics kinetics1;
    importKinetics(thermo1.xml(), phases1, &kinetics1);

    IdealGasPhase thermo2("../data/frac.xml", "gas");
    std::vector<ThermoPhase*> phases2;
    phases2.push_back(&thermo2);
    GasKinetics kinetics2;
    importKinetics(thermo2.xml(), phases2, &kinetics2);

    ASSERT_EQ(thermo2.nSpecies(), thermo1.nSpecies());
    ASSERT_EQ(kinetics2.nReactions(), kinetics1.nReactions());

    for (size_t i = 0; i < kinetics1.nReactions(); i++) {
        for (size_t k = 0; k < thermo1.nSpecies(); k++) {
            EXPECT_EQ(kinetics1.reactantStoichCoeff(k,i),
                      kinetics2.reactantStoichCoeff(k,i));
            EXPECT_EQ(kinetics1.productStoichCoeff(k,i),
                      kinetics2.productStoichCoeff(k,i));
        }
    }
}
#endif

class FracCoeffTest : public testing::Test
{
public:
    FracCoeffTest() :
        therm("../data/frac.xml", "gas")
    {
        std::vector<ThermoPhase*> phases;
        phases.push_back(&therm);
        importKinetics(therm.xml(), phases, &kin);
        therm.setState_TPX(2000, 4*OneAtm,
                            "H2O:0.5, OH:.05, H:0.1, O2:0.15, H2:0.2");
        kH2O = therm.speciesIndex("H2O");
        kH = therm.speciesIndex("H");
        kOH = therm.speciesIndex("OH");
        kO2 = therm.speciesIndex("O2");
        kH2 = therm.speciesIndex("H2");
    }
    IdealGasPhase therm;
    GasKinetics kin;
    size_t kH2O, kH, kOH, kO2, kH2;
};

TEST_F(FracCoeffTest, StoichCoeffs)
{
    EXPECT_DOUBLE_EQ(1.0, kin.reactantStoichCoeff(kH2O, 0));
    EXPECT_DOUBLE_EQ(1.4, kin.productStoichCoeff(kH, 0));
    EXPECT_DOUBLE_EQ(0.6, kin.productStoichCoeff(kOH, 0));
    EXPECT_DOUBLE_EQ(0.2, kin.productStoichCoeff(kO2, 0));

    EXPECT_DOUBLE_EQ(0.7, kin.reactantStoichCoeff(kH2, 1));
    EXPECT_DOUBLE_EQ(0.6, kin.reactantStoichCoeff(kOH, 1));
    EXPECT_DOUBLE_EQ(0.2, kin.reactantStoichCoeff(kO2, 1));
    EXPECT_DOUBLE_EQ(1.0, kin.productStoichCoeff(kH2O, 1));
}

TEST_F(FracCoeffTest, RateConstants)
{
    vector_fp kf(kin.nReactions(), 0.0);
    vector_fp kr(kin.nReactions(), 0.0);
    kin.getFwdRateConstants(&kf[0]);
    kin.getRevRateConstants(&kr[0]);

    // sum of reaction orders is 1.0; kf has units of 1/s
    EXPECT_DOUBLE_EQ(1e13, kf[0]);

    // sum of reaction orders is 3.8.
    // kf = 1e13 (mol/cm^3)^-2.8 s^-1 = 1e13*1000^-2.8 (kmol/m^3)^-2.8 s^-1
    EXPECT_NEAR(1e13*pow(1e3, -2.8), kf[1], 1e-2);

    // Reactions are irreversible
    EXPECT_DOUBLE_EQ(0.0, kr[0]);
    EXPECT_DOUBLE_EQ(0.0, kr[1]);
}

TEST_F(FracCoeffTest, RatesOfProgress)
{
    vector_fp kf(kin.nReactions(), 0.0);
    vector_fp conc(therm.nSpecies(), 0.0);
    vector_fp ropf(kin.nReactions(), 0.0);
    therm.getConcentrations(&conc[0]);
    kin.getFwdRateConstants(&kf[0]);
    kin.getFwdRatesOfProgress(&ropf[0]);

    EXPECT_DOUBLE_EQ(conc[kH2O]*kf[0], ropf[0]);
    EXPECT_DOUBLE_EQ(pow(conc[kH2], 0.8)*conc[kO2]*pow(conc[kOH],2)*kf[1],
                     ropf[1]);
}

TEST_F(FracCoeffTest, CreationDestructionRates)
{
    vector_fp ropf(kin.nReactions(), 0.0);
    vector_fp cdot(therm.nSpecies(), 0.0);
    vector_fp ddot(therm.nSpecies(), 0.0);
    kin.getFwdRatesOfProgress(&ropf[0]);
    kin.getCreationRates(&cdot[0]);
    kin.getDestructionRates(&ddot[0]);

    EXPECT_DOUBLE_EQ(ropf[0], ddot[kH2O]);
    EXPECT_DOUBLE_EQ(1.4*ropf[0], cdot[kH]);
    EXPECT_DOUBLE_EQ(0.6*ropf[0], cdot[kOH]);
    EXPECT_DOUBLE_EQ(0.2*ropf[0], cdot[kO2]);

    EXPECT_DOUBLE_EQ(0.7*ropf[1]+ropf[2], ddot[kH2]);
    EXPECT_DOUBLE_EQ(0.6*ropf[1], ddot[kOH]);
    EXPECT_DOUBLE_EQ(0.2*ropf[1]+0.5*ropf[2], ddot[kO2]);
    EXPECT_DOUBLE_EQ(ropf[1]+ropf[2], cdot[kH2O]);

    EXPECT_DOUBLE_EQ(0.0, cdot[therm.speciesIndex("O")]);
    EXPECT_DOUBLE_EQ(0.0, ddot[therm.speciesIndex("O")]);
}

TEST_F(FracCoeffTest, EquilibriumConstants)
{
    vector_fp Kc(kin.nReactions(), 0.0);
    vector_fp mu0(therm.nSpecies(), 0.0);

    kin.getEquilibriumConstants(&Kc[0]);
    therm.getGibbs_ref(&mu0[0]); // at pRef

    double deltaG0_0 = 1.4 * mu0[kH] + 0.6 * mu0[kOH] + 0.2 * mu0[kO2] - mu0[kH2O];
    double deltaG0_1 = mu0[kH2O] - 0.7 * mu0[kH2] - 0.6 * mu0[kOH] - 0.2 * mu0[kO2];

    double pRef = therm.refPressure();
    double RT = GasConstant * therm.temperature();

    // Net stoichiometric coefficients are 1.2 and -0.5
    EXPECT_NEAR(exp(-deltaG0_0/RT) * pow(pRef/RT, 1.2), Kc[0], 1e-13 * Kc[0]);
    EXPECT_NEAR(exp(-deltaG0_1/RT) * pow(pRef/RT, -0.5), Kc[1], 1e-13 * Kc[1]);
}

TEST(InterfaceReaction, CoverageDependency) {
    IdealGasPhase gas("ptcombust.cti", "gas");
    SurfPhase surf("ptcombust.cti", "Pt_surf");
    std::vector<ThermoPhase*> phases;
    phases.push_back(&gas);
    phases.push_back(&surf);
    shared_ptr<Kinetics> kin(newKineticsMgr(surf.xml(), phases));
    ASSERT_EQ(kin->nReactions(), 25);

    double T = 500;
    surf.setState_TP(T, 101325);
    surf.setCoveragesByName("PT(S):0.7, H(S):0.3");
    vector_fp kf(kin->nReactions());
    kin->getFwdRateConstants(&kf[0]);
    EXPECT_NEAR(kf[0], 4.4579e7 * pow(T, 0.5), 1e-14*kf[0]);
    // Energies in XML file are converted from J/mol to J/kmol
    EXPECT_NEAR(kf[1], 3.7e20 * exp(-(67.4e6-6e6*0.3)/(GasConstant*T)), 1e-14*kf[1]);
}

}
