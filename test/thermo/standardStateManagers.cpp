#include "gtest/gtest.h"
#include "cantera/thermo/HMWSoln.h"

using namespace Cantera;

TEST(HMW, VPSSMgrGeneral_vs_VPSSMgrWater_ConstVol)
{
    // Calculations should give the same result using either the generic
    // VPSSMgr_General class or one of the more specialized classes such as
    // VPSSMgr_Water_ConstVol.
    HMWSoln p1("../data/HMW_NaCl.xml", "water_constvol");
    HMWSoln p2("../data/HMW_NaCl.xml", "general");
    size_t n = p1.nSpecies();
    vector_fp molalities(n);
    p1.getMolalities(molalities.data());
    molalities[2] = 2.1628E-9;
    molalities[3] = 6.0997;
    molalities[4] = 1.3977E-6;
    molalities[1] = molalities[2] + molalities[3] - molalities[4];
    p1.setMolalities(molalities.data());
    p2.setMolalities(molalities.data());
    p1.setState_TP(310.15, 201325);
    p2.setState_TP(310.15, 201325);

    vector_fp v1(n);
    vector_fp v2(n);
    p1.getStandardVolumes(v1.data());
    p2.getStandardVolumes(v2.data());
    for (size_t i = 0; i < n; i++) {
        EXPECT_NEAR(v1[i], v2[i], 1e-9) << p1.speciesName(i);
    }

    p1.getCp_R(v1.data());
    p2.getCp_R(v2.data());
    for (size_t i = 0; i < n; i++) {
        EXPECT_NEAR(v1[i], v2[i], 1e-10) << p1.speciesName(i);
    }

    p1.getEntropy_R(v1.data());
    p2.getEntropy_R(v2.data());
    for (size_t i = 0; i < n; i++) {
        EXPECT_NEAR(v1[i], v2[i], 1e-10) << p1.speciesName(i);
    }

    p1.getEnthalpy_RT(v1.data());
    p2.getEnthalpy_RT(v2.data());
    for (size_t i = 0; i < n; i++) {
        EXPECT_NEAR(v1[i], v2[i], 1e-10) << p1.speciesName(i);
    }

    p1.getChemPotentials_RT(v1.data());
    p2.getChemPotentials_RT(v2.data());
    for (size_t i = 0; i < n; i++) {
        EXPECT_NEAR(v1[i], v2[i], 1e-10) << p1.speciesName(i);
    }

    EXPECT_NEAR(p1.entropy_mole(), p2.entropy_mole(), 1e-7);
    EXPECT_NEAR(p1.enthalpy_mole(), p2.enthalpy_mole(), 1e-4);
}
