#include "gtest/gtest.h"
#include "cantera/thermo/CoverageDependentSurfPhase.h"
#include "cantera/thermo/ThermoFactory.h"


namespace Cantera
{

class CoverageDependentSurfPhase_Test: public testing::Test
{
public:
    CoverageDependentSurfPhase_Test() {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", "coverages_0"));
    }
    std::unique_ptr<ThermoPhase> test_phase;
    CoverageDependentSurfPhase* covdepsurf;
};

TEST_F(CoverageDependentSurfPhase_Test, construct_from_file)
{
    CoverageDependentSurfPhase* coveragedependentsurf_phase = dynamic_cast<CoverageDependentSurfPhase*>(test_phase.get());
    EXPECT_TRUE(coveragedependentsurf_phase != NULL);
}

TEST_F(CoverageDependentSurfPhase_Test, reference_enthalpies_RT)
{
    const double expected_result[10] = {
        4.895991737500073e-05,
        4.617646334628971e-05,
        5.963833377856196e-05,
        5.692092800682202e-05,
        -9.815993386192746e-05,
        6.644423752958034e-05,
        5.2492790438840376e-05,
        6.225436283744746e-05,
        -0.0035000916056617164,
        5.94838246127799e-05
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp enthalpies_RT_ref(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getEnthalpy_RT_ref(&enthalpies_RT_ref[0]);
        EXPECT_NEAR(enthalpies_RT_ref[0], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, reference_entropies_R)
{
    const double expected_result[10] = {
        17.30376125156932,
        18.782927044715542,
        16.336043707069322,
        19.692734172570507,
        13.687172677897806,
        21.46362710555691,
        16.985994108560934,
        20.063327006288937,
        9.714347806292915,
        15.544459107650459
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp entropies_R_ref(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getEntropy_R_ref(&entropies_R_ref[0]);
        EXPECT_NEAR(entropies_R_ref[2], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, reference_cp_R)
{
    const double expected_result[10] = {
        2.3911891432499996,
        2.6169441249003906,
        2.193392723419327,
        2.710068009615135,
        1.4902296533617059,
        2.80663483290828,
        2.3304807267205296,
        2.7387444952834965,
        0.37775961493433574,
        2.0046505063190665
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp cps_R_ref(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getCp_R_ref(&cps_R_ref[0]);
        EXPECT_NEAR(cps_R_ref[3], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, reference_gibbs_RT)
{
    const double expected_result[10] = {
        -35.75544525672188,
        -29.27051136018905,
        -41.389632097578534,
        -26.191806426905902,
        -66.68729427601727,
        -21.54814904209595,
        -37.46076589723343,
        -25.089092853962732,
        -172.1295537915689,
        -47.14174199946188
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp gibbs_RT_ref(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getGibbs_RT_ref(&gibbs_RT_ref[0]);
        EXPECT_NEAR(gibbs_RT_ref[4], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_enthalpies_RT)
{
    const double expected_result[10] = {
        4.895991737500073e-05,
        4.617646334628971e-05,
        5.963833377856196e-05,
        5.692092800682202e-05,
        -9.815993386192746e-05,
        6.644423752958034e-05,
        5.2492790438840376e-05,
        6.225436283744746e-05,
        -0.0035000916056617164,
        5.94838246127799e-05
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp enthalpies_RT(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getEnthalpy_RT(&enthalpies_RT[0]);
        EXPECT_NEAR(enthalpies_RT[0], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_entropies_R)
{
    const double expected_result[10] = {
        7.262785047892177,
        -211.88704374580803,
        -104.49781910237925,
        -264.41047258705385,
        -564.1004252619438,
        -668.1159683605337,
        -191.9493546909492,
        -212.12619995755966,
        -143.45731583540848,
        -390.8228284931619
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp entropies_R(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getEntropy_R(&entropies_R[0]);
        EXPECT_NEAR(entropies_R[1], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_cp_R)
{
    const double expected_result[10] = {
        2.3911891432499996,
        3.1575801771450513,
        2.2353566973964605,
        2.7752564954964587,
        2.452060193173945,
        2.8440413820654813,
        2.42326223769133,
        2.9858976609594317,
        0.5113627141915665,
        2.218066003239782
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp cps_R(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getCp_R(&cps_R[0]);
        EXPECT_NEAR(cps_R[3], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_gibbs_RT)
{
    const double expected_result[10] = {
        11.8756216200221,
        2202.4396179996957,
        -477.3453796516724,
        -316.51179756543485,
        3907.895845540286,
        -1324.569049068033,
        -13.293752588541661,
        922.4725337308008,
        1369.7551786229228,
        795.2766107464173
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp gibbs_RT(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getGibbs_RT(&gibbs_RT[0]);
        EXPECT_NEAR(gibbs_RT[3], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_gibbs)
{
    const double expected_result[10] = {
        -148643656.49121934,
        -228451988.4293846,
        -403796393.1146004,
        -634513292.1906644,
        -355109634.2098712,
        -4640855820.12355,
        -660054096.9197942,
        -744616318.8616358,
        -161367473.42965755,
        -1000153036.6544077
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp gibbs(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getPureGibbs(&gibbs[0]);
        EXPECT_NEAR(gibbs[4], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_chempotentials)
{
    const double expected_result[10] = {
        -274870452.00677145,
        856519725.4301237,
        145886101.04112163,
        1332157691.3398027,
        960092570.6206956,
        4779252095.107842,
        508906879.58705205,
        1080088611.4436927,
        -145155906.70495453,
        952738760.7658446
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp chempotentials_st(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getStandardChemPotentials(&chempotentials_st[0]);
        EXPECT_NEAR(chempotentials_st[1], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, partial_molar_enthalpies)
{
    const double expected_result[10] = {
        -244677374.6145806,
        -244559589.87957713,
        -225109949.88732946,
        -230926742.67340738,
        -259357323.04774666,
        -236935699.79515216,
        -249650880.85754883,
        -239699592.35303646,
        -264552232.7353466,
        -256393416.2056975
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp partialmolar_enthalpies(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getPartialMolarEnthalpies(&partialmolar_enthalpies[0]);
        EXPECT_NEAR(partialmolar_enthalpies[1], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, partial_molar_entropies)
{
    const double expected_result[10] = {
        66149.30110602578,
        -1744763.5646839596,
        -865388.4229791872,
        -2181467.65036173,
        -4668081.580555751,
        -5530117.33952245,
        -1584429.4415132466,
        -1746752.0200664676,
        -1178512.8623860686,
        -3211192.2824656786
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp partialmolar_entropies(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getPartialMolarEntropies(&partialmolar_entropies[0]);
        EXPECT_NEAR(partialmolar_entropies[1], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, partial_molar_cp)
{
    const double expected_result[10] = {
        52706.65227577072,
        57448.35040870828,
        49270.05260539342,
        59840.103776506076,
        40036.68657757771,
        63083.2060679071,
        51595.35078424445,
        60671.1013954003,
        30734.96527336805,
        46407.802588887964
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp partialmolar_cps(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getPartialMolarCp(&partialmolar_cps[0]);
        EXPECT_NEAR(partialmolar_cps[2], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, chemical_potentials)
{
    const double expected_result[10] = {
        -2822343946.108881,
        11440797412.804848,
        -1701950970.6310093,
        -1882300719.0459917,
        8447100310.246291,
        -9961383229.998562,
        -59098087.381852254,
        5731657082.504713,
        1139068374.937362,
        2456711091.0340676
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};
    vector_fp chempotentials(5);

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        test_phase->getChemPotentials(&chempotentials[0]);
        EXPECT_NEAR(chempotentials[3], expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, enthalpy_mole)
{
    const double expected_result[10] = {
        -122338585.53843959,
        -59309366.65989761,
        -193558162.68747127,
        -199642681.1436115,
        40537434.57875116,
        -142791039.41545668,
        -142518074.5167915,
        -108002792.86461554,
        -97942472.24175003,
        -137479276.4811465
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        EXPECT_NEAR(test_phase->enthalpy_mole(), expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, entropy_mole)
{
    const double expected_result[10] = {
        35956.532917352175,
        -10176083.957783429,
        -1043623.7221215261,
        -12891288.513604177,
        -22018125.06619569,
        4101773.9070109557,
        -2257818.910081016,
        -5537608.314706989,
        -5452194.621248438,
        -5486169.579651305
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        EXPECT_NEAR(test_phase->entropy_mole(), expected_result[i], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, cp_mole)
{
    const double expected_result[10] = {
        20304.914543820392,
        29562.94897980166,
        36080.77976856744,
        44311.85551592095,
        19165.080781805933,
        28655.920516293252,
        28534.196335297413,
        32061.519786347304,
        10104.170506513487,
        28017.603867332116
    };

    std::vector<std::string> names = {"coverages_0", "coverages_1",
                                      "coverages_2", "coverages_3",
                                      "coverages_4", "coverages_5",
                                      "coverages_6", "coverages_7",
                                      "coverages_8", "coverages_9"};

    for (size_t i = 0; i < 10; i++) {
        test_phase.reset(newPhase("copt_covdepsurf_example.yaml", names[i]));
        EXPECT_NEAR(test_phase->cp_mole(), expected_result[i], 1.e-6);
    }
}

}