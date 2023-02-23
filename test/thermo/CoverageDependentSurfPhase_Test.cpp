#include "gtest/gtest.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/base/utilities.h"
#include "cantera/thermo/CoverageDependentSurfPhase.h"
#include "cantera/thermo/ThermoFactory.h"


namespace Cantera
{

class CoverageDependentSurfPhase_Test: public testing::Test
{
public:
    CoverageDependentSurfPhase_Test():
        // Array of ten different temperatures to test if
        // temperature dependence is implemented properly
        test_Ts({500.0, 625.0, 427.0, 711.0, 260.0,
                 903.0, 475.3, 748.3, 100.1, 372.1}),
        // Array of ten different surface coverages to test if
        // coverage dependence is implemented properly
        test_covs({
                  {0.00, 0.10, 0.50, 0.10, 0.30},
                  {0.12, 0.07, 0.21, 0.17, 0.43},
                  {0.00, 0.71, 0.08, 0.07, 0.14},
                  {0.60, 0.10, 0.10, 0.10, 0.10},
                  {0.20, 0.23, 0.57, 0.00, 0.00}
                  })
    {
        // CoverageDependentSurfPhase method values will be
        // compared against SurfPhase method values or explicitly calculated values
        covdepsurf_phase.reset(newPhase("copt_covdepsurf_example.yaml", "covdep"));
        idealsurf_phase.reset(newPhase("copt_covdepsurf_example.yaml", "ideal"));
        surfphase_ptr = dynamic_cast<SurfPhase*>(covdepsurf_phase.get());
    }
    std::unique_ptr<ThermoPhase> covdepsurf_phase;
    std::unique_ptr<ThermoPhase> idealsurf_phase;
    // A pointer to call setCoverage method
    SurfPhase* surfphase_ptr;
    // To call unit convert methods
    UnitSystem us;

    vector_fp test_Ts;
    std::vector<vector_fp> test_covs;
};

TEST_F(CoverageDependentSurfPhase_Test, construct_from_yaml)
{
    // This test checks if CoverageDependentSurfPhase object is
    // properly created without errors
    CoverageDependentSurfPhase* coveragedependentsurf_phase =
        dynamic_cast<CoverageDependentSurfPhase*>(covdepsurf_phase.get());
    EXPECT_TRUE(coveragedependentsurf_phase != NULL);
}

TEST_F(CoverageDependentSurfPhase_Test, reference_enthalpies_RT)
{
    // This test checks if CoverageDependentSurfPhase reference enthalpy is
    // properly matches with SurfPhase reference enthalpy
    vector_fp enthalpies_RT_ref(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 10; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        // Get CoverageDependentSurfPhase reference enthalpy method value
        covdepsurf_phase->getEnthalpy_RT_ref(&enthalpies_RT_ref[0]);
        idealsurf_phase->setTemperature(test_Ts[i]);
        // Get SurfPhase reference enthalpy method value
        idealsurf_phase->getEnthalpy_RT_ref(&expected_result[0]);
        EXPECT_NEAR(enthalpies_RT_ref[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, reference_entropies_R)
{
    // This test checks if CoverageDependentSurfPhase reference entropy is
    // properly matches with SurfPhase reference entropy
    vector_fp entropies_R_ref(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 10; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        // Get CoverageDependentSurfPhase reference entropy method value
        covdepsurf_phase->getEntropy_R_ref(&entropies_R_ref[0]);
        idealsurf_phase->setTemperature(test_Ts[i]);
        // Get SurfPhase reference entropy method value
        idealsurf_phase->getEntropy_R_ref(&expected_result[0]);
        EXPECT_NEAR(entropies_R_ref[2], expected_result[2], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, reference_cp_R)
{
    // This test checks if CoverageDependentSurfPhase reference heat capacity is
    // properly matches with SurfPhase reference heat capacity
    vector_fp cps_R_ref(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 10; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        // Get CoverageDependentSurfPhase reference heat capacity method value
        covdepsurf_phase->getCp_R_ref(&cps_R_ref[0]);
        idealsurf_phase->setTemperature(test_Ts[i]);
        // Get SurfPhase reference heat capacity method value
        idealsurf_phase->getCp_R_ref(&expected_result[0]);
        EXPECT_NEAR(cps_R_ref[3], expected_result[3], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, reference_gibbs_RT)
{
    // This test checks if CoverageDependentSurfPhase reference gibbs free energy is
    // properly matches with SurfPhase reference gibbs free energy
    vector_fp gibbs_RT_ref(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 10; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        // Get CoverageDependentSurfPhase reference gibbs free energy method value
        covdepsurf_phase->getGibbs_RT_ref(&gibbs_RT_ref[0]);
        idealsurf_phase->setTemperature(test_Ts[i]);
        // Get SurfPhase reference gibbs free energy method value
        idealsurf_phase->getGibbs_RT_ref(&expected_result[0]);
        EXPECT_NEAR(gibbs_RT_ref[4], expected_result[4], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_enthalpies_RT)
{
    // This test checks if CoverageDependentSurfPhase standard enthalpy is
    // properly matches with explicitly calculated standard enthalpy

    // Define coverage-dependent parameters same as given in
    // copt_covdepsurf_example.yaml for explicit external calculations
    // Parameters for piecewise linear dependency
    double h_slope = us.convertFrom(0.48, "eV/molec");
    double h_low = us.convertFrom(0.5e2, "kJ/mol");
    double h_high = us.convertFrom(1.0e2, "kJ/mol");
    // Parameters for polynomial dependency
    vector_fp h_coeffs {0.0, 0.0, -3.86e4, 0.0, 4.2e5};
    for (size_t i = 0; i < h_coeffs.size(); i++) {
        h_coeffs[i] = us.convertFrom(h_coeffs[i], "J/mol");
    }
    // Parameter for interpolative dependency
    double h_int = us.convertFrom(0.75, "kcal/mol");
    // Parameters for dependency on heat capacity
    double cp_a = us.convertFrom(0.02e-1, "kJ/mol/K");
    double cp_b = us.convertFrom(-0.156e-1, "kJ/mol/K");
    double int_cp_tnow = test_Ts[5] * (cp_a * log(test_Ts[5]) - cp_a + cp_b);
    double int_cp_298 = 298.15 * (cp_a * log(298.15) - cp_a + cp_b);

    vector_fp enthalpies_RT(5);
    vector_fp expected_result(5);

    // Get CoverageDependentSurfPhase standard enthalpy method value
    covdepsurf_phase->setTemperature(test_Ts[5]);
    surfphase_ptr->setCoverages(test_covs[0].data());
    covdepsurf_phase->getEnthalpy_RT(&enthalpies_RT[0]);
    // Externally calculating value by adding coverage-dependent terms
    // to SurfPhase standard enthalpy method value
    idealsurf_phase->setTemperature(test_Ts[5]);
    idealsurf_phase->getEnthalpy_RT(&expected_result[0]);
    double RT = idealsurf_phase->RT();
    expected_result[1] += (h_slope * test_covs[0][1]) / RT;
    expected_result[1] += (h_low * 0.4) / RT;
    expected_result[1] += (h_high * (test_covs[0][2] - 0.4)) / RT;
    expected_result[1] += ((int_cp_tnow - int_cp_298)
                            * test_covs[0][2] * test_covs[0][2]) / RT;
    expected_result[1] += poly4(test_covs[0][3], h_coeffs.data()) / RT;
    expected_result[1] += h_int / RT;
    EXPECT_NEAR(enthalpies_RT[1], expected_result[1], 1.e-6);
}

TEST_F(CoverageDependentSurfPhase_Test, standard_entropies_R)
{
    // This test checks if CoverageDependentSurfPhase standard entropy is
    // properly matches with explicitly calculated standard entropy

    // Define coverage-dependent parameters same as given in
    // copt_covdepsurf_example.yaml for explicit external calculations
    // Parameters for piecewise linear dependency
    double s_slope = us.convertFrom(-0.031, "eV/molec/K");
    double s_low = us.convertFrom(0.1e2, "kJ/mol/K");
    double s_high = us.convertFrom(-0.2e2, "kJ/mol/K");
    // Parameters for polynomial dependency
    vector_fp s_coeffs {0.0, 0.8e3, 0.0, -1.26e4, 0.0};
    for (size_t i = 0; i < s_coeffs.size(); i++) {
        s_coeffs[i] = us.convertFrom(s_coeffs[i], "J/mol/K");
    }
    // Parameters for interpolative dependency
    double s_int = us.convertFrom(-0.42, "kcal/mol/K");
    // Parameters for dependency on heat capacity
    double cp_a = us.convertFrom(0.02e-1, "kJ/mol/K");
    double cp_b = us.convertFrom(-0.156e-1, "kJ/mol/K");
    double int_cp_T_tnow = log(test_Ts[3]) * (cp_a * log(test_Ts[3]) + 2 * cp_b);
    double int_cp_T_298 = log(298.15) * (cp_a * log(298.15) + 2 * cp_b);
    // Parameter for reference coverage correction
    double ref_cov = 0.22;

    vector_fp entropies_R(5);
    vector_fp expected_result(5);

    // Get CoverageDependentSurfPhase standard entropy method value
    covdepsurf_phase->setTemperature(test_Ts[3]);
    surfphase_ptr->setCoverages(test_covs[0].data());
    covdepsurf_phase->getEntropy_R(&entropies_R[0]);
    // Externally calculate value by adding coverage-dependent terms
    // to SurfPhase standard entropy method value
    idealsurf_phase->setTemperature(test_Ts[3]);
    idealsurf_phase->getEntropy_R(&expected_result[0]);
    expected_result[1] += (s_slope * test_covs[0][1]) / GasConstant;
    expected_result[1] += (s_low * 0.4) / GasConstant;
    expected_result[1] += (s_high * (test_covs[0][2] - 0.4)) / GasConstant;
    expected_result[1] += poly4(test_covs[0][3], s_coeffs.data()) / GasConstant;
    expected_result[1] += s_int / GasConstant;
    expected_result[1] += 0.5 * (int_cp_T_tnow - int_cp_T_298)
                          * test_covs[0][2] * test_covs[0][2] / GasConstant;
    expected_result[1] -= -log(ref_cov);
    EXPECT_NEAR(entropies_R[1], expected_result[1], 1.e-6);
}

TEST_F(CoverageDependentSurfPhase_Test, standard_cp_R)
{
    // This test checks if CoverageDependentSurfPhase standard heat capacity is
    // properly matches with explicitly calculated standard heat capacity

    // Define coverage-dependent parameters same as given in
    // copt_covdepsurf_example.yaml for explicit external calculations
    // Parameters for dependency on heat capacity
    double cp_a = us.convertFrom(0.02e-1, "kJ/mol/K");
    double cp_b = us.convertFrom(-0.156e-1, "kJ/mol/K");

    vector_fp cps_R(5);
    vector_fp expected_result(5);

    // Get CoverageDependentSurfPhase standard heat capacity method value
    covdepsurf_phase->setTemperature(test_Ts[4]);
    surfphase_ptr->setCoverages(test_covs[0].data());
    covdepsurf_phase->getCp_R(&cps_R[0]);
    // Externally calculating value by adding coverage-dependent terms
    // to SurfPhase standard heat capacity method value
    idealsurf_phase->setTemperature(test_Ts[4]);
    idealsurf_phase->getCp_R(&expected_result[0]);
    expected_result[1] += (cp_a * log(test_Ts[4]) + cp_b)
                           * test_covs[0][2] * test_covs[0][2] / GasConstant;
    EXPECT_NEAR(cps_R[1], expected_result[1], 1.e-6);
}

TEST_F(CoverageDependentSurfPhase_Test, standard_gibbs_RT)
{
    // This test checks if standard gibbs free energy is
    // properly matches value calculated with enthalpy and entropy
    vector_fp gibbs_RT(5);
    vector_fp expected_result(5);
    vector_fp entropies_R(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get standard gibbs free energy method value
        covdepsurf_phase->getGibbs_RT(&gibbs_RT[0]);
        // Get standard gibbs free energy with enthalpy and entropy
        covdepsurf_phase->getEnthalpy_RT(&expected_result[0]);
        covdepsurf_phase->getEntropy_R(&entropies_R[0]);
        expected_result[1] -= entropies_R[1];
        EXPECT_NEAR(gibbs_RT[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_gibbs)
{
    // This test checks if standard gibbs free energy is
    // correctly calculated from dimensionless gibbs free energy
    vector_fp gibbs(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get standard gibbs free energy method value
        covdepsurf_phase->getPureGibbs(&gibbs[0]);
        // Get standard gibbs free energy from a dimensionless value
        covdepsurf_phase->getGibbs_RT(&expected_result[0]);
        expected_result[1] *= covdepsurf_phase->RT();
        EXPECT_NEAR(gibbs[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, standard_chempotentials)
{
    // This test checks if standard chemical potential array
    // is equivalent to gibbs free energy array
    vector_fp chempotentials_st(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get standard chemical potential method value
        covdepsurf_phase->getStandardChemPotentials(&chempotentials_st[0]);
        // Get standard gibbs free energy method value
        covdepsurf_phase->getPureGibbs(&expected_result[0]);
        EXPECT_NEAR(chempotentials_st[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, partial_molar_enthalpies)
{
    // This test checks if partial molar enthalpy is correctly
    // calculated from standard enthalpy
    vector_fp partialmolar_enthalpies(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get partial molar enthalpy method value
        covdepsurf_phase->getPartialMolarEnthalpies(&partialmolar_enthalpies[0]);
        // Get standard enthalpy method value
        covdepsurf_phase->getEnthalpy_RT(&expected_result[0]);
        expected_result[1] *= covdepsurf_phase->RT();
        EXPECT_NEAR(partialmolar_enthalpies[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, partial_molar_entropies)
{
    // This test checks if partial molar entropy is correctly
    // calculated from standard entropy
    vector_fp partialmolar_entropies(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get partial molar entropy method value
        covdepsurf_phase->getPartialMolarEntropies(&partialmolar_entropies[0]);
        // Get standard entropy method value
        covdepsurf_phase->getEntropy_R(&expected_result[0]);
        expected_result[1] -= log(test_covs[i][1]);
        expected_result[1] *= GasConstant;
        EXPECT_NEAR(partialmolar_entropies[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, partial_molar_cp)
{
    // This test checks if partial molar heat capacity is correctly
    // calculated from standard heat capacity
    vector_fp partialmolar_cps(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get partial molar heat capacity method value
        covdepsurf_phase->getPartialMolarCp(&partialmolar_cps[0]);
        // Get standard heat capacity method value
        covdepsurf_phase->getCp_R(&expected_result[0]);
        expected_result[1] *= GasConstant;
        EXPECT_NEAR(partialmolar_cps[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, chemical_potentials)
{
    // This test checks if chemical potential is correctly
    // calculated from standard chemical potential
    vector_fp chempotentials(5);
    vector_fp expected_result(5);

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get chemical potential method value
        covdepsurf_phase->getChemPotentials(&chempotentials[0]);
        // Get standard chemical potential method value
        covdepsurf_phase->getStandardChemPotentials(&expected_result[0]);
        expected_result[1] += log(test_covs[i][1]) * covdepsurf_phase->RT();
        EXPECT_NEAR(chempotentials[1], expected_result[1], 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, enthalpy_mole)
{
    // This test checks if phase enthalpy is correctly
    // calculated from partial molar enthalpy
    vector_fp partialmolar_enthalpies(5);
    double expected_result;

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get partial molar enthalpy method value
        covdepsurf_phase->getPartialMolarEnthalpies(&partialmolar_enthalpies[0]);
        expected_result = covdepsurf_phase->mean_X(partialmolar_enthalpies);
        // Compare it with phase enthalpy method value
        EXPECT_NEAR(covdepsurf_phase->enthalpy_mole(), expected_result, 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, entropy_mole)
{
    // This test checks if phase entropy is correctly
    // calculated from partial molar entropy
    vector_fp partialmolar_entropies(5);
    double expected_result;

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get partial molar entropy method value
        covdepsurf_phase->getPartialMolarEntropies(&partialmolar_entropies[0]);
        expected_result = covdepsurf_phase->mean_X(partialmolar_entropies);
        // Compare it with phase entropy method value
        EXPECT_NEAR(covdepsurf_phase->entropy_mole(), expected_result, 1.e-6);
    }
}

TEST_F(CoverageDependentSurfPhase_Test, cp_mole)
{
    // This test checks if phase heat capacity is correctly
    // calculated from partial molar heat capacity
    vector_fp partialmolar_cps(5);
    double expected_result;

    for (size_t i = 0; i < 5; i++) {
        covdepsurf_phase->setTemperature(test_Ts[i]);
        surfphase_ptr->setCoverages(test_covs[i].data());
        // Get partial molar heat capacity method value
        covdepsurf_phase->getPartialMolarCp(&partialmolar_cps[0]);
        expected_result = covdepsurf_phase->mean_X(partialmolar_cps);
        // Compare it with phase heat capacity method value
        EXPECT_NEAR(covdepsurf_phase->cp_mole(), expected_result, 1.e-6);
    }
}

};
