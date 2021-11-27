#include "gtest/gtest.h"
#include "cantera/thermo/PengRobinson.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class cubicSolver_Test : public testing::Test
{
public:
    cubicSolver_Test() {
        test_phase.reset(newPhase("../data/co2_PR_example.yaml"));
    }

    //vary the composition of a co2-h2 mixture
    void set_r(const double r) {
        vector_fp moleFracs(7);
        moleFracs[0] = r;
        moleFracs[2] = 1-r;
        test_phase->setMoleFractions(&moleFracs[0]);
    }

    std::unique_ptr<ThermoPhase> test_phase;
};

#ifdef __MINGW32__
TEST_F(cubicSolver_Test, solve_cubic_DISABLED)
{
    // the following test fails on mingw: EXPECT_NEAR(nSolnValues, -2, 1.e-6);
    // where the positive root is found instead
}
#else
TEST_F(cubicSolver_Test, solve_cubic)
{
    /* This tests validates the cubic solver by considering CO2 as an example.
    *  Values of a_coeff, b_coeff and acentric factor are hard-coded.
    *  The temperature dependent parameter in P-R EoS is calculated as
    *       \alpha = [1 + \kappa(1 - sqrt{T/T_crit}]^2
    *  kappa is a function calculated based on the acentric factor.
    *
    * Three different states are considered as follows:
    * 1. T = 300 T, P = 1 bar => Vapor (1 real root of the cubic equation)
    * 2. T = 300 K, P = 80 bar => Supercritical (1 real root of the cubic equation)
    * 3. T = Tc, P = Pc => Near critical region
    */

    // Define a Peng-Robinson phase
    PengRobinson* peng_robinson_phase = dynamic_cast<PengRobinson*>(test_phase.get());
    EXPECT_TRUE(peng_robinson_phase != NULL);

    set_r(1.0);
    double a_coeff = 3.958134E+5;
    double b_coeff = 26.6275/1000;
    double acc_factor = 0.228;
    const double Tcrit = test_phase->critTemperature();
    const double pCrit = test_phase->critPressure();
    double kappa = 0.37464 + 1.54226 * acc_factor - 0.26992 * acc_factor * acc_factor;
    double temp, pres, alpha, rho, p;
    int nSolnValues;
    double Vroot[3];

    const double expected_result[3] = {
        24.809417072270659,
        0.063638901847459045,
        0.10521518550521104
    };

    //Vapor phase -> nSolnValues = 1
    temp = 300;
    pres = 1e5;
    //calculate alpha
    alpha = pow(1 + kappa * (1 - sqrt(temp / Tcrit)), 2);
    //Find cubic roots
    nSolnValues = peng_robinson_phase->solveCubic(temp, pres, a_coeff, b_coeff, alpha * a_coeff, Vroot);
    EXPECT_NEAR(expected_result[0], Vroot[0], 1.e-6);
    EXPECT_NEAR(nSolnValues, 1 , 1.e-6);

    // Obtain pressure using EoS and compare against the given pressure value
    set_r(1.0);
    rho = test_phase->meanMolecularWeight()/Vroot[0];
    peng_robinson_phase->setState_TR(temp, rho);
    p = peng_robinson_phase->pressure();
    EXPECT_NEAR(p, pres, 1);

    //Liquid phase, supercritical -> nSolnValues = -1
    temp = 300;
    pres = 80e5;
    //calculate alpha
    alpha = pow(1 + kappa * (1 - sqrt(temp / Tcrit)), 2);
    //Find cubic roots
    nSolnValues = peng_robinson_phase->solveCubic(temp, pres, a_coeff, b_coeff, alpha * a_coeff, Vroot);
    EXPECT_NEAR(expected_result[1], Vroot[0], 1.e-6);
    EXPECT_NEAR(nSolnValues, -1, 1.e-6);

    // Obtain pressure using EoS and compare against the given pressure value
    set_r(1.0);
    rho = test_phase->meanMolecularWeight()/Vroot[0];
    peng_robinson_phase->setState_TR(temp, rho);
    p = peng_robinson_phase->pressure();
    EXPECT_NEAR(p, pres, 1);

    //Near critical point -> nSolnValues = -2
    temp = Tcrit;
    //calculate alpha
    alpha = pow(1 + kappa * (1 - sqrt(temp / Tcrit)), 2);
    //Find cubic roots
    nSolnValues = peng_robinson_phase->solveCubic(Tcrit, pCrit, a_coeff, b_coeff, alpha * a_coeff, Vroot);
    EXPECT_NEAR(expected_result[2], Vroot[0], 1.e-6);
    EXPECT_NEAR(nSolnValues, -2, 1.e-6);

    // Obtain pressure using EoS and compare against the given pressure value
    set_r(1.0);
    rho = test_phase->meanMolecularWeight()/Vroot[0];
    peng_robinson_phase->setState_TR(Tcrit, rho);
    p = peng_robinson_phase->pressure();
    EXPECT_NEAR(p, pCrit, 1);
}
#endif
};
