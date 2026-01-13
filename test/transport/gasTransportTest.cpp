// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "gtest/gtest.h"

#include "cantera/core.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/transport/MultiTransport.h"

using namespace Cantera;

class GasTransportTest : public testing::Test
{
public:
    GasTransportTest() {
        nsp = s_thermo->nSpecies();
        s_thermo->setState_TPX(T0, P0, X0.data());
    }

    static void SetUpTestCase() {
        s_soln = newSolution("gri30.yaml", "", "none");
        s_thermo = s_soln->thermo();
        s_mix = newTransport(s_thermo, "mixture-averaged");
        s_multi = std::dynamic_pointer_cast<MultiTransport>(
            newTransport(s_thermo, "multicomponent"));
    }

    auto SetUpFluxes() {
        vector<double> X2(nsp), X3(nsp);
        X2[0] = 0.25; X2[5] = 0.17; X2[14] = 0.15; X2[15] = 0.05; X2[47] = 0.38;
        // sum(X3) == 1.02, but should still lead to net zero mass flux
        X3[0] = 0.27; X3[5] = 0.15; X3[14] = 0.18; X3[15] = 0.06; X3[47] = 0.36;

        vector<double> grad_T(2, 0.0);
        Array2D grad_X(nsp, 2, 0.0);
        for (size_t k = 0; k < nsp; k++) {
            grad_X(k,0) = (X2[k] - X0[k]) / dist;
            grad_X(k,1) = (X3[k] - X0[k]) / dist;
        }

        grad_T[0] = (T2 - T1) / dist;
        grad_T[1] = (T3 - T1) / dist;
        return std::make_pair(grad_T, grad_X);
    }

    static shared_ptr<Solution> s_soln;
    static shared_ptr<ThermoPhase> s_thermo;
    static shared_ptr<Transport> s_mix;
    static shared_ptr<MultiTransport> s_multi;
    size_t nsp;

    double T0 = 1500;
    double P0 = 1e5;
    double T1 = 1300; // used for flux tests and diffusion coefficient tests
    const vector<double> X0 = {
        0.269205, 0.000107082, 1.36377e-09, 4.35475e-10, 4.34036e-06, 0.192249,
        3.59356e-13, 2.78061e-12, 4.7406e-18, 4.12955e-17, 2.58549e-14, 8.96502e-16,
        6.09056e-11, 7.56752e-09, 0.192253, 0.0385036, 1.49596e-08, 2.22378e-08,
        1.43096e-13, 1.45312e-15, 1.96948e-12, 8.41937e-19, 3.18852e-13, 7.93625e-18,
        3.20653e-15, 1.15149e-19, 1.61189e-18, 1.4719e-15, 5.24728e-13, 6.90582e-17,
        6.37248e-12, 5.93728e-11, 2.71219e-09, 2.66645e-06, 6.57142e-11, 9.52453e-08,
        1.26006e-14, 3.49802e-12, 1.19232e-11, 7.17782e-13, 1.85347e-07, 8.25325e-14,
        5.00914e-20, 1.54407e-16, 3.07176e-11, 4.93198e-08, 4.84792e-12, 0.307675,
        0.0, 6.21649e-29, 8.42393e-28, 6.77865e-18, 2.19225e-16
    };
    double dist = 0.1;
    double T2 = 1000;
    double T3 = 1200;
};

shared_ptr<Solution> GasTransportTest::s_soln;
shared_ptr<ThermoPhase> GasTransportTest::s_thermo;
shared_ptr<Transport> GasTransportTest::s_mix;
shared_ptr<MultiTransport> GasTransportTest::s_multi;

TEST_F(GasTransportTest, mixDiffCoeffs)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double Dref[] = {
        1.610684e-03, 2.257324e-03, 6.399886e-04, 4.226893e-04, 6.291336e-04,
        5.542400e-04, 4.202332e-04, 4.176794e-04, 5.919444e-04, 6.798972e-04,
        4.605880e-04, 4.605880e-04, 4.513537e-04, 4.501499e-04, 4.024747e-04,
        3.314755e-04, 3.617492e-04, 3.590862e-04, 3.501083e-04, 3.501083e-04,
        3.517251e-04, 3.522411e-04, 3.490302e-04, 3.460247e-04, 3.460823e-04,
        3.174043e-04, 3.150784e-04, 5.305271e-04, 3.062054e-04, 3.062054e-04,
        5.649743e-04, 6.736707e-04, 6.600519e-04, 5.019766e-04, 4.049890e-04,
        4.119510e-04, 3.700861e-04, 3.351239e-04, 4.185845e-04, 4.067501e-04,
        3.583467e-04, 3.554188e-04, 5.305177e-04, 3.364528e-04, 3.364528e-04,
        3.364528e-04, 3.378655e-04, 3.916823e-04, 4.122727e-04, 2.465436e-04,
        2.455565e-04, 3.049268e-04, 3.037003e-04
    };

    vector<double> mixDiffs(nsp, 0.0);
    s_mix->getMixDiffCoeffs(mixDiffs.data());
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(mixDiffs[k], Dref[k], 1e-9) << k;
    }
}

TEST_F(GasTransportTest, speciesViscosities)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double viscRef[] = {
        2.555827e-05, 3.012044e-05, 7.284512e-05, 6.239122e-05, 7.510484e-05,
        5.323262e-05, 6.336632e-05, 6.432665e-05, 4.461934e-05, 6.571173e-05,
        3.273503e-05, 3.273503e-05, 3.389082e-05, 3.612525e-05, 5.309983e-05,
        5.429542e-05, 4.165683e-05, 4.237417e-05, 4.194489e-05, 4.194489e-05,
        4.327915e-05, 3.541213e-05, 3.611814e-05, 3.681062e-05, 3.795644e-05,
        3.356043e-05, 3.413748e-05, 1.285435e-04, 4.239037e-05, 4.239037e-05,
        4.818437e-05, 7.599594e-05, 7.850542e-05, 4.717766e-05, 5.229843e-05,
        5.589102e-05, 6.635946e-05, 5.291370e-05, 5.949548e-05, 4.769577e-05,
        3.789946e-05, 3.859977e-05, 1.285498e-04, 5.231643e-05, 5.231643e-05,
        5.231643e-05, 5.169996e-05, 5.400395e-05, 7.252699e-05, 3.016719e-05,
        3.051801e-05, 4.289560e-05, 4.339494e-05
    };

    vector<double> specVisc(nsp, 0.0);
    s_mix->getSpeciesViscosities(specVisc.data());
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(specVisc[k], viscRef[k], 1e-9) << k;
    }
}

TEST_F(GasTransportTest, speciesMobilities)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double mobRef[] = {
        1.133773e-02, 1.586484e-02, 4.502536e-03, 2.971472e-03, 4.426166e-03,
        3.878312e-03, 2.954904e-03, 2.936948e-03, 4.165152e-03, 4.783312e-03,
        3.236468e-03, 3.236468e-03, 3.171577e-03, 3.162222e-03, 2.828705e-03,
        2.324660e-03, 2.528445e-03, 2.509832e-03, 2.445929e-03, 2.445929e-03,
        2.458890e-03, 2.472322e-03, 2.449785e-03, 2.428690e-03, 2.426302e-03,
        2.226214e-03, 2.209901e-03, 3.727740e-03, 2.141911e-03, 2.141911e-03,
        3.975376e-03, 4.739516e-03, 4.640102e-03, 3.501594e-03, 2.849627e-03,
        2.896544e-03, 2.598048e-03, 2.351319e-03, 2.942800e-03, 2.861817e-03,
        2.502430e-03, 2.481983e-03, 3.727674e-03, 2.360642e-03, 2.360642e-03,
        2.360642e-03, 2.370553e-03, 2.751910e-03, 2.897403e-03, 1.728760e-03,
        1.721839e-03, 2.132969e-03, 2.124390e-03
     };

    vector<double> mob(nsp, 0.0);
    s_thermo->setState_TP(T1, P0);
    s_mix->getMobilities(mob.data());
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(mob[k], mobRef[k], 1e-7) << k;
    }
}

TEST_F(GasTransportTest, mixtureViscosity)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double viscRef[] = {
        1.982378e-05, 2.362091e-05, 2.716680e-05, 3.051153e-05, 3.369042e-05,
        3.672940e-05, 3.964816e-05, 4.246202e-05, 4.518316e-05, 4.782149e-05
    };
    for (size_t i = 0; i < 10; i++) {
        double T = 400. + 100. * i;
        s_thermo->setState_TP(T, P0);
        EXPECT_NEAR(s_mix->viscosity(), viscRef[i], 1e-9) << T;
    }
}

TEST_F(GasTransportTest, thermalDiffCoeffs)
{
    // regression test based on case from legacy multiGasTransport.cpp
    const double thermalDiffRef[] = {
        -1.749421e-06, -5.795523e-10, -1.036581e-14,  2.976174e-15, -2.877580e-11,
        -8.681002e-07,  2.942992e-18,  2.437224e-17, -3.075226e-22, -7.242761e-22,
        -1.227900e-19, -4.193045e-21, -2.609672e-16, -3.238708e-14,  7.585519e-07,
         7.037944e-07,  7.190017e-14,  1.204366e-13,  8.710654e-19,  8.966277e-21,
         1.419255e-17, -2.828291e-22,  1.092672e-18,  2.067526e-23,  1.607623e-20,
         -4.376819e-22, -2.041513e-22,  1.158367e-20,  8.472501e-18,  1.273332e-21,
         -4.408996e-17, -4.943295e-16, -2.097570e-14, -1.207527e-11,  3.627187e-16,
          5.204528e-13,  2.546007e-19,  6.566709e-17,  7.782264e-17,  2.117999e-18,
          6.059605e-13,  3.289296e-19,  5.920676e-22,  2.653811e-21,  5.506497e-16,
          8.840154e-13,  8.195568e-17,  1.155794e-06,  2.282194e-23,  7.479437e-23,
         -1.400893e-22,  1.968555e-22,  4.414863e-21
    };

    const double thermalDiffMixRef[] = {
        -1.97226421719082e-06, -6.04696504247212e-10, -9.90802378755412e-15,
        3.14480069520489e-15, -2.96155643142378e-11, -6.82205028945992e-07,
        2.92077699110767e-18, 2.49573967024886e-17, -3.02014062767274e-23,
        -3.43728222424916e-22, -8.19843353222891e-20, -2.84275400736814e-21,
        -1.64208862608778e-16, -1.80369284235688e-14, 8.24896551635327e-07,
        7.04151357884373e-07, 9.34047159936078e-14, 1.54627741091289e-13,
        1.09208834358759e-18, 1.10900054077964e-20, 1.68248638493913e-17,
        4.05329445771533e-24, 1.7718934165487e-18, 5.01179394798185e-23,
        2.16210022278212e-20, 9.85061024456701e-25, 1.50682299818537e-23,
        1.35163391781345e-20, 9.31081556209504e-18, 1.22537421912008e-21,
        -3.60737893615715e-17, -4.86699553776885e-16, -2.10048949693204e-14,
        -8.80399546006993e-12, 3.98149356084823e-16, 5.99167005188409e-13,
        2.51475683431014e-19, 6.8424918923382e-17, 8.00157517467182e-17,
        2.85203772258056e-18, 9.4009635515829e-13, 4.74141387733766e-19,
        4.60157223656553e-25, 2.88317189129165e-21, 5.73575815122357e-16,
        9.20926220995577e-13, 8.61597096394322e-17, 1.12606178968437e-06,
        0, 2.00216138684567e-33, 2.79494070813552e-32,
        1.26003408367763e-22, 4.26137699636566e-21
    };

    vector<double> thermalDiff(nsp, -1);
    s_mix->getThermalDiffCoeffs(thermalDiff.data());
    for (size_t k = 0; k < nsp; k++) {
        double tol = std::max(1e-16, 1e-5 * fabs(thermalDiffMixRef[k]));
        EXPECT_NEAR(thermalDiff[k], thermalDiffMixRef[k], tol) << k;
    }

    s_multi->getThermalDiffCoeffs(thermalDiff.data());
    for (size_t k = 0; k < nsp; k++) {
        double tol = std::max(1e-16, 1e-5 * fabs(thermalDiffRef[k]));
        EXPECT_NEAR(thermalDiff[k], thermalDiffRef[k], tol) << k;
    }
}

TEST_F(GasTransportTest, thermalConductivity_mix)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double condRef[] = {
        6.441981e-02, 7.678686e-02, 8.886259e-02, 1.008372e-01, 1.127636e-01,
        1.246465e-01, 1.364744e-01, 1.482310e-01, 1.598998e-01, 1.714656e-01
    };
    for (size_t i = 0; i < 10; i++) {
        double T = 400. + 100. * i;
        s_thermo->setState_TP(T, P0);
        EXPECT_NEAR(s_mix->thermalConductivity(), condRef[i], 1e-6) << T;
    }
}

TEST_F(GasTransportTest, thermalConductivity_multi)
{
    // regression test based on case from legacy multiGasTransport.cpp
    const double condRef[] = {
        6.336474e-02, 7.586649e-02, 8.779171e-02, 9.964744e-02, 1.116838e-01,
        1.237784e-01, 1.353296e-01, 1.472053e-01, 1.589833e-01, 1.706627e-01
    };
    for (size_t i = 0; i < 10; i++) {
        double T = 400. + 100. * i;
        s_thermo->setState_TP(T, P0);
        EXPECT_NEAR(s_multi->thermalConductivity(), condRef[i], 1e-6) << T;
    }
}

TEST_F(GasTransportTest, binaryDiffCoeffs)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double BdiffRef_H2[] = {
        1.687881e-03, 2.572028e-03, 1.268285e-03, 9.517971e-04, 1.264072e-03,
        1.131811e-03, 9.509353e-04, 9.501239e-04, 1.083528e-03, 1.284426e-03,
        8.691328e-04, 8.691328e-04, 8.654639e-04, 8.775252e-04, 9.071478e-04,
        8.058646e-04, 8.099309e-04, 8.090473e-04, 7.956857e-04, 7.956857e-04,
        8.008189e-04, 7.516299e-04, 7.505447e-04, 7.495389e-04, 7.586945e-04,
        6.959224e-04, 6.951653e-04, 1.275518e-03, 7.237840e-04, 7.237840e-04,
        1.072375e-03, 1.319044e-03, 1.314123e-03, 1.031986e-03, 8.882841e-04,
        9.135785e-04, 8.863859e-04, 7.935054e-04, 9.365126e-04, 8.730697e-04,
        7.929474e-04, 7.919572e-04, 1.275515e-03, 7.939043e-04, 7.939043e-04,
        7.939043e-04, 7.943305e-04, 9.156211e-04, 9.669845e-04, 5.725485e-04,
        5.722560e-04, 7.233961e-04, 7.230257e-04
    };
    Array2D Bdiff(nsp, nsp, 0.0);
    s_thermo->setState_TP(T1, P0);
    size_t kH2 = s_thermo->speciesIndex("H2");
    s_mix->getBinaryDiffCoeffs(nsp, &Bdiff(0,0));
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(Bdiff(k, kH2), BdiffRef_H2[k], 1e-8) << k;
        EXPECT_DOUBLE_EQ(Bdiff(kH2, k), Bdiff(k, kH2)) << k;
    }
}

TEST_F(GasTransportTest, getSpeciesFluxes_mix)
{
    // regression test based on case from legacy mixGasTransport.cpp
    const double refFlux0[] = {
        4.9516e-06, 1.7828e-08, 1.0344e-12, 4.3946e-13, 3.4409e-09, 1.8677e-05,
        3.7203e-16, 2.9490e-15, 2.5002e-21, 2.7051e-20, 1.2434e-17, 4.3113e-19,
        3.0779e-14, 4.0689e-12, 3.8686e-05, -1.0935e-05, 1.1697e-11, 1.7864e-11,
        1.1587e-16, 1.1766e-18, 1.6550e-15, 5.5560e-22, 2.1695e-16, 5.5619e-21,
        2.3287e-18, 7.9692e-23, 1.1461e-21, 2.3779e-18, 5.0601e-16, 6.6695e-20,
        3.7438e-15, 4.4451e-14, 2.1221e-12, 1.6817e-09, 5.7739e-14, 8.7914e-11,
        1.6039e-17, 3.8660e-15, 1.1552e-14, 5.6775e-16, 1.3363e-10, 6.1232e-17,
        8.0932e-23, 1.6746e-19, 3.3314e-14, 5.3489e-11, 5.1555e-15, -5.1403e-05,
        0.0, 4.9978e-32, 6.9043e-31, 6.6666e-21, 2.1979e-19
    };
    const double refFlux1[] = {
        4.7572e-07, 1.788e-08, 1.0449e-12, 4.4619e-13, 3.4766e-09, 3.4831e-05,
        3.7776e-16, 2.9947e-15, 2.5277e-21, 2.7311e-20, 1.2609e-17, 4.3721e-19,
        3.1221e-14, 4.1275e-12, 1.6651e-05, -2.072e-05, 1.1907e-11, 1.8186e-11,
        1.1801e-16, 1.1984e-18, 1.6854e-15, 5.6577e-22, 2.2095e-16, 5.6656e-21,
        2.3722e-18, 8.1308e-23, 1.1695e-21, 2.4071e-18, 5.1666e-16, 6.7997e-20,
        3.7869e-15, 4.4882e-14, 2.1431e-12, 1.7036e-09, 5.866e-14, 8.9294e-11,
        1.6319e-17, 3.9403e-15, 1.1731e-14, 5.7677e-16, 1.3605e-10, 6.2349e-17,
        8.1925e-23, 1.7067e-19, 3.3952e-14, 5.4514e-11, 5.2538e-15, -3.1262e-05,
        0.0, 5.1271e-32, 7.0836e-31, 6.8075e-21, 2.2445e-19
    };

    auto [grad_T, grad_X] = SetUpFluxes();
    Array2D fluxes(nsp, 2, 0.0);
    s_thermo->setState_TPX(T1, P0, X0.data());
    s_mix->getSpeciesFluxes(2, grad_T.data(), nsp, &grad_X(0, 0), nsp, &fluxes(0, 0));

    double netFlux0 = 0.0, netFlux1 = 0.0;
    for (size_t k = 0; k < nsp; k++) {
        double tol = std::max(1e-14, 1e-4 * fabs(refFlux0[k]));
        EXPECT_NEAR(fluxes(k, 0), refFlux0[k], tol) << k;
        tol = std::max(1e-14, 1e-4 * fabs(refFlux1[k]));
        EXPECT_NEAR(fluxes(k, 1), refFlux1[k], tol) << k;
        netFlux0 += fluxes(k, 0);
        netFlux1 += fluxes(k, 1);
    }
    EXPECT_NEAR(netFlux0, 0.0, 1e-19);
    EXPECT_NEAR(netFlux1, 0.0, 1e-19);
}

TEST_F(GasTransportTest, getSpeciesFluxes_multi)
{
    // regression test based on case from legacy multiGasTransport.cpp
    const double refFlux0[] = {
        1.128904e-06, 1.674139e-08, 1.018666e-12, 4.476504e-13, 3.399815e-09,
        1.747460e-05, 3.797483e-16, 3.011704e-15, 0.000000e+00, 2.639948e-20,
        1.232246e-17, 4.270911e-19, 3.058024e-14, 4.044714e-12, 4.232824e-05,
       -9.596638e-06, 1.193098e-11, 1.823651e-11, 1.186113e-16, 1.203028e-18,
        1.695087e-15, 0.000000e+00, 2.211883e-16, 0.000000e+00, 2.382811e-18,
        0.000000e+00, 0.000000e+00, 2.413676e-18, 5.262102e-16, 6.925377e-20,
        3.683223e-15, 4.367280e-14, 2.090469e-12, 1.673417e-09, 5.884649e-14,
        8.947782e-11, 1.660660e-17, 4.017245e-15, 1.177065e-14, 5.761799e-16,
        1.359860e-10, 6.239818e-17, 0.000000e+00, 1.757578e-19, 3.459268e-14,
        5.554134e-11, 5.347453e-15,-5.135724e-05, 0.000000e+00, 0.000000e+00,
        0.000000e+00, 0.000000e+00, 2.291112e-19
    };
    const double refFlux1[] = {
       -1.010313e-06, 1.741584e-08, 1.022405e-12, 4.382750e-13, 3.404361e-09,
        3.132854e-05, 3.711876e-16, 2.942438e-15, 0.000000e+00, 2.667083e-20,
        1.231510e-17, 4.269575e-19, 3.050538e-14, 4.033802e-12, 1.379673e-05,
       -2.196399e-05, 1.168527e-11, 1.784817e-11, 1.160699e-16, 1.178189e-18,
        1.654891e-15, 0.000000e+00, 2.165422e-16, 0.000000e+00, 2.327262e-18,
        0.000000e+00, 0.000000e+00, 2.361602e-18, 5.079424e-16, 6.684925e-20,
        3.699479e-15, 4.391219e-14, 2.098664e-12, 1.672448e-09, 5.759425e-14,
        8.766254e-11, 1.606061e-17, 3.876206e-15, 1.152131e-14, 5.657122e-16,
        1.334565e-10, 6.117334e-17, 0.000000e+00, 1.685016e-19, 3.339743e-14,
        5.362242e-11, 5.166923e-15,-2.217377e-05, 0.000000e+00, 0.000000e+00,
        0.000000e+00, 0.000000e+00, 2.207645e-19
    };

    auto [grad_T, grad_X] = SetUpFluxes();
    Array2D fluxes(nsp, 2, 0.0);
    s_thermo->setState_TPX(T1, P0, X0.data());
    s_multi->getSpeciesFluxes(2, grad_T.data(), nsp, &grad_X(0, 0), nsp, &fluxes(0, 0));

    double netFlux0 = 0.0, netFlux1 = 0.0;
    for (size_t k = 0; k < nsp; k++) {
        double tol = std::max(1e-14, 1e-4 * fabs(refFlux0[k]));
        EXPECT_NEAR(fluxes(k, 0), refFlux0[k], tol) << k;
        tol = std::max(1e-14, 1e-4 * fabs(refFlux1[k]));
        EXPECT_NEAR(fluxes(k, 1), refFlux1[k], tol) << k;
        netFlux0 += fluxes(k, 0);
        netFlux1 += fluxes(k, 1);
    }
    EXPECT_NEAR(netFlux0, 0.0, 1e-19);
    EXPECT_NEAR(netFlux1, 0.0, 1e-19);
}

TEST_F(GasTransportTest, multicomponentDiffusionCoefficients)
{
    const double D_H2_X_ref[] = {
        0.000000e+00, 2.069166e-02, 1.533272e-03, 7.984950e-04, 1.450324e-03,
        1.363695e-03, 7.766991e-04, 7.561488e-04, 1.968227e-03, 1.851414e-03,
        1.681716e-03, 1.681716e-03, 1.576561e-03, 1.486114e-03, 8.951919e-04,
        5.900643e-04, 8.577657e-04, 8.314705e-04, 8.055783e-04, 8.055783e-04,
        7.831462e-04, 9.755145e-04, 9.407016e-04, 9.084341e-04, 8.794914e-04,
        8.454421e-04, 8.193710e-04, 6.691079e-04, 6.070352e-04, 6.070352e-04,
        1.707391e-03, 1.629845e-03, 1.536190e-03, 1.424771e-03, 8.652237e-04,
        8.420368e-04, 5.745784e-04, 5.889965e-04, 8.196214e-04, 9.534937e-04,
        9.135834e-04, 8.834312e-04, 6.690525e-04, 6.009906e-04, 6.009906e-04,
        6.009906e-04, 6.137995e-04, 8.959026e-04, 6.575936e-04, 5.802367e-04,
        5.681411e-04, 5.942564e-04, 5.820537e-04
    };
    const double D_X_H2_ref[] = {
        0.000000e+00, 1.737090e-02, 4.930150e-03, 3.253679e-03, 4.846513e-03,
        3.623689e-03, 3.235537e-03, 3.215876e-03, 4.560725e-03, 5.237593e-03,
        3.543842e-03, 3.543842e-03, 3.472788e-03, 3.462544e-03, 2.547296e-03,
        2.446430e-03, 2.768576e-03, 2.748195e-03, 2.678224e-03, 2.678224e-03,
        2.692416e-03, 2.707123e-03, 2.682446e-03, 2.659347e-03, 2.656733e-03,
        2.437642e-03, 2.419780e-03, 4.081771e-03, 2.345333e-03, 2.345333e-03,
        4.352925e-03, 5.189638e-03, 5.080782e-03, 3.834139e-03, 3.120262e-03,
        3.171634e-03, 2.844790e-03, 2.574628e-03, 3.222283e-03, 3.133609e-03,
        2.740090e-03, 2.717702e-03, 4.081698e-03, 2.584837e-03, 2.584837e-03,
        2.584836e-03, 2.595690e-03, 2.178005e-03, 3.172576e-03, 1.892944e-03,
        1.885366e-03, 2.335541e-03, 2.326147e-03
    };

    Array2D multiDiff(nsp, nsp);
    s_thermo->setState_TPX(T1, P0, X0.data());
    s_multi->getMultiDiffCoeffs(nsp, &multiDiff(0,0));
    size_t kH2 = s_thermo->speciesIndex("H2");
    for (size_t k = 0; k < nsp; k++) {
        EXPECT_NEAR(multiDiff(k, kH2), D_X_H2_ref[k], 1e-8) << k;
        EXPECT_NEAR(multiDiff(kH2, k), D_H2_X_ref[k], 1e-8) << k;
    }
}

TEST_F(GasTransportTest, getFluxes_multi)
{
    vector<double> fluxS(nsp), fluxMass(nsp), fluxMole(nsp);
    vector<double> state2(s_thermo->stateSize());
    vector<double> state3(s_thermo->stateSize());
    vector<double> X1(nsp), X2(nsp), X3(nsp), grad_X(nsp);

    s_thermo->setState_TPX(T2, P0,
        "H2:0.25, H:0.0001, H2O:0.17, CO:0.15, CO2:0.05, NO:0.001, N2: 0.38");
    s_thermo->saveState(state2);
    s_thermo->getMoleFractions(X2);
    s_thermo->setState_TPX(T3, P0, "H2:0.27, H2O:0.18, CO:0.13, CO2:0.04, N2: 0.38");
    s_thermo->saveState(state3);
    s_thermo->getMoleFractions(X3);
    double grad_T = (T3 - T2) / dist;
    for (size_t k = 0; k < nsp; k++) {
        X1[k] = 0.5 * (X2[k] + X3[k]);
        grad_X[k] = (X3[k] - X2[k]) / dist;
    }

    s_thermo->setState_TPX(0.5 * (T2 + T3), P0, X1.data());
    s_multi->getSpeciesFluxes(1, &grad_T, nsp, grad_X.data(), nsp, fluxS.data());
    s_multi->getMassFluxes(state2.data(), state3.data(), dist, fluxMass.data());
    s_multi->getMolarFluxes(state2.data(), state3.data(), dist, fluxMole.data());

    double netFlux = 0.0;
    for (size_t k = 0; k < nsp; k++) {
        double Wk = s_thermo->molecularWeight(k);
        double tol = std::max(1e-14, 1e-4 * fabs(fluxS[k]));
        EXPECT_NEAR(fluxMass[k], fluxS[k], tol) << k;
        EXPECT_NEAR(fluxMole[k] * Wk, fluxS[k], tol) << k;
        netFlux += fluxMass[k];
    }
    EXPECT_NEAR(netFlux, 0.0, 1e-19);
}
