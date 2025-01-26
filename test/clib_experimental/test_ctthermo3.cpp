#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera/clib/clib_defs.h"
#include "cantera/clib_experimental/ct3.h"
#include "cantera/clib_experimental/ctsol3.h"
#include "cantera/clib_experimental/ctthermo3.h"

using namespace Cantera;


TEST(ctthermo3, thermo)
{
    int ret;
    int sol = sol3_newSolution("gri30.yaml", "gri30", "none");
    int thermo = sol3_thermo(sol);
    ASSERT_GE(thermo, 0);
    size_t nsp = thermo3_nSpecies(thermo);
    ASSERT_EQ(nsp, 53u);

    ret = thermo3_setTemperature(thermo, 500);
    ASSERT_EQ(ret, 0);
    ret = thermo3_setPressure(thermo, 5 * 101325);
    ASSERT_EQ(ret, 0);
    ret = thermo3_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    ASSERT_EQ(ret, 0);

    ret = thermo3_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    ASSERT_EQ(ret, 0);
    double T = thermo3_temperature(thermo);
    ASSERT_GT(T, 2200);
    ASSERT_LT(T, 2300);

    size_t ns = thermo3_nSpecies(thermo);
    vector<double> work(ns);
    vector<double> X(ns);
    thermo3_getMoleFractions(thermo, ns, X.data());

    thermo3_getPartialMolarEnthalpies(thermo, ns, work.data());
    double prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo3_enthalpy_mole(thermo), 1e-6);

    thermo3_getPartialMolarEntropies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo3_entropy_mole(thermo), 1e-6);

    thermo3_getPartialMolarIntEnergies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo3_intEnergy_mole(thermo), 1e-6);

    thermo3_getPartialMolarCp(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo3_cp_mole(thermo), 1e-6);

    thermo3_getPartialMolarVolumes(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, 1./thermo3_molarDensity(thermo), 1e-6);
}

TEST(ctthermo3, atomicWeights)
{
    int ret;
    int sol = sol3_newSolution("h2o2.yaml", "", "none");
    int thermo = sol3_thermo(sol);

    auto cxx_sol = newSolution("h2o2.yaml", "", "none");
    auto cxx_thermo = cxx_sol->thermo();
    auto cxx_weights = cxx_thermo->atomicWeights();
    ASSERT_EQ(cxx_weights.size(), cxx_thermo->nElements());

    int buflen = thermo3_nElements(thermo);
    ASSERT_EQ(buflen, cxx_thermo->nElements());  // 4
    vector<double> buf(buflen);
    ret = thermo3_getAtomicWeights(thermo, buflen, buf.data());
    ASSERT_EQ(ret, buflen);

    for (size_t i = 0; i < buflen; i++) {
        ASSERT_NEAR(buf[i], cxx_weights[i], 1e-6);
    }
}
