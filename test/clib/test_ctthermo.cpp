#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <fstream>

#include "cantera/core.h"
#include "cantera_clib/ct.h"
#include "cantera_clib/ctsol.h"
#include "cantera_clib/ctthermo.h"

using namespace Cantera;


TEST(ctthermo, thermo)
{
    int32_t ret;
    int32_t sol = sol_newSolution("gri30.yaml", "gri30", "none");
    int32_t thermo = sol_thermo(sol);
    ASSERT_GE(thermo, 0);
    size_t nsp = thermo_nSpecies(thermo);
    ASSERT_EQ(nsp, 53u);

    ret = thermo_setTemperature(thermo, 500);
    ASSERT_EQ(ret, 0);
    ret = thermo_setPressure(thermo, 5 * 101325);
    ASSERT_EQ(ret, 0);
    ret = thermo_setMoleFractionsByName(thermo, "CH4:1.0, O2:2.0, N2:7.52");
    ASSERT_EQ(ret, 0);

    ret = thermo_equilibrate(thermo, "HP", "auto", 1e-9, 50000, 1000, 0);
    ASSERT_EQ(ret, 0);
    double T = thermo_temperature(thermo);
    ASSERT_GT(T, 2200);
    ASSERT_LT(T, 2300);

    size_t ns = thermo_nSpecies(thermo);
    vector<double> work(ns);
    vector<double> X(ns);
    thermo_getMoleFractions(thermo, ns, X.data());

    thermo_getPartialMolarEnthalpies(thermo, ns, work.data());
    double prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_enthalpy_mole(thermo), 1e-6);

    thermo_getPartialMolarEntropies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_entropy_mole(thermo), 1e-6);

    thermo_getPartialMolarIntEnergies(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_intEnergy_mole(thermo), 1e-6);

    thermo_getPartialMolarCp(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, thermo_cp_mole(thermo), 1e-6);

    thermo_getPartialMolarVolumes(thermo, ns, work.data());
    prod = std::inner_product(X.begin(), X.end(), work.begin(), 0.0);
    ASSERT_NEAR(prod, 1./thermo_molarDensity(thermo), 1e-6);
}

TEST(ctthermo, atomicWeights)
{
    int32_t ret;
    int32_t sol = sol_newSolution("h2o2.yaml", "", "none");
    int32_t thermo = sol_thermo(sol);

    auto cxx_sol = newSolution("h2o2.yaml", "", "none");
    auto cxx_thermo = cxx_sol->thermo();
    auto cxx_weights = cxx_thermo->atomicWeights();
    ASSERT_EQ(cxx_weights.size(), cxx_thermo->nElements());

    int32_t buflen = thermo_nElements(thermo);
    ASSERT_EQ(buflen, static_cast<int32_t>(cxx_thermo->nElements()));  // 4
    vector<double> buf(buflen);
    ret = thermo_atomicWeights(thermo, buflen, buf.data());
    ASSERT_EQ(ret, buflen);

    for (int32_t i = 0; i < buflen; i++) {
        ASSERT_NEAR(buf[i], cxx_weights[i], 1e-6);
    }
}
