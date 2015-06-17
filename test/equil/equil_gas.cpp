#include "gtest/gtest.h"

#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/IdealGasPhase.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/base/global.h"
#include "cantera/base/utilities.h"

using namespace Cantera;

bool double_close(double expected, double actual, double tol)
{
    return std::abs(expected-actual) / (std::abs(expected) + tol) < tol;
}

#define EXPECT_CLOSE(a,b,tol) EXPECT_PRED3(double_close, a,b,tol)

class OverconstrainedEquil : public testing::Test
{
public:
    OverconstrainedEquil() {}
    void setup(const std::string& elements="H C O N Ar") {
        XML_Node* phase = get_XML_from_string(
            "ideal_gas(elements='" + elements + "', species='gri30: CH C2H2')");
        gas.reset(newPhase(*phase->findByName("phase")));
        gas->setState_TPX(1000, 1e5, "C2H2:0.9, CH:0.1");
    }

    shared_ptr<ThermoPhase> gas;
};

TEST_F(OverconstrainedEquil, ChemEquil)
{
    setup();
    gas->equilibrate("TP", "element_potential");
    EXPECT_NEAR(gas->moleFraction("C2H2"), 1.0, 1e-10);
    EXPECT_NEAR(gas->moleFraction("CH"), 0.0, 1e-10);
    vector_fp mu(2);
    gas->getChemPotentials(&mu[0]);
    EXPECT_NEAR(2*mu[0], mu[1], 1e-7*std::abs(mu[0]));
}

TEST_F(OverconstrainedEquil, VcsNonideal)
{
    setup();
    gas->equilibrate("TP", "vcs");
    EXPECT_NEAR(gas->moleFraction("C2H2"), 1.0, 1e-10);
    EXPECT_NEAR(gas->moleFraction("CH"), 0.0, 1e-10);
    vector_fp mu(2);
    gas->getChemPotentials(&mu[0]);
    EXPECT_NEAR(2*mu[0], mu[1], 1e-7*std::abs(mu[0]));
}

TEST_F(OverconstrainedEquil, DISABLED_MultiphaseEquil)
{
    setup();
    gas->equilibrate("TP", "gibbs");
    EXPECT_NEAR(gas->moleFraction("C2H2"), 1.0, 1e-10);
    EXPECT_NEAR(gas->moleFraction("CH"), 0.0, 1e-10);
    vector_fp mu(2);
    gas->getChemPotentials(&mu[0]);
    EXPECT_NEAR(2*mu[0], mu[1], 1e-7*std::abs(mu[0]));
}

TEST_F(OverconstrainedEquil, BasisOptimize)
{
    setup();
    MultiPhase mphase;
    mphase.addPhase(gas.get(), 10.0);
    mphase.init();
    int usedZeroedSpecies = 0;
    std::vector<size_t> orderVectorSpecies;
    std::vector<size_t> orderVectorElements;

    bool doFormMatrix = true;
    vector_fp formRxnMatrix;

    size_t nc = BasisOptimize(&usedZeroedSpecies, doFormMatrix, &mphase,
                              orderVectorSpecies, orderVectorElements,
                              formRxnMatrix);
    ASSERT_EQ(1, (int) nc);
}

TEST_F(OverconstrainedEquil, DISABLED_BasisOptimize2)
{
    setup("O H C N Ar");
    MultiPhase mphase;
    mphase.addPhase(gas.get(), 10.0);
    mphase.init();
    int usedZeroedSpecies = 0;
    std::vector<size_t> orderVectorSpecies;
    std::vector<size_t> orderVectorElements;

    bool doFormMatrix = true;
    vector_fp formRxnMatrix;

    size_t nc = BasisOptimize(&usedZeroedSpecies, doFormMatrix, &mphase,
                              orderVectorSpecies, orderVectorElements,
                              formRxnMatrix);
    ASSERT_EQ(1, (int) nc);
}

class GriMatrix : public testing::Test
{
public:
    GriMatrix() : gas("gri30.xml", "gri30") {
        X.resize(gas.nSpecies());
        Yelem.resize(gas.nElements());
    };

    void save_elemental_mole_fractions() {
        for (size_t i = 0; i < gas.nElements(); i++) {
            Yelem[i] = gas.elementalMassFraction(i);
        }
    }

    void check(double T, double P) {
        EXPECT_CLOSE(gas.temperature(), T, 1e-9);
        EXPECT_CLOSE(gas.pressure(), P, 1e-9);

        for (size_t i = 0; i < gas.nElements(); i++) {
            EXPECT_CLOSE(Yelem[i], gas.elementalMassFraction(i), 1e-8);
        }

        vector_fp mu(gas.nSpecies());
        gas.getChemPotentials(&mu[0]);
        double mu_C = mu[gas.speciesIndex("C")];
        double mu_H = mu[gas.speciesIndex("H")];
        double mu_O = mu[gas.speciesIndex("O")];
        double mu_N = mu[gas.speciesIndex("N")];
        double mu_Ar = mu[gas.speciesIndex("AR")];

        gas.getMoleFractions(&X[0]);
        for (size_t k = 0; k < gas.nSpecies(); k++) {
            if (X[k] < 1e-15) {
                continue;
            }
            shared_ptr<Species> s = gas.species(k);
            double muk = mu_C * getValue(s->composition, std::string("C"), 0.0) +
                         mu_H * getValue(s->composition, std::string("H"), 0.0) +
                         mu_O * getValue(s->composition, std::string("O"), 0.0) +
                         mu_N * getValue(s->composition, std::string("N"), 0.0) +
                         mu_Ar * getValue(s->composition, std::string("AR"), 0.0);
            EXPECT_CLOSE(muk, mu[k], 1e-7);
        }
    }

    void check_CH4_N2(const std::string& solver) {
        for (int i = 0; i < 5; i++) {
            double T = 500 + 300 * i;
            gas.setState_TPX(T, OneAtm, "CH4:3, N2:2");
            save_elemental_mole_fractions();
            gas.equilibrate("TP", solver);
            check(T, OneAtm);
        }
    }

    void check_O2_N2(const std::string& solver) {
        for (int i = 0; i < 5; i++) {
            double T = 500 + 300 * i;
            gas.setState_TPX(T, OneAtm, "O2:3, N2:2");
            save_elemental_mole_fractions();
            gas.equilibrate("TP", solver);
            check(T, OneAtm);
        }
    }

    void check_CH4_O2_N2(const std::string& solver) {
        for (int i = 0; i < 6; i++) {
            double T = 500 + 300 * i;
            gas.setState_TPX(T, OneAtm, "CH4:3, O2:3, N2:4");
            save_elemental_mole_fractions();
            gas.equilibrate("TP", solver);
            check(T, OneAtm);
        }
    }

    void check_CH4_O2(const std::string& solver) {
        for (int i = 0; i < 5; i++) {
            compositionMap comp;
            comp["CH4"] = i * 0.6 / 5.0;
            comp["O2"] = 1.0 - i * 0.6 / 5.0;
            comp["N2"] = 0.2;
            for (int j = 0; j < 8; j++) {
                double P = std::pow(10.0, j) * 1e-2;
                for (int k = 0; k < 10; k++) {
                    double T = 300 + 250 * k;
                    gas.setState_TPX(T, P, "CH4:1, O2:1");
                    save_elemental_mole_fractions();
                    gas.equilibrate("TP", solver);
                    check(T, P);
                }
            }
        }
    }

    IdealGasPhase gas;
    vector_fp X;
    vector_fp Yelem;
};

TEST_F(GriMatrix, ChemEquil_CH4_N2) { check_CH4_N2("element_potential"); }
TEST_F(GriMatrix, ChemEquil_O2_N2) { check_O2_N2("element_potential"); }
TEST_F(GriMatrix, ChemEquil_CH4_O2_N2) { check_CH4_O2_N2("element_potential"); }
TEST_F(GriMatrix, ChemEquil_CH4_O2) { check_CH4_O2("element_potential"); }
TEST_F(GriMatrix, MultiPhase_CH4_N2) { check_CH4_N2("gibbs"); }
TEST_F(GriMatrix, MultiPhase_O2_N2) { check_O2_N2("gibbs"); }
TEST_F(GriMatrix, MultiPhase_CH4_O2_N2) { check_CH4_O2_N2("gibbs"); }
TEST_F(GriMatrix, DISABLED_MultiPhase_CH4_O2) { check_CH4_O2("gibbs"); }
TEST_F(GriMatrix, VcsNonideal_CH4_N2) { check_CH4_N2("vcs"); }
TEST_F(GriMatrix, VcsNonideal_O2_N2) { check_O2_N2("vcs"); }
TEST_F(GriMatrix, VcsNonideal_CH4_O2_N2) { check_CH4_O2_N2("vcs"); }
TEST_F(GriMatrix, VcsNonideal_CH4_O2) { check_CH4_O2("vcs"); }

int main(int argc, char** argv)
{
    printf("Running main() from equil_gas.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    appdelete();
    return result;
}
