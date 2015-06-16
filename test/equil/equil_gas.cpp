#include "gtest/gtest.h"

#include "cantera/thermo/ThermoFactory.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/base/global.h"

using namespace Cantera;

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

int main(int argc, char** argv)
{
    printf("Running main() from equil_gas.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    appdelete();
    return result;
}
