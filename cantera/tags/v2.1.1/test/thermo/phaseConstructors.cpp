#include "gtest/gtest.h"
#include "cantera/thermo/FixedChemPotSSTP.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class FixedChemPotSstpConstructorTest : public testing::Test
{
};

TEST_F(FixedChemPotSstpConstructorTest, fromXML)
{
    ThermoPhase* p = newPhase("../data/LiFixed.xml", "");
    ASSERT_EQ((int) p->nSpecies(), 1);
    double mu;
    p->getChemPotentials(&mu);
    ASSERT_FLOAT_EQ(-2.3e7, mu);
    delete p;
}

TEST_F(FixedChemPotSstpConstructorTest, SimpleConstructor)
{
    FixedChemPotSSTP p("Li", -2.3e7);
    ASSERT_EQ((int) p.nSpecies(), 1);
    double mu;
    p.getChemPotentials(&mu);
    ASSERT_FLOAT_EQ(-2.3e7, mu);
}

} // namespace Cantera
