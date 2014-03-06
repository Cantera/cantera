#include "gtest/gtest.h"
#include "cantera/thermo/ThermoPhase.h"
#include <vector>

namespace Cantera
{

class ThermoPhase_Fixture : public testing::Test
{
protected:
    ThermoPhase test_phase;
public:
    ThermoPhase_Fixture() {}

    ~ThermoPhase_Fixture() {}

    void initializeElements()
    {
      test_phase.addElement("A", 1.);
      test_phase.addElement("B", 2.);
      test_phase.addElement("C", 3.);
    }
};

TEST_F(ThermoPhase_Fixture, SetAndGetElementPotentials)
{
  initializeElements();

  // Check that getElementPotentials returns false if no element potentials have been set yet.
  std::vector<double> getLambda(3);
  EXPECT_FALSE(test_phase.getElementPotentials(&getLambda[0]));

  std::vector<double> tooSmall(2);
  EXPECT_THROW(test_phase.setElementPotentials(tooSmall), CanteraError);

  std::vector<double> setLambda(3);
  setLambda[0] = 1.;
  setLambda[1] = 2.;
  setLambda[2] = 3.;
  test_phase.setElementPotentials(setLambda);

  EXPECT_TRUE(test_phase.getElementPotentials(&getLambda[0]));
  EXPECT_DOUBLE_EQ(setLambda[0], getLambda[0]);
  EXPECT_DOUBLE_EQ(setLambda[1], getLambda[1]);
  EXPECT_DOUBLE_EQ(setLambda[2], getLambda[2]);
}

}

