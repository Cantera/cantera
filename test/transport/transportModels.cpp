#include "gtest/gtest.h"

#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace Cantera;

class WaterTransportTest : public testing::Test
{
public:
    WaterTransportTest() {
        phase = newPhase("liquid-water.xml");
        tran = newDefaultTransportMgr(phase);
    }

    void check_viscosity(double T, double P, double mu_expected) {
        phase->setState_TP(T + 273.15, P);
        EXPECT_NEAR(tran->viscosity(), mu_expected, 2e-7);
    }

    void check_thermal_conductivity(double T, double P, double lambda_expected) {
        phase->setState_TP(T + 273.15, P);
        EXPECT_NEAR(tran->thermalConductivity(), lambda_expected, 3e-4);
    }

    ThermoPhase* phase;
    Transport* tran;
};

// Reference values taken from tables in the Sengers and Watson paper
// (doi:10.1063/1.555763) which is the source of the interpolating equations.
TEST_F(WaterTransportTest, viscosity)
{
    check_viscosity(25, 1e5, 890.5e-6);
    check_viscosity(100, 5e5, 281.9e-6);
    check_viscosity(100, 1e7, 284.5e-6);
    check_viscosity(250, 5e6, 106.4e-6);
    check_viscosity(250, 5e7, 117.5e-6);
    check_viscosity(350, 1.75e7, 67.0e-6);
    check_viscosity(400, 1.5e7, 24.93e-6);
}

TEST_F(WaterTransportTest, thermal_conductivity)
{
    check_thermal_conductivity(25, 1e5, 0.6072);
    check_thermal_conductivity(100, 5e5, 0.6793);
    check_thermal_conductivity(100, 1e7, 0.6845);
    check_thermal_conductivity(250, 5e6, 0.6227);
    check_thermal_conductivity(250, 5e7, 0.6721);
    check_thermal_conductivity(350, 1.75e7, 0.4523);
    check_thermal_conductivity(400, 1.5e7, 0.08068);
}
