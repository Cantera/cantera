#include "gtest/gtest.h"

#include "cantera/base/Solution.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo/ThermoFactory.h"

using namespace Cantera;

class WaterTransportTest : public testing::Test
{
public:
    WaterTransportTest() {
        phase = newPhase("thermo-models.yaml", "liquid-water");
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

class NoTransportTest : public testing::Test
{
public:
    NoTransportTest() {}

    static void SetUpTestCase() {
        soln_ = newSolution("h2o2.yaml", "", "None");
    }

    static void TearDownTestCase() {
        soln_.reset();
    }

protected:
    static shared_ptr<Solution> soln_;
};

shared_ptr<Solution> NoTransportTest::soln_;

TEST_F(NoTransportTest, check_type)
{
    auto tr = soln_->transport();
    ASSERT_EQ(tr->transportType(), "None");
}

TEST_F(NoTransportTest, check_exceptions_scalar)
{
    // scalar quantities
    auto tr = soln_->transport();
    ASSERT_THROW(tr->viscosity(), CanteraError);
    ASSERT_THROW(tr->bulkViscosity(), CanteraError);
    ASSERT_THROW(tr->ionConductivity(), CanteraError);
    ASSERT_THROW(tr->thermalConductivity(), CanteraError);
    ASSERT_THROW(tr->electricalConductivity(), CanteraError);
    ASSERT_THROW(tr->getElectricConduct(), CanteraError);
    ASSERT_THROW(tr->CKMode(), CanteraError);
}

TEST_F(NoTransportTest, check_exceptions_vector)
{
    // vector quantities
    auto tr = soln_->transport();
    vector_fp out(soln_->thermo()->nSpecies());
    ASSERT_THROW(tr->getSpeciesViscosities(out.data()), CanteraError);
    ASSERT_THROW(tr->getSpeciesIonConductivity(out.data()), CanteraError);
    ASSERT_THROW(tr->mobilityRatio(out.data()), CanteraError);
    ASSERT_THROW(tr->getMobilities(out.data()), CanteraError);
    ASSERT_THROW(tr->getFluidMobilities(out.data()), CanteraError);
    ASSERT_THROW(tr->getThermalDiffCoeffs(out.data()), CanteraError);
    ASSERT_THROW(tr->getMixDiffCoeffs(out.data()), CanteraError);
    ASSERT_THROW(tr->getMixDiffCoeffsMole(out.data()), CanteraError);
    ASSERT_THROW(tr->getMixDiffCoeffsMass(out.data()), CanteraError);
}

class DefaultTransportTest : public testing::Test
{
public:
    DefaultTransportTest() {}

    static void SetUpTestCase() {
        soln_ = newSolution("h2o2.yaml");
    }

    static void TearDownTestCase() {
        soln_.reset();
    }

protected:
    static shared_ptr<Solution> soln_;
};

shared_ptr<Solution> DefaultTransportTest::soln_;

TEST_F(DefaultTransportTest, check_type)
{
    auto tr = soln_->transport();
    ASSERT_EQ(tr->transportType(), "Mix");
}

TEST_F(DefaultTransportTest, check_scalar)
{
    // scalar quantities
    auto tr = soln_->transport();
    EXPECT_GE(tr->viscosity(), 0.);
    EXPECT_GE(tr->thermalConductivity(), 0.);
    EXPECT_FALSE(tr->CKMode());
}
