#include "gtest/gtest.h"
#include "cantera/base/Units.h"

using namespace Cantera;

TEST(Units, convert_to_base_units) {
    UnitSystem U;
    EXPECT_DOUBLE_EQ(U.convert(1.0, "Pa", "kg/m/s^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "J", "kg*m^2/s^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "ohm", "kg*m^2/s^3/A^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "V", "kg*m^2/A*s^-3"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "coulomb", "A*s"), 1.0);
}

TEST(Units, notation) {
    UnitSystem U;
    EXPECT_DOUBLE_EQ(U.convert(2.0, "m^2", "m*m"), 2.0);
    EXPECT_DOUBLE_EQ(U.convert(3.0, "", "kg/kg"), 3.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "1/m^2", "m^-2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(4.0, "s^3", "s^5/s^2"), 4.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "kg * m/s ^2", "s^-2*kg*m"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "  kg*m / s^ 2", "s ^-2 * kg*m"), 1.0);
}

TEST(Units, basic_conversions) {
    UnitSystem U;
    EXPECT_DOUBLE_EQ(U.convert(100, "cm", "m"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(2, "kmol", "mol"), 2000);
    EXPECT_DOUBLE_EQ(U.convert(1000, "cal", "J"), 4184);
    EXPECT_DOUBLE_EQ(U.convert(2, "m^3", "l"), 2000);
    EXPECT_DOUBLE_EQ(U.convert(1, "atm", "Pa"), 101325);
}

TEST(Units, prefixes) {
    UnitSystem U;
    EXPECT_DOUBLE_EQ(U.convert(1.0, "MJ", "J"), 1e6);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "nm", "cm"), 1e-7);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "m^2", "cm^2"), 1e4);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "m/s", "km/hr"), 3.6);
}

TEST(Units, with_defaults) {
    UnitSystem U({"cm", "g", "mol", "atm"});
    EXPECT_DOUBLE_EQ(U.convert(1.0, "m"), 0.01);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "kmol/m^3"), 1000);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "kg/kmol"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "cm^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "Pa"), 101325);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "hPa"), 1013.25);
    EXPECT_DOUBLE_EQ(U.convert(1.0, "Pa*m^6/kmol"), 101325*1e-12*1000);
}
