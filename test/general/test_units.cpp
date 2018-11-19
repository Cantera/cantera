#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/base/AnyMap.h"

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

TEST(Units, activation_energies) {
    UnitSystem U;
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(1000, "J/kmol", "J/mol"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(100, "K", "K"), 100);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(500, "K", "J/kmol"), 500 * GasConstant);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(3, "J/mol", "K"), 3000 / GasConstant);

    U.setDefaults({"cm", "g"});
    U.setDefaultMolarEnergy("cal/mol");
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(1000, "cal/mol"), 1000);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(1000, "J/kmol"), 4184e3);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(1000, "K"), 4184e3 / GasConstant);

    U.setDefaultMolarEnergy("K");
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(2000, "K"), 2000);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(2000, "J/kmol"), 2000 * GasConstant);
}

TEST(Units, from_anymap) {
    AnyMap m = AnyMap::fromYamlString(
        "{p: 12 bar, v: 10, A: 1 cm^2, V: 1,"
        " k1: [5e2, 2, 29000], k2: [1e14, -1, 1300 cal/kmol]}");
    UnitSystem U({"mm", "min", "atm"});
    EXPECT_DOUBLE_EQ(U.convert(m["p"], "Pa"), 12e5);
    EXPECT_DOUBLE_EQ(U.convert(m["v"], "cm/min"), 1.0);
    EXPECT_DOUBLE_EQ(U.convert(m["A"], "mm^2"), 100);
    EXPECT_DOUBLE_EQ(U.convert(m["V"], "m^3"), 1e-9);
    auto k1 = m["k1"].asVector<AnyValue>();
    EXPECT_DOUBLE_EQ(U.convert(k1[0], "m^3/kmol"), 1e-9*5e2);
    EXPECT_DOUBLE_EQ(U.convertMolarEnergy(k1[2], "J/kmol"), 29000);
}
