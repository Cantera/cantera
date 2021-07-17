#include "gtest/gtest.h"
#include "cantera/base/Units.h"
#include "cantera/base/AnyMap.h"

using namespace Cantera;

TEST(Units, from_string) {
    EXPECT_EQ(Units("").str(), "1.0");
    EXPECT_EQ(Units("1.").str(), "1.0");
    EXPECT_EQ(Units("2").str(), "2.0");
    EXPECT_EQ(Units("kg").str(), "1.0 kg");
    EXPECT_EQ(Units("1.0 kg^0.5").str(), "1.0 kg^0.5");
    EXPECT_EQ(Units("kg / m^3").str(), "1.0 kg / m^3");
    EXPECT_EQ(Units("0.001 m^3").str(), "0.001 m^3");
}

TEST(Units, from_string_short) {
    EXPECT_EQ(Units("").unit_str(), "1");
    EXPECT_EQ(Units("").unit_str(false), "");
    EXPECT_EQ(Units("1.").unit_str(), "1");
    EXPECT_EQ(Units("kg").unit_str(), "kg");
    EXPECT_EQ(Units("1.0 kg^0.5").unit_str(), "kg^0.5");
    EXPECT_EQ(Units("kg / m^3").unit_str(), "kg / m^3");
    EXPECT_EQ(Units("1 / s").unit_str(), "1 / s");
}

TEST(Units, copy_construct) {
    EXPECT_EQ(Units(Units(1.)).str(), "1.0");
    EXPECT_EQ(Units(Units(1., 1.)).str(), "1.0 kg");
    EXPECT_EQ(Units(Units(1., 2.)).str(), "1.0 kg^2");
    EXPECT_EQ(Units(Units(1., .5)).str(), "1.0 kg^0.5");
    EXPECT_EQ(Units(Units(1., -1.)).str(), "1.0 / kg");
    EXPECT_EQ(Units(Units(1., 1., -3.)).str(), "1.0 kg / m^3");
    EXPECT_EQ(Units(Units(.001, 0., 3.)).str(), "0.001 m^3");
}

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

TEST(Units, with_defaults1) {
    UnitSystem U({"cm", "g", "mol", "atm", "kcal"});
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "m"), 0.01);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "kmol/m^3"), 1000);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "kg/kmol"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "cm^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "Pa"), 101325);
    EXPECT_DOUBLE_EQ(U.convertFrom(101325, "Pa"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "hPa"), 1013.25);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "Pa*m^6/kmol"), 101325*1e-12*1000);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "J"), 4184);
}

TEST(Units, with_defaults2) {
    UnitSystem U({"dyn/cm^2", "K"});
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "Pa"), 0.1);
    EXPECT_DOUBLE_EQ(U.convertFrom(1.0, "Pa"), 10);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "N/m^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(300.0, "K"), 300.0);
}

TEST(Units, with_defaults_map) {
    std::map<std::string, std::string> defaults{
        {"length", "cm"}, {"mass", "g"}, {"quantity", "mol"},
        {"pressure", "atm"}, {"energy", "J"}
    };
    UnitSystem U;
    U.setDefaults(defaults);
    EXPECT_DOUBLE_EQ(U.convertFrom(0.01, "m"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "kmol/m^3"), 1000);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "kg/kmol"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "cm^2"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "Pa"), 101325);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "hPa"), 1013.25);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "Pa*m^6/kmol"), 101325*1e-12*1000);
    EXPECT_DOUBLE_EQ(U.convertTo(1.0, "J/cm^3"), 1.0);
}

TEST(Units, bad_defaults) {
    UnitSystem U;
    std::map<std::string, std::string> bad_key{{"length", "m"}, {"joy", "MJ"}};
    EXPECT_THROW(U.setDefaults(bad_key), CanteraError);
    std::map<std::string, std::string> bad_value{{"length", "m"}, {"time", "J"}};
    EXPECT_THROW(U.setDefaults(bad_value), CanteraError);
}


TEST(Units, activation_energies1) {
    UnitSystem U;
    EXPECT_DOUBLE_EQ(U.convertActivationEnergy(1000, "J/kmol", "J/mol"), 1.0);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergy(100, "K", "K"), 100);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergy(500, "K", "J/kmol"), 500 * GasConstant);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergy(3, "J/mol", "K"), 3000 / GasConstant);
}

TEST(Units, activation_energies2) {
    UnitSystem U;
    U.setDefaultActivationEnergy("cal/mol");
    U.setDefaults({"cm", "g", "J"});
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1000, "cal/mol"), 1000);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1000, "J/kmol"), 4184e3);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyFrom(4184e3, "J/kmol"), 1000);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1000, "K"), 4184e3 / GasConstant);
}

TEST(Units, activation_energies3) {
    UnitSystem U({"cal", "mol"});
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1000, "cal/mol"), 1000);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1000, "J/kmol"), 4184e3);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1000, "K"), 4184e3 / GasConstant);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyFrom(4184e3 / GasConstant, "K"), 1000);
}

TEST(Units, activation_energies4) {
    UnitSystem U;
    U.setDefaultActivationEnergy("K");
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(2000, "K"), 2000);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(2000, "J/kmol"), 2000 * GasConstant);
}

TEST(Units, activation_energies5) {
    UnitSystem U;
    std::map<std::string, std::string> defaults{
        {"quantity", "mol"}, {"energy", "cal"}, {"activation-energy", "K"}
    };
    U.setDefaults(defaults);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(2000, "K"), 2000);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(2000, "J/kmol"), 2000 * GasConstant);
}

TEST(Units, activation_energies6) {
    UnitSystem U;
    std::map<std::string, std::string> defaults{
        {"activation-energy", "eV"}
    };
    U.setDefaults(defaults);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1, "J/kmol"), ElectronCharge * Avogadro);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergyTo(1, "eV"), 1.0);
}

TEST(Units, activation_energies_bad) {
    UnitSystem U;
    EXPECT_THROW(U.convertActivationEnergyTo(1000, "kg"), CanteraError);
    EXPECT_THROW(U.convertActivationEnergyFrom(1000, "K^2"), CanteraError);
}

TEST(Units, from_anymap) {
    AnyMap m = AnyMap::fromYamlString(
        "{p: 12 bar, v: 10, A: 1 cm^2, V: 1,"
        " k1: [5e2, 2, 29000], k2: [1e14, -1, 1300 cal/kmol]}");
    UnitSystem U({"mm", "min", "atm"});
    m.setUnits(U);
    m.applyUnits();
    EXPECT_DOUBLE_EQ(m.convert("p", "Pa"), 12e5);
    EXPECT_DOUBLE_EQ(m.convert("v", "cm/min"), 1.0);
    EXPECT_DOUBLE_EQ(m.convert("A", "mm^2"), 100);
    EXPECT_DOUBLE_EQ(m.convert("V", "m^3"), 1e-9);
    auto k1 = m["k1"].asVector<AnyValue>();
    EXPECT_DOUBLE_EQ(U.convert(k1[0], "m^3/kmol"), 1e-9*5e2);
    EXPECT_DOUBLE_EQ(U.convertActivationEnergy(k1[2], "J/kmol"), 29000);

    // calling applyUnits again should not affect results
    m.setUnits(U);
    m.applyUnits();
    EXPECT_DOUBLE_EQ(m.convert("p", "Pa"), 12e5);
    EXPECT_DOUBLE_EQ(U.convert(k1[0], "m^3/kmol"), 1e-9*5e2);
}

TEST(Units, from_anymap_default) {
    AnyMap m = AnyMap::fromYamlString("{p0: 10 atm, h0: 10 cal/kmol}");
    EXPECT_DOUBLE_EQ(m.convert("p0", "Pa", 999), 10*OneAtm);
    EXPECT_DOUBLE_EQ(m.convert("p1", "Pa", 999), 999);
    EXPECT_DOUBLE_EQ(m.convert("h0", "J/kmol", 999), 41.84);
    EXPECT_DOUBLE_EQ(m.convert("h1", "J/kmol", 999), 999);
}

TEST(Units, to_anymap) {
    UnitSystem U{"kcal", "mol", "cm"};
    AnyMap m;
    m["h0"].setQuantity(90, "kJ/kg");
    m["density"].setQuantity({10, 20}, "kg/m^3");
    m.setUnits(U);
    m.applyUnits();
    EXPECT_DOUBLE_EQ(m["h0"].asDouble(), 90e3 / 4184);
    EXPECT_DOUBLE_EQ(m["density"].asVector<double>()[1], 20.0 * 1e-6);
}

TEST(Units, anymap_quantities) {
    AnyMap m;
    std::vector<AnyValue> values(2);
    values[0].setQuantity(8, "kg/m^3");
    values[1].setQuantity(12, "mg/cl");
    m["a"] = values;
    values.emplace_back("hello");
    m["b"] = values;
    m.applyUnits();
    EXPECT_TRUE(m["a"].is<vector_fp>());
    m.applyUnits();
    EXPECT_TRUE(m["a"].is<vector_fp>());
    auto converted = m["a"].asVector<double>();
    EXPECT_DOUBLE_EQ(converted[0], 8.0);
    EXPECT_DOUBLE_EQ(converted[1], 1.2);
    EXPECT_FALSE(m["b"].is<vector_fp>());
}

TEST(Units, to_anymap_nested) {
    UnitSystem U1{"g", "cm", "mol"};
    UnitSystem U2{"mg", "km"};
    for (int i = 0; i < 4; i++) {
        AnyMap m;
        m["A"].setQuantity(90, "kg/m");
        m["nested"]["B"].setQuantity(12, "m^2");
        auto C = std::vector<AnyMap>(2);
        C[0]["foo"].setQuantity(17, "m^2");
        C[1]["bar"].setQuantity(19, "kmol");
        m["nested"]["C"] = C;
        // Test different orders of setting units, and repeated calls to setUnits
        if (i == 0) {
            m.setUnits(U1);
            m["nested"].as<AnyMap>().setUnits(U2);
        } else if (i == 1) {
            m["nested"].as<AnyMap>().setUnits(U2);
            m.setUnits(U1);
        } else if (i == 2) {
            m.setUnits(U1);
            m["nested"].as<AnyMap>().setUnits(U2);
            m.setUnits(U1);
        } else if (i == 3) {
            m["nested"].as<AnyMap>().setUnits(U2);
            m.setUnits(U1);
            m["nested"].as<AnyMap>().setUnits(U2);
        }
        m.applyUnits();
        EXPECT_DOUBLE_EQ(m["A"].asDouble(), 900) << "case " << i;
        EXPECT_DOUBLE_EQ(m["nested"]["B"].asDouble(), 12e-6) << "case " << i;
        EXPECT_DOUBLE_EQ(m["nested"]["C"].asVector<AnyMap>()[0]["foo"].asDouble(), 17e-6);
        EXPECT_DOUBLE_EQ(m["nested"]["C"].asVector<AnyMap>()[1]["bar"].asDouble(), 19000);
    }
}

TEST(Units, from_yaml) {
    AnyMap m = AnyMap::fromYamlString(
        "units: {length: km}\n"
        "foo:\n"
        "- units: {length: cm}\n" // applies to items in foo
        "- bar: 0.6\n"
        "- baz: 0.2\n"
        "  units: {length: mm}\n" // applies to just this entry (with "baz")
        "spam:\n"
        "- eggs: 3\n"
        "- ham: [0.1, 0.3, 0.5]\n"
    );

    EXPECT_FALSE(m.hasKey("units"));
    EXPECT_DOUBLE_EQ(m.units().convertTo(1, "m"), 1000);
    auto& foo = m["foo"].asVector<AnyMap>();
    EXPECT_DOUBLE_EQ(foo[0].units().convertTo(1, "m"), 0.01);
    EXPECT_DOUBLE_EQ(foo[1].units().convertTo(1, "m"), 0.001);
    EXPECT_DOUBLE_EQ(foo[0].convert("bar", "m"), 0.006);
    auto& spam = m["spam"].asVector<AnyMap>();
    EXPECT_DOUBLE_EQ(spam[0].convert("eggs", "m"), 3000);
    EXPECT_DOUBLE_EQ(spam[1].convertVector("ham", "m")[2], 500);
}

TEST(Units, act_energy_from_yaml) {
    AnyMap m = AnyMap::fromYamlString(
        "units: {energy: J, quantity: mol, activation-energy: K}\n"
        "foo:\n"
        "- units: {quantity: kmol}\n" // applies to items in foo
        "- bar: 0.6\n"
        "- baz: 0.2\n"
        "  units: {energy: kJ}\n" // applies to just this entry (with "baz")
    );
    auto& foo = m["foo"].asVector<AnyMap>();
    EXPECT_DOUBLE_EQ(foo[0].units().convertActivationEnergy(foo[0]["bar"], "K"), 0.6);
    EXPECT_DOUBLE_EQ(foo[1].units().convertActivationEnergy(foo[1]["baz"], "K"), 0.2);
    EXPECT_DOUBLE_EQ(foo[0].convert("bar", "J/mol"), 0.0006);
    EXPECT_DOUBLE_EQ(foo[1].convert("baz", "J/mol"), 0.2);
}
