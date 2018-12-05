#include "gtest/gtest.h"
#include "cantera/base/AnyMap.h"

using namespace Cantera;

TEST(AnyValue, is_copyable) {
    AnyMap map1, map2;
    map1["key"] = "1";
    map2["key"] = "2";
    AnyValue value1 = map1["key"];
    AnyValue value2 = map2["key"];
    value2 = value1;
    EXPECT_EQ(value1.asString(), "1");
    EXPECT_EQ(value2.asString(), "1");
}

TEST(AnyValue, is_moveable) {
    AnyMap map1, map2;
    map1["key"] = "1";
    map2["key"] = "2";
    AnyValue value1 = map1["key"];
    AnyValue value2 = map2["key"];
    value2 = std::move(value1);
    EXPECT_EQ(value2.asString(), "1");
}

TEST(AnyMap, paths) {
    AnyMap m;
    m["simple"] = "qux";
    m["compound/first"] = "bar";
    m["compound/second"] = "baz";
    m["several"]["layers"]["deep"] = "foo";

    EXPECT_TRUE(m.hasKey("simple"));
    EXPECT_TRUE(m.hasKey("compound"));
    EXPECT_TRUE(m.hasKey("compound/first"));
    EXPECT_TRUE(m["compound"].hasKey("second"));

    EXPECT_EQ(m["compound/first"].asString(), "bar");
    EXPECT_EQ(m["compound/second"].as<std::string>(), "baz");
    EXPECT_EQ(m["several/layers/deep"].asString(), "foo");
    EXPECT_THROW(m.at("missing"), std::exception);
    EXPECT_FALSE(m.hasKey("missing"));
    EXPECT_FALSE(m.hasKey("compound/missing"));
    EXPECT_THROW(m["missing"].asString(), std::exception);
}

TEST(AnyMap, map_conversion) {
    AnyMap m;
    m["compound/first"] = "bar";
    m["compound/second"] = "baz";

    auto x = m["compound"].asMap<std::string>();
    EXPECT_EQ(x.size(), (size_t) 2);
    EXPECT_EQ(x["first"], "bar");
    EXPECT_EQ(x["second"], "baz");

    std::map<std::string, double> zz{{"a", 9.0}, {"b", 13.5}};
    m["foo"] = zz;
    EXPECT_TRUE(m.hasKey("foo/a"));
    EXPECT_DOUBLE_EQ(m["foo"]["b"].asDouble(), 13.5);

    EXPECT_THROW(m["foo"]["a"].asString(), CanteraError);
    EXPECT_THROW(m["foo"]["b"].asVector<double>(), CanteraError);

    m["qux"]["c"] = 3;
    m["qux"]["d"] = 3.14;
    auto y = m["qux"].asMap<double>();
    EXPECT_DOUBLE_EQ(y["c"], 3.0);
    EXPECT_DOUBLE_EQ(y["d"], 3.14);
    EXPECT_THROW(m["qux"].asMap<long int>(), CanteraError);
}

TEST(AnyMap, nested)
{
    AnyMap m;
    AnyMap b;
    b["foo"] = "bar";
    b["qux"] = 0.5;
    m["baz"] = b;
    b["late"] = "nope"; // Should have added a copy
    EXPECT_EQ(m["baz/foo"].asString(), "bar");
    EXPECT_EQ(m["baz"]["qux"].asDouble(), 0.5);
    EXPECT_FALSE(m.hasKey("baz/late"));
}

TEST(AnyMap, vector)
{
    AnyMap m;
    vector_fp yy{9.6, 14.4, 28.8};
    m["nested"]["item"] = yy;
    vector_fp& yref = m["nested/item"].asVector<double>();
    yref.push_back(-1);
    EXPECT_EQ(yy.size(), (size_t) 3); // Should have added a copy
    // Should have modified the copy in the map
    EXPECT_EQ(m["nested/item"].asVector<double>().size(), (size_t) 4);
}

TEST(AnyMap, getters_with_defaults)
{
    AnyMap m;
    std::map<std::string, double> zz{{"a", 9.0}, {"b", 13.5}};
    m["foo"] = zz;
    m["foo/c"] = 4;
    m["bar"] = "baz";
    m["qux"] = false;
    EXPECT_FALSE(m.getBool("qux", true));
    EXPECT_TRUE(m.getBool("missing", true));
    EXPECT_EQ(m.getString("missing", "hi"), "hi");
    EXPECT_EQ(m.getString("bar", "hi"), "baz");
    EXPECT_EQ(m.getInt("missing", 3), 3);
    EXPECT_EQ(m.getInt("foo/c", 3), 4);
    EXPECT_EQ(m.getDouble("foo/missing", 3), 3);
    EXPECT_EQ(m.getDouble("foo/b", 3), 13.5);
}

TEST(AnyMap, conversion_to_double)
{
    AnyMap m = AnyMap::fromYamlString(
        "{scalar: 8, list: [7, 4, 1], nested: [[3, 4, 5], [12, 22, 91]]}");
    const AnyMap n = m;
    EXPECT_EQ(m["scalar"].asDouble(), 8);
    EXPECT_EQ(m["list"].asVector<double>()[1], 4);
    EXPECT_EQ(m["nested"].asVector<vector_fp>()[1][2], 91);
    EXPECT_EQ(n.at("scalar").asDouble(), 8);
    EXPECT_EQ(n.at("list").asVector<double>()[0], 7);
    EXPECT_EQ(n.at("nested").asVector<vector_fp>()[0][2], 5);
}

TEST(AnyMap, conversion_to_anyvalue)
{
    AnyMap m = AnyMap::fromYamlString(
        "{floats: [7.5, 40, -3.14], strings: [foo, bar]}");
    const AnyMap n = m;
    EXPECT_EQ(m["floats"].asVector<AnyValue>()[0].asDouble(), 7.5);
    EXPECT_EQ(m["strings"].asVector<AnyValue>()[1].asString(), "bar");
    EXPECT_EQ(n.at("floats").asVector<AnyValue>()[2].asDouble(), -3.14);
    EXPECT_EQ(n.at("strings").asVector<AnyValue>()[0].asString(), "foo");
}


TEST(AnyMap, loadYaml)
{
    AnyMap m = AnyMap::fromYamlString(
        "name: NO2\n"
        "composition: {N: 1, O: 2}\n"
        "thermo:\n"
        "  model: NASA7\n"
        "  reference-pressure: 1 atm\n"
        "  temperature-ranges: [200, 1000, 6000]\n"
        "  data:\n"
        "  - [3.944031200E+00, -1.585429000E-03, 1.665781200E-05, -2.047542600E-08,\n"
        "     7.835056400E-12, 2.896617900E+03, 6.311991700E+00]\n"
        "  - [4.884754200E+00, 2.172395600E-03, -8.280690600E-07, 1.574751000E-10,\n"
        "     -1.051089500E-14, 2.316498300E+03, -1.174169500E-01]\n"
        "transport: # this is a comment\n"
        "  model: gas\n"
        "  geometry: nonlinear\n"
        "  flag: true\n");

    EXPECT_EQ(m["name"].asString(), "NO2");
    EXPECT_EQ(m["composition/N"].asInt(), 1);
    EXPECT_EQ(m["thermo/reference-pressure"].asString(), "1 atm");
    EXPECT_EQ(m["thermo/temperature-ranges"].asVector<long int>()[0], 200);
    EXPECT_EQ(m["transport/geometry"].asString(), "nonlinear");
    EXPECT_TRUE(m["transport/flag"].asBool());
    auto coeffs = m["thermo/data"].asVector<vector_fp>();
    EXPECT_EQ(coeffs.size(), (size_t) 2);
    EXPECT_EQ(coeffs[0].size(), (size_t) 7);
    EXPECT_DOUBLE_EQ(coeffs[1][2], -8.280690600E-07);
}
