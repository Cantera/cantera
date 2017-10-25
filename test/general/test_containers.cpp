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
