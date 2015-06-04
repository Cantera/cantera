#include "gtest/gtest.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

TEST(parseCompString, space_separated)
{
    compositionMap c = parseCompString("foo:1.0  bar:2   baz:1e-4");
    ASSERT_EQ((size_t) 3, c.size());
    ASSERT_DOUBLE_EQ(1.0, c["foo"]);
    ASSERT_DOUBLE_EQ(2.0, c["bar"]);
    ASSERT_DOUBLE_EQ(1e-4, c["baz"]);
}

TEST(parseCompString, extra_spaces)
{
    compositionMap c = parseCompString("foo: 1.0  bar:  2   baz : 1e-4");
    ASSERT_EQ((size_t) 3, c.size());
    ASSERT_DOUBLE_EQ(1.0, c["foo"]);
    ASSERT_DOUBLE_EQ(2.0, c["bar"]);
    ASSERT_DOUBLE_EQ(1e-4, c["baz"]);
}

TEST(parseCompString, default_values)
{
    std::vector<std::string> x;
    x.push_back("foo");
    x.push_back("bar");
    x.push_back("baz");
    compositionMap c = parseCompString("foo:1.0  baz:2", x);
    ASSERT_EQ((size_t) 3, c.size());
    ASSERT_FALSE(c.find("bar") == c.end());
    ASSERT_DOUBLE_EQ(0.0, c["bar"]);
}

TEST(parseCompString, delimiters)
{
    compositionMap c = parseCompString("\nfoo:1.0\tbar:2;baz:1e-4,qux:-1  ");
    ASSERT_EQ((size_t) 4, c.size());
    ASSERT_DOUBLE_EQ(1.0, c["foo"]);
    ASSERT_DOUBLE_EQ(2.0, c["bar"]);
    ASSERT_DOUBLE_EQ(1e-4, c["baz"]);
    ASSERT_DOUBLE_EQ(-1, c["qux"]);
}

TEST(parseCompString, missing_element)
{
    std::vector<std::string> x;
    x.push_back("foo");
    x.push_back("bar");
    ASSERT_THROW(parseCompString("foo:1.0  bar:2   baz:1e-4", x),
                 CanteraError);
}

TEST(parseCompString, not_a_number)
{
    ASSERT_THROW(parseCompString("foo:1.0  bar:five   baz:1e-4"),
                 CanteraError);
}

TEST(parseCompString, missing_value)
{
    ASSERT_THROW(parseCompString("foo:1.0  bar:   baz:1e-4"),
                 CanteraError);
}

TEST(parseCompString, missing_last_value)
{
    ASSERT_THROW(parseCompString("foo:1.0  bar"),
                 CanteraError);
}

}

int main(int argc, char** argv)
{
    printf("Running main() from string_processing.cpp\n");
    testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    Cantera::appdelete();
    return result;
}
