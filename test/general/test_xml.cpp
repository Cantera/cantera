#include "gtest/gtest.h"
#include "cantera/base/xml.h"
#include <fstream>

namespace Cantera
{

TEST(XML_Node, read_write)
{
    // Check that we can re-read our own XML output
    XML_Node node1, node2;
    std::stringstream out;
    node1.build("../data/air-no-reactions.xml");
    node1.write(out);

    node2.build(out);
    ASSERT_EQ(node2.name(), "ctml");
    ASSERT_TRUE(node2.hasChild("speciesData"));

    std::vector<XML_Node*> species1, species2;
    species1 = node1.child("speciesData").getChildren("species");
    species2 = node2.child("speciesData").getChildren("species");
    ASSERT_EQ(species1.size(), species2.size());
    for (size_t i = 0; i < species1.size(); i++) {
        ASSERT_EQ(species1[i]->attrib("name"), species2[i]->attrib("name"));
    }
}

}
