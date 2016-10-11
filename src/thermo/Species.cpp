// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/Species.h"
#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/transport/TransportData.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include <iostream>
#include <limits>

namespace Cantera {

Species::Species()
    : charge(0.0)
    , size(1.0)
{
}

Species::Species(const std::string& name_, const compositionMap& comp_,
                 double charge_, double size_)
    : name(name_)
    , composition(comp_)
    , charge(charge_)
    , size(size_)
{
}

Species::~Species()
{
}

Species::Species(const Species& other)
    : name(other.name)
    , composition(other.composition)
    , charge(other.charge)
    , size(other.size)
    , transport(other.transport)
{
    warn_deprecated("Species copy constructor",
                    "To be removed after Cantera 2.3.");
    if (other.thermo) {
        thermo.reset(other.thermo->duplMyselfAsSpeciesThermoInterpType());
    }
}

Species& Species::operator=(const Species& other)
{
    warn_deprecated("Species assignment operator",
                    "To be removed after Cantera 2.3.");
    if (this == &other) {
        return *this;
    }
    name = other.name;
    composition = other.composition;
    charge = other.charge;
    size = other.size;
    transport = other.transport;
    if (other.thermo) {
        thermo.reset(other.thermo->duplMyselfAsSpeciesThermoInterpType());
    }
    return *this;
}

shared_ptr<Species> newSpecies(const XML_Node& species_node)
{
    std::string name = species_node["name"];
    compositionMap comp = parseCompString(species_node.child("atomArray").value());
    auto s = make_shared<Species>(name, comp);
    if (species_node.hasChild("charge")) {
        s->charge = getFloat(species_node, "charge");
    }
    if (species_node.hasChild("size")) {
        s->size = getFloat(species_node, "size");
    }
    s->thermo.reset(newSpeciesThermoInterpType(species_node.child("thermo")));

    // Read transport data, if provided
    if (species_node.hasChild("transport")) {
        s->transport = newTransportData(species_node.child("transport"));
        s->transport->validate(*s);
    }

    return s;
}

std::vector<shared_ptr<Species> > getSpecies(const XML_Node& node)
{
    std::vector<shared_ptr<Species> > all_species;
    for (const auto& spnode : node.child("speciesData").getChildren("species")) {
        all_species.push_back(newSpecies(*spnode));
    }
    return all_species;
}

}
