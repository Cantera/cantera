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
    if (other.thermo) {
        thermo.reset(other.thermo->duplMyselfAsSpeciesThermoInterpType());
    }
}

Species& Species::operator=(const Species& other)
{
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
    shared_ptr<Species> s(new Species(name, comp));
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
    std::vector<XML_Node*> species_nodes =
        node.child("speciesData").getChildren("species");

    for (std::vector<XML_Node*>::iterator iter = species_nodes.begin();
         iter != species_nodes.end();
         ++iter)
    {
        all_species.push_back(newSpecies(**iter));
    }
    return all_species;
}

}
