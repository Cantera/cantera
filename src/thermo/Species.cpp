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

    // Extra data used for some electrolyte species
    if (species_node.hasChild("stoichIsMods")) {
        s->extra["weak_acid_charge"] = getFloat(species_node, "stoichIsMods");
    }

    if (species_node.hasChild("electrolyteSpeciesType")) {
        s->extra["electrolyte_species_type"] = species_node.child("electrolyteSpeciesType").value();
    }

    // Extra data optionally used by LatticePhase
    const XML_Node* stdstate = species_node.findByName("standardState");
    if (stdstate && stdstate->findByName("molarVolume")) {
        s->extra["molar_volume"] = getFloat(*stdstate, "molarVolume", "toSI");
    }

    // Extra data possibly used by IonsFromNeutralVPSSTP
    const XML_Node* thermo = species_node.findByName("thermo");
    if (thermo && thermo->attrib("model") == "IonFromNeutral") {
        if (thermo->hasChild("specialSpecies")) {
            s->extra["special_species"] = true;
        }
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
