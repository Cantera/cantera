#include "cantera/thermo/Species.h"

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include <iostream>
#include <limits>

namespace Cantera {

Species::Species()
    : charge(0.0)
    , size(0.0)
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

}
