#include "cantera/thermo/Species.h"

#include "cantera/thermo/SpeciesThermoInterpType.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include <iostream>
#include <limits>

namespace Cantera {

Species::Species()
    : charge(std::numeric_limits<double>::quiet_NaN())
    , size(std::numeric_limits<double>::quiet_NaN())
    , thermo_(0)
{
}

Species::Species(const std::string& name_, const compositionMap& comp_,
                 SpeciesThermoInterpType* therm, double charge_, double size_)
    : name(name_)
    , composition(comp_)
    , charge(charge_)
    , size(size_)
    , thermo_(therm)
{
}

Species::~Species()
{
    delete thermo_;
}

Species::Species(const Species& other)
    : name(other.name)
    , composition(other.composition)
    , charge(other.charge)
    , size(other.size)
{
    if (other.thermo_) {
        thermo_ = other.thermo_->duplMyselfAsSpeciesThermoInterpType();
    } else {
        thermo_ = 0;
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
    delete thermo_;
    if (other.thermo_) {
        thermo_ = other.thermo_->duplMyselfAsSpeciesThermoInterpType();
    } else {
        thermo_ = 0;
    }
    return *this;
}


const SpeciesThermoInterpType& Species::thermo() const
{
    if (thermo_) {
        return *thermo_;
    } else {
        throw CanteraError("Species::thermo",
                           "No thermo for species " + name);
    }
}

}
