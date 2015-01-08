//! @file TransportData.cpp

#include "cantera/transport/TransportData.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{

TransportData::TransportData(const std::string& name_)
    : name(name_)
{
}

GasTransportData::GasTransportData()
    : diameter(0.0)
    , well_depth(0.0)
    , dipole(0.0)
    , polarizability(0.0)
    , rotational_relaxation(0.0)
    , acentric_factor(0.0)
{
}

GasTransportData::GasTransportData(
        const std::string& name_, const std::string& geometry_,
        double diameter_, double well_depth_, double dipole_,
        double polarizability_, double rot_relax, double acentric)
    : TransportData(name_)
    , geometry(geometry_)
    , diameter(diameter_)
    , well_depth(well_depth_)
    , dipole(dipole_)
    , polarizability(polarizability_)
    , rotational_relaxation(rot_relax)
    , acentric_factor(acentric)
{
}

void GasTransportData::setCustomaryUnits(
        const std::string& name_, const std::string& geometry_,
        double diameter_, double well_depth_, double dipole_,
        double polarizability_, double rot_relax, double acentric)
{
    name = name_;
    geometry = geometry_;
    diameter = 1e-10 * diameter_; // convert from Angstroms to m
    well_depth = Boltzmann * well_depth_; // convert from K to J
    dipole = 1e-21 / lightSpeed * dipole_; // convert from Debye to Coulomb-m
    polarizability = 1e-30 * polarizability_; // convert from Angstroms^3 to m^3
    rotational_relaxation = rot_relax; // pure number
    acentric_factor = acentric; // dimensionless
}

void GasTransportData::validate(const Species& sp)
{
    double nAtoms = 0;
    for (compositionMap::const_iterator iter = sp.composition.begin();
         iter != sp.composition.end();
         ++iter) {
        nAtoms += iter->second;
    }

    if (geometry == "atom") {
        if (nAtoms != 1) {
            throw CanteraError("GasTransportData::validate",
                "invalid geometry for species '" + name + "'. 'atom' specified,"
                " but species contains multiple atoms.");
        }
    } else if (geometry == "linear") {
        if (nAtoms == 1) {
            throw CanteraError("GasTransportData::validate",
                "invalid geometry for species '" + name + "'. 'linear'"
                " specified, but species only contains one atom.");
        }
    } else if (geometry == "nonlinear") {
        if (nAtoms < 3) {
            throw CanteraError("GasTransportData::validate",
                "invalid geometry for species '" + name + "'. 'nonlinear'"
                " specified, but species only contains " + fp2str(nAtoms) +
                " atoms.");
        }
    } else {
        throw CanteraError("GasTransportData::validate",
                           "invalid geometry for species '" + name + "': '" +
                           geometry + "'.");
    }

    if (well_depth < 0.0) {
        throw CanteraError("GasTransportData::validate",
                           "negative well depth for species '" + name + "'.");
    }

    if (diameter <= 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative or zero diameter for species '" + name + "'.");
    }

    if (dipole < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative dipole moment for species '" + name + "'.");
    }

    if (polarizability < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative polarizability for species '" + name + "'.");
    }

    if (rotational_relaxation < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative rotation relaxation number for species '" + name + "'");
    }
}

}
