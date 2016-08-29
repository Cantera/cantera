//! @file TransportData.cpp

#include "cantera/transport/TransportData.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
using namespace std;
#include <iostream>

namespace Cantera
{
GasTransportData::GasTransportData()
    : diameter(0.0)
    , well_depth(0.0)
    , dipole(0.0)
    , polarizability(0.0)
    , rotational_relaxation(0.0)
    , dispersion(0.0)
    , quadrupole_polarizability(0.0)
    , acentric_factor(0.0)
{
}

GasTransportData::GasTransportData(
        const std::string& geometry_,
        double diameter_, double well_depth_, double dipole_,
        double polarizability_, double rot_relax, double dispersion_, 
        double quadrupole_polarizability_, double acentric)
    : geometry(geometry_)
    , diameter(diameter_)
    , well_depth(well_depth_)
    , dipole(dipole_)
    , polarizability(polarizability_)
    , rotational_relaxation(rot_relax)
    , dispersion(dispersion_)
    , quadrupole_polarizability(quadrupole_polarizability_)
    , acentric_factor(acentric)
{
}

void GasTransportData::setCustomaryUnits(
        const std::string& geometry_,
        double diameter_, double well_depth_, double dipole_,
        double polarizability_, double rot_relax, double dispersion_, 
        double quadrupole_polarizability_, double acentric)
{
    geometry = geometry_;
    diameter = 1e-10 * diameter_; // convert from Angstroms to m
    well_depth = Boltzmann * well_depth_; // convert from K to J
    dipole = 1e-21 / lightSpeed * dipole_; // convert from Debye to Coulomb-m
    polarizability = 1e-30 * polarizability_; // convert from Angstroms^3 to m^3
    rotational_relaxation = rot_relax; // pure number
    acentric_factor = acentric; // dimensionless
    dispersion = dispersion_;
    quadrupole_polarizability = quadrupole_polarizability_;
}

void GasTransportData::validate(const Species& sp)
{
    double nAtoms = 0;
    for (const auto& elem : sp.composition) {
        nAtoms += elem.second;
    }

    // geometry validation
    if (geometry == "atom") {
        if (nAtoms != 1) {
            throw CanteraError("GasTransportData::validate",
                "invalid geometry for species '{}'. 'atom' specified, but "
                "species contains multiple atoms.", sp.name);
        }
    } else if (geometry == "linear") {
        if (nAtoms == 1) {
            throw CanteraError("GasTransportData::validate",
                "invalid geometry for species '{}'. 'linear' specified, but "
                "species only contains one atom.", sp.name);
        }
    } else if (geometry == "nonlinear") {
        if (nAtoms < 3) {
            throw CanteraError("GasTransportData::validate",
                "invalid geometry for species '{}'. 'nonlinear' specified, but "
                "species only contains {} atoms.", sp.name, nAtoms);
        }
    } else if (geometry == "positive ion") {
        //do nothing for now
        
    } else if (geometry == "negative ion") {
        //do nothing for now
        
    } else if (geometry == "Electron") {
        //do nothing for now
        
    } else {
        throw CanteraError("GasTransportData::validate",
            "invalid geometry for species '{}': '{}'.", sp.name, geometry);
    }
    // others validation
    if (well_depth < 0.0) {
        throw CanteraError("GasTransportData::validate",
                           "negative well depth for species '{}'.", sp.name);
    }

//    if (diameter <= 0.0) {
//        throw CanteraError("GasTransportData::validate",
//            "negative or zero diameter for species '{}'.", sp.name);
//    }
    if (diameter < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative diameter for species '{}'.", sp.name);
    }

    if (dipole < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative dipole moment for species '{}'.", sp.name);
    }

    if (polarizability < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative polarizability for species '{}'.", sp.name);
    }

    if (rotational_relaxation < 0.0) {
        throw CanteraError("GasTransportData::validate",
            "negative rotation relaxation number for species '{}'.", sp.name);
    }
}

void setupGasTransportData(GasTransportData& tr, const XML_Node& tr_node)
{
    std::string geometry, dummy;
    getString(tr_node, "geometry", geometry, dummy);

    double diam = getFloat(tr_node, "LJ_diameter");
    double welldepth = getFloat(tr_node, "LJ_welldepth");

    double dipole = 0.0;
    getOptionalFloat(tr_node, "dipoleMoment", dipole);

    double polar = 0.0;
    getOptionalFloat(tr_node, "polarizability", polar);

    double rot = 0.0;
    getOptionalFloat(tr_node, "rotRelax", rot);

    double dispersion = 0.0;
    getOptionalFloat(tr_node, "dispersion_coefficient", dispersion);

    double qua_polar = 0.0;
    getOptionalFloat(tr_node, "quadrupole_polarizability", qua_polar);

    double acentric = 0.0;
    getOptionalFloat(tr_node, "acentric_factor", acentric);

    tr.setCustomaryUnits(geometry, diam, welldepth, dipole, polar,
                         rot, dispersion, qua_polar, acentric);
}

shared_ptr<TransportData> newTransportData(const XML_Node& transport_node)
{
    std::string model = transport_node["model"];
    if (model == "gas_transport") {
        auto tr = make_shared<GasTransportData>();
        setupGasTransportData(*tr, transport_node);
        return tr;
    } else {
        // Transport model not handled here
        return make_shared<TransportData>();
    }
}

}
