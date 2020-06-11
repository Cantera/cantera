//! @file TransportData.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_TRANSPORTDATA_H
#define CT_TRANSPORTDATA_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class Species;
class XML_Node;

//! Base class for transport data for a single species
class TransportData
{
public:
    TransportData() {}
    virtual ~TransportData() {}

    virtual void validate(const Species& species) {}

    //! Input data used for specific models
    AnyMap input;
};

//! Transport data for a single gas-phase species which can be used in
//! mixture-averaged or multicomponent transport models.
class GasTransportData : public TransportData
{
public:
    GasTransportData();

    //! Construct a GasTransportData object using MKS units for all parameters.
    GasTransportData(const std::string& geometry, double diameter,
                     double well_depth, double dipole=0.0,
                     double polarizability=0.0, double rot_relax=0.0,
                     double acentric=0.0, double dispersion=0.0,
                     double quad_polar=0.0);

    //! Set the parameters using "customary" units: diameter in Angstroms, well
    //! depth in Kelvin, dipole in Debye, and polarizability in Angstroms^3.
    //! These are the units used in in CK-style input files.
    void setCustomaryUnits(const std::string& geometry, double diameter,
                           double well_depth, double dipole=0.0,
                           double polarizability=0.0, double rot_relax=0.0,
                           double acentric=0.0, double dispersion=0.0,
                           double quad_polar=0.0);

    //! Check transport data for invalid parameters such as a geometry
    //! inconsistent with the atomic composition, non-positive diameter, or
    //! negative values for well depth, dipole, polarizability, or
    //! rotational relaxation number.
    virtual void validate(const Species& species);

    //! A string specifying the molecular geometry. One of `atom`, `linear`, or
    //! `nonlinear`.
    std::string geometry;

    //! The Lennard-Jones collision diameter [m]
    double diameter;

    //! The Lennard-Jones well depth [J]
    double well_depth;

    //! The permanent dipole moment of the molecule [Coulomb-m]. Default 0.0.
    double dipole;

    //! The polarizability of the molecule [m^3]. Default 0.0.
    double polarizability;

    //! The rotational relaxation number (the number of collisions it takes to
    //! equilibrate the rotational degrees of freedom with the temperature).
    //! Default 0.0.
    double rotational_relaxation;

    //! Pitzer's acentric factor [dimensionless]. Default 0.0.
    double acentric_factor;

    //! dispersion normalized by e^2. [m^5] Default 0.0.
    double dispersion_coefficient;

    //! quadrupole. Default 0.0.
    double quadrupole_polarizability;
};

//! Create a new TransportData object from a 'transport' XML_Node.
//!
//! @deprecated The XML input format is deprecated and will be removed in
//!     Cantera 3.0.
shared_ptr<TransportData> newTransportData(const XML_Node& transport_node);

//! Create a new TransportData object from an AnyMap specification
unique_ptr<TransportData> newTransportData(const AnyMap& node);

}

#endif
