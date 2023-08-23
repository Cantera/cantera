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

//! Base class for transport data for a single species
class TransportData
{
public:
    TransportData() = default;
    virtual ~TransportData() = default;

    virtual void validate(const Species& species) {}

    //! Return the parameters such that an identical species transport object
    //! could be reconstructed using the newTransportData() function. Behavior
    //! specific to derived classes is handled by the getParameters() method.
    //! @param withInput  If true, include additional input data fields associated
    //!   with the object, such as user-defined fields from a YAML input file, as
    //!   stored in the #input attribute.
    AnyMap parameters(bool withInput) const;

    //! Input data used for specific models
    AnyMap input;

protected:
    //! Store the parameters needed to reconstruct a TransportData object. Does
    //! not include user-defined fields available in #input.
    virtual void getParameters(AnyMap& transportNode) const;
};

//! Transport data for a single gas-phase species which can be used in
//! mixture-averaged or multicomponent transport models.
class GasTransportData : public TransportData
{
public:
    GasTransportData() = default;

    //! Construct a GasTransportData object using MKS units for all parameters.
    GasTransportData(const string& geometry, double diameter,
                     double well_depth, double dipole=0.0,
                     double polarizability=0.0, double rot_relax=0.0,
                     double acentric=0.0, double dispersion=0.0,
                     double quad_polar=0.0);

    //! Set the parameters using "customary" units: diameter in Angstroms, well
    //! depth in Kelvin, dipole in Debye, and polarizability in Angstroms^3.
    //! These are the units used in in CK-style input files.
    void setCustomaryUnits(const string& geometry, double diameter,
                           double well_depth, double dipole=0.0,
                           double polarizability=0.0, double rot_relax=0.0,
                           double acentric=0.0, double dispersion=0.0,
                           double quad_polar=0.0);

    //! Check transport data for invalid parameters such as a geometry
    //! inconsistent with the atomic composition, non-positive diameter, or
    //! negative values for well depth, dipole, polarizability, or
    //! rotational relaxation number.
    void validate(const Species& species) override;

    void getParameters(AnyMap& transportNode) const override;

    //! A string specifying the molecular geometry. One of `atom`, `linear`, or
    //! `nonlinear`.
    string geometry;

    //! The Lennard-Jones collision diameter [m]
    double diameter = 0.0;

    //! The Lennard-Jones well depth [J]
    double well_depth = 0.0;

    //! The permanent dipole moment of the molecule [Coulomb-m]. Default 0.0.
    double dipole = 0.0;

    //! The polarizability of the molecule [m^3]. Default 0.0.
    double polarizability = 0.0;

    //! The rotational relaxation number (the number of collisions it takes to
    //! equilibrate the rotational degrees of freedom with the temperature).
    //! Default 0.0.
    double rotational_relaxation = 0.0;

    //! Pitzer's acentric factor [dimensionless]. Default 0.0.
    double acentric_factor = 0.0;

    //! dispersion normalized by e^2. [m^5] Default 0.0.
    double dispersion_coefficient = 0.0;

    //! quadrupole. Default 0.0.
    double quadrupole_polarizability = 0.0;
};

//! Create a new TransportData object from an AnyMap specification
unique_ptr<TransportData> newTransportData(const AnyMap& node);

}

#endif
