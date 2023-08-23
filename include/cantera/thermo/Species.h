//! @file Species.h Declaration for class Cantera::Species.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_SPECIES_H
#define CT_SPECIES_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"

namespace Cantera
{

class SpeciesThermoInterpType;
class TransportData;
class ThermoPhase;

//! Contains data about a single chemical species
/*!
 *  This class stores the data about a species which may be needed to add it to
 *  a ThermoPhase or Transport object.
 */
class Species
{
public:
    Species() = default;

    //! Constructor
    Species(const string& name, const Composition& comp,
            double charge=0.0, double size=1.0);

    //! Species objects are not copyable or assignable
    Species(const Species&) = delete;
    Species& operator=(const Species& other) = delete;
    ~Species() = default;

    AnyMap parameters(const ThermoPhase* phase=0, bool withInput=true) const;

    //! The name of the species
    string name;

    //! The elemental composition of the species. Keys are element names; values
    //! are the corresponding atomicities.
    Composition composition;

    //! The electrical charge on the species, in units of the elementary charge.
    double charge = 0.0;

    //! The effective size of the species. Currently used only for surface
    //! species, where it represents the number of sites occupied.
    double size = 1.0;

    //! The molecular weight [amu] of the species.
    /*!
     * Calculates and sets the molecular weight from the elemental composition of the
     * species and element definitions in Elements.cpp, if the molecular weight is
     * Undef.
     *
     * @since New in version 3.0
     */
    double molecularWeight();

    //! Set the molecular weight of the species.
    /*!
     * Since phases can have custom element weights, the phase will always call this
     * method when a species is added to that phase. The species may also call this
     * method the first time the molecularWeight() method is called if the species has
     * not been added to a phase.
     *
     * @param weight: The weight of this species to assign
     *
     * @since New in version 3.0
     */
    void setMolecularWeight(double weight);

    shared_ptr<TransportData> transport;

    //! Thermodynamic data for the species
    shared_ptr<SpeciesThermoInterpType> thermo;

    //! Input parameters used to define a species, for example from a YAML input file.
    AnyMap input;

protected:

    //! The molecular weight of the species, in atomic mass units. Includes
    //! electron mass for charged species.
    double m_molecularWeight = Undef;
};

//! Create a new Species object from an AnyMap specification
unique_ptr<Species> newSpecies(const AnyMap& node);

//! Generate Species objects for each item (an AnyMap) in `items`.
vector<shared_ptr<Species>> getSpecies(const AnyValue& items);

}

#endif
