//! @file Species.h Declaration for class Cantera::Species.

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_SPECIES_H
#define CT_SPECIES_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class SpeciesThermoInterpType;
class TransportData;
class XML_Node;

//! Contains data about a single chemical species
/*!
 *  This class stores the data about a species which may be needed to add it to
 *  a ThermoPhase or Transport object.
 */
class Species
{
public:
    Species();

    //! Constructor
    Species(const std::string& name, const compositionMap& comp,
            double charge=0.0, double size=1.0);

    //! @deprecated To be removed after Cantera 2.3.
    Species(const Species& other);
    //! @deprecated To be removed after Cantera 2.3.
    Species& operator=(const Species& other);
    ~Species();

    //! The name of the species
    std::string name;

    //! The elemental composition of the species. Keys are element names; values
    //! are the corresponding atomicities.
    compositionMap composition;

    //! The electrical charge on the species, in units of the elementary charge.
    double charge;

    //! The effective size of the species. Currently used only for surface
    //! species, where it represents the number of sites occupied.
    double size;

    shared_ptr<TransportData> transport;

    //! Thermodynamic data for the species
    shared_ptr<SpeciesThermoInterpType> thermo;
};

//! Create a new Species object from a 'species' XML_Node.
shared_ptr<Species> newSpecies(const XML_Node& species_node);

//! Generate Species objects for all `<species>` nodes in an XML document.
//!
//! The `<species>` nodes are assumed to be children of the `<speciesData>` node
//! in an XML document with a `<ctml>` root node, as in the case of XML files
//! produced by conversion from CTI files.
//!
//! This function can be used in combination with get_XML_File and
//! get_XML_from_string to get Species objects from either a file or a string,
//! respectively, where the string or file is formatted as either CTI or XML.
std::vector<shared_ptr<Species> > getSpecies(const XML_Node& node);

}

#endif
