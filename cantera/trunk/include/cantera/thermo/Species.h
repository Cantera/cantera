//! @file Species.h Declaration for class Cantera::Species.

#ifndef CT_SPECIES_H
#define CT_SPECIES_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class SpeciesThermoInterpType;

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
    /*!
     *  The Species object takes ownership of the SpeciesThermoInterpType
     *  object, if provided.
     */
    Species(const std::string& name, const compositionMap& comp,
            SpeciesThermoInterpType* thermo=0, double charge=0.0,
            double size=0.0);

    Species(const Species& other);
    Species& operator=(const Species& other);
    ~Species();

    //! Access the thermodynamic parameterization for the species
    const SpeciesThermoInterpType& thermo() const;

    //! The name of the species
    std::string name;

    //! The elemental composition of the species. Keys are element names; values
    //! are the corresponding atomicities.
    compositionMap composition;

    //! The electrical charge on the species, in units of the elementary charge.
    double charge;

    //! The effective size [m] of the species
    double size;

protected:
    //! Thermodynamic data for the species
    SpeciesThermoInterpType* thermo_;
};

}

#endif
