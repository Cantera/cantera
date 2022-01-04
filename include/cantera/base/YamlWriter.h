//! @file YamlWriter.h Declaration for class Cantera::YamlWriter.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_YAMLWRITER_H
#define CT_YAMLWRITER_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/AnyMap.h"
#include "cantera/base/Units.h"

namespace Cantera
{

class Solution;
class ThermoPhase;
class Kinetics;
class Transport;

//! A class for generating full YAML input files from multiple data sources
class YamlWriter
{
public:
    YamlWriter();

    //! Include top-level information used in YAML header block
    void setHeader(const AnyMap& header);

    //! Include a phase definition for the specified Solution object
    void addPhase(shared_ptr<Solution> soln, bool includeAdjacent=true);

    //! Include a phase definition using the specified ThermoPhase, (optional)
    //! Kinetics, and (optional) Transport objects
    void addPhase(shared_ptr<ThermoPhase> thermo, shared_ptr<Kinetics> kin={},
                  shared_ptr<Transport> tran={});

    //! Return a YAML string that contains the definitions for the added phases,
    //! species, and reactions
    std::string toYamlString() const;

    //! Write the definitions for the added phases, species and reactions to
    //! the specified file.
    void toYamlFile(const std::string& filename) const;

    //! For output floating point values, set the maximum number of digits to
    //! the right of the decimal point. The default is 15 digits.
    void setPrecision(long int n) {
        m_float_precision = n;
    }

    //! By default user-defined data present in the input is preserved on
    //! output. This method can be used to skip output of user-defined data
    //! fields which are not directly used by Cantera.
    void skipUserDefined(bool skip=true) {
        m_skip_user_defined = skip;
    }

    //! Set the units to be used in the output file. Dimensions not specified
    //! will use Cantera's defaults.
    //! @param units  A map where keys are dimensions (mass, length, time,
    //!     quantity, pressure, energy, activation-energy) and the values are
    //!     corresponding units supported by the UnitSystem class.
    void setUnits(const std::map<std::string, std::string>& units={});

    //! Set the units to be used in the output file. Dimensions not specified
    //! will use Cantera's defaults.
    //! @param units  A UnitSystem object specifying dimensions (mass, length, time,
    //!     quantity, pressure, energy, activation-energy).
    void setUnitSystem(const UnitSystem& units=UnitSystem());

protected:
    //! Top-level information used in YAML header block
    AnyMap m_header;

    std::vector<shared_ptr<Solution>> m_phases;

    //! @see setPrecision()
    long int m_float_precision;

    //! @see skipUserDefined()
    bool m_skip_user_defined;

    //! Top-level units directive for the output file. Defaults to Cantera's
    //! native SI+kmol system.
    UnitSystem m_output_units;
};

}

#endif
