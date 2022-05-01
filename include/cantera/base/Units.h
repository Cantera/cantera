/**
 * @file Units.h
 * Header for unit conversion utilities, which are used to translate
 * user input from input files (See \ref inputfiles and
 * class \link Cantera::Units Units\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITS_H
#define CT_UNITS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

class AnyValue;
class AnyMap;

//! A representation of the units associated with a dimensional quantity.
/*!
 * Used for converting quantities between unit systems and checking for
 * dimensional consistency. Units objects are mainly used within UnitSystem
 * class to convert values from a user-specified Unit system to Cantera's
 * base units (SI + kmol).
 */
class Units
{
public:
    //! Create a Units object with the specified dimensions.
    explicit Units(double factor=1.0, double mass=0, double length=0,
                   double time=0, double temperature=0, double current=0,
                   double quantity=0);

    //! Create an object with the specified dimensions
    //! @param units        A string representation of the units. See UnitSystem
    //!                     for a description of the formatting options.
    //! @param force_unity  ensure that conversion factor is equal to one
    explicit Units(const std::string& units, bool force_unity=false);

    //! Returns `true` if the specified Units are dimensionally consistent
    bool convertible(const Units& other) const;

    //! Return the factor for converting from this unit to Cantera's base
    //! units.
    double factor() const { return m_factor; }

    //! Multiply two Units objects, combining their conversion factors and
    //! dimensions
    Units& operator*=(const Units& other);

    //! Provide a string representation of these Units
    //! @param skip_unity  do not print '1' if conversion factor is equal to one
    std::string str(bool skip_unity=true) const;

    //! Raise these Units to a power, changing both the conversion factor and
    //! the dimensions of these Units.
    Units pow(double exponent) const;

    bool operator==(const Units& other) const;

    //! Return dimension of primary unit component
    //! ("mass", "length", "time", "temperature", "current" or "quantity")
    double dimension(const std::string& primary) const;

private:
    //! Scale the unit by the factor `k`
    void scale(double k) { m_factor *= k; }

    double m_factor; //!< conversion factor to Cantera base units
    double m_mass_dim;
    double m_length_dim;
    double m_time_dim;
    double m_temperature_dim;
    double m_current_dim;
    double m_quantity_dim;
    double m_pressure_dim; //!< pseudo-dimension to track explicit pressure units
    double m_energy_dim; //!< pseudo-dimension to track explicit energy units

    friend class UnitSystem;
};


//! Unit aggregation utility
/*!
 *  Provides functions for updating and calculating effective units from a stack
 *  of unit-exponent pairs. Matching units are aggregated, where a standard unit
 *  simplifies access when joining exponents. The utility is used in the context
 *  of effective reaction rate units.
 *
 * @internal Helper utility class.
 *
 * @warning This class is an experimental part of the %Cantera API and
 *    may be changed or removed without notice.
 */
struct UnitStack
{
    UnitStack(const Units& standardUnits) {
        stack.reserve(2); // covers memory requirements for most applications
        stack.emplace_back(standardUnits, 0.);
    }

    //! Alternative constructor allows for direct assignment of vector
    UnitStack(std::initializer_list<std::pair<Units, double>> units)
        : stack(units) {}

    //! Size of UnitStack
    size_t size() const { return stack.size(); }

    //! Get standard unit used by UnitStack
    Units standardUnits() const;

    //! Set standard units
    void setStandardUnits(Units& standardUnits);

    //! Effective exponent of standard unit
    double standardExponent() const;

    //! Join (update) exponent of standard units, where the updated exponent is
    //! the sum of the pre-existing exponent and the exponent passed as the argument.
    void join(double exponent);

    //! Update exponent of item with matching units; if it does not exist,
    //! add unit-exponent pair at end of stack
    void update(const Units& units, double exponent);

    //! Calculate product of units-exponent stack
    Units product() const;

    std::vector<std::pair<Units, double>> stack; //!< Stack uses vector of pairs
};


//! Unit conversion utility
/*!
 * Provides functions for converting dimensional values from a given unit system.
 * The main use is for converting values specified in input files to Cantera's
 * native unit system, which is SI units except for the use of kmol as the base
 * unit of quantity, that is, kilogram, meter, second, kelvin, ampere, and kmol.
 *
 * String representations of units can be written using multiplication,
 * division, and exponentiation. Spaces are ignored. Positive, negative, and
 * decimal exponents are permitted. Examples:
 *
 *     kg*m/s^2
 *     J/kmol
 *     m*s^-2
 *     J/kg/K
 *
 * Metric prefixes are recognized for all units, such as nm, hPa, mg, EJ, mL, kcal.
 *
 * Special functions for converting activation energies allow these values to be
 * expressed as either energy per quantity, energy (for example, eV), or temperature by
 * applying a factor of the Avogadro number or the gas constant where needed.
 *
 * @ingroup inputfiles
 */
class UnitSystem
{
public:
    //! Create a unit system with the specified default units
    UnitSystem(std::initializer_list<std::string> units);

    //! Default constructor for unit system (needed as VS2019 does not
    //! recognize an optional argument with a default value)
    UnitSystem() : UnitSystem({}) {}

    //! Return default units used by the unit system
    std::map<std::string, std::string> defaults() const;

    //! Set the default units to convert from when explicit units are not
    //! provided. Defaults can be set for mass, length, time, quantity, energy,
    //! and pressure. Conversion using the pressure or energy units is done only
    //! when the target units explicitly contain pressure or energy units.
    //!
    //! * To use SI+kmol: `setDefaults({"kg", "m", "s", "Pa", "J", "kmol"});`
    //! * To use CGS+mol: `setDefaults({"cm", "g", "dyn/cm^2", "erg", "mol"});`
    void setDefaults(std::initializer_list<std::string> units);

    //! Set the default units using a map of dimension to unit pairs.
    //!
    //! Defaults for dimensions not specified will be left unchanged. To use
    //! Cantera's default units:
    //! ```
    //! UnitSystem system;
    //! std::map<string, string> defaults{
    //!     {"length", "m"}, {"mass", "kg"}, {"time", "s"},
    //!     {"quantity", "kmol"}, {"pressure", "Pa"}, {"energy", "J"},
    //!     {"activation-energy", "J/kmol"}
    //! };
    //! setDefaults(defaults);
    //! ```
    void setDefaults(const std::map<std::string, std::string>& units);

    //! Set the default units to convert from when using the
    //! `convertActivationEnergy` function.
    void setDefaultActivationEnergy(const std::string& e_units);

    //! Convert `value` from the units of `src` to the units of `dest`.
    double convert(double value, const std::string& src,
                   const std::string& dest) const;
    double convert(double value, const Units& src, const Units& dest) const;

    //! Convert `value` to the specified `dest` units from the appropriate units
    //! for this unit system (defined by `setDefaults`)
    double convertTo(double value, const std::string& dest) const;
    double convertTo(double value, const Units& dest) const;

    //! Convert `value` from the specified `src` units to units appropriate for
    //! this unit system (defined by `setDefaults`)
    double convertFrom(double value, const std::string& src) const;
    double convertFrom(double value, const Units& src) const;

    //! Convert a generic AnyValue node to the units specified in `dest`. If the
    //! input is a double, convert it using the default units. If the input is a
    //! string, treat this as a dimensioned value, such as '988 kg/m^3' and convert
    //! from the specified units.
    double convert(const AnyValue& val, const std::string& dest) const;
    double convert(const AnyValue& val, const Units& dest) const;

    //! Convert an array of AnyValue nodes to the units specified in `dest`. For
    //! each node, if the value is a double, convert it using the default units,
    //! and if it is a string, treat it as a value with the given dimensions.
    vector_fp convert(const std::vector<AnyValue>& vals,
                      const std::string& dest) const;
    vector_fp convert(const std::vector<AnyValue>& vals,
                      const Units& dest) const;

    //! Convert `value` from the units of `src` to the units of `dest`, allowing
    //! for the different dimensions that can be used for activation energies
    double convertActivationEnergy(double value, const std::string& src,
                                   const std::string& dest) const;

    //! Convert `value` to the units specified by `dest` from the default
    //! activation energy units
    double convertActivationEnergyTo(double value, const std::string& dest) const;
    double convertActivationEnergyTo(double value, const Units& dest) const;

    //! Convert `value` from the units specified by `src` to the default
    //! activation energy units
    double convertActivationEnergyFrom(double value, const std::string& src) const;

    //! Convert a generic AnyValue node to the units specified in `dest`. If the
    //! input is a double, convert it using the default units. If the input is a
    //! string, treat this as a dimensioned value, such as '2.7e4 J/kmol', and
    //! convert from the specified units.
    double convertActivationEnergy(const AnyValue& val,
                                   const std::string& dest) const;

    //! Get the changes to the defaults from `other` to this UnitSystem
    AnyMap getDelta(const UnitSystem& other) const;

private:
    //! Factor to convert mass from this unit system to kg
    double m_mass_factor;

    //! Factor to convert length from this unit system to meters
    double m_length_factor;

    //! Factor to convert time from this unit system to seconds
    double m_time_factor;

    //! Factor to convert pressure from this unit system to Pa
    double m_pressure_factor;

    //! Factor to convert energy from this unit system to J
    double m_energy_factor;

    //! Factor to convert activation energy from this unit system to J/kmol
    double m_activation_energy_factor;

    //! Factor to convert quantity from this unit system to kmol
    double m_quantity_factor;

    //! True if activation energy units are set explicitly, rather than as a
    //! combination of energy and quantity units
    bool m_explicit_activation_energy;

    //! Map of dimensions (mass, length, etc.) to names of specified default
    //! units
    std::map<std::string, std::string> m_defaults;
};

}

#endif
