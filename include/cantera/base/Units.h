/**
 * @file Units.h
 * Header for unit conversion utilities, which are used to translate
 * user input from input files (See \ref inputfiles and
 * class \link Cantera::Units Units\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITS_H
#define CT_UNITS_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

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
    explicit Units(const std::string& name);

    //! Returns `true` if the specified Units are dimensionally consistent
    bool convertible(const Units& other) const;

    //! Return the factor for converting from this unit to Cantera's base
    //! units.
    double factor() const { return m_factor; }

    //! Multiply two Units objects, combining their conversion factors and
    //! dimensions
    Units& operator*=(const Units& other);

    //! Provide a string representation of these Units
    std::string str() const;

private:
    //! Scale the unit by the factor `k`
    void scale(double k) { m_factor *= k; }

    //! Raise these Units to a power, changing both the conversion factor and
    //! the dimensions of these Units.
    Units pow(double expoonent) const;

    double m_factor; //!< conversion factor to Cantera base units
    double m_mass_dim;
    double m_length_dim;
    double m_time_dim;
    double m_temperature_dim;
    double m_current_dim;
    double m_quantity_dim;
    double m_pressure_dim; //!< pseudo-dimension to track explicit pressure units

    friend class UnitSystem;
};


//! Unit conversion utility
/*!
 * Provides functions for converting dimensional values from a given unit system.
 * The main use is for converting values specified in input files to Cantera's
 * native unit system, which is SI units except for the use of kmol as the base
 * unit of quantity, i.e. kilogram, meter, second, kelvin, ampere, and kmol.
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
 * Metric prefixes are recognized for all units, e.g. nm, hPa, mg, EJ, mL, kcal.
 *
 * @ingroup inputfiles
 */
class UnitSystem
{
public:
    //! Create a unit system with the specified default units
    UnitSystem(std::initializer_list<std::string> units={});

    //! Set the default units to convert from when explicit units are not
    //! provided. Defaults can be set for mass, length, time, quantity, and
    //! pressure. Conversion using the pressure unit is done only when the
    //! target units explicitly contain pressure units.
    //!
    //! * To use SI+kmol: `setDefaults({"kg", "m", "s", "kmol"});`
    //! * To use CGS+mol: `setDefaults({"cm", "g", "mol"});`
    void setDefaults(std::initializer_list<std::string> units);

    //! Convert `value` from the units of `src` to the units of `dest`.
    double convert(double value, const std::string& src,
                   const std::string& dest) const;
    double convert(double value, const Units& src, const Units& dest) const;
    //! Convert `value` from this unit system (defined by `setDefaults`) to the
    //! specified units.
    double convert(double value, const std::string& dest) const;
    double convert(double value, const Units& dest) const;

private:
    //! Factor to convert mass from this unit system to kg
    double m_mass_factor;

    //! Factor to convert length from this unit system to meters
    double m_length_factor;

    //! Factor to convert time from this unit system to seconds
    double m_time_factor;

    //! Factor to convert energy from this unit system to Joules
    double m_energy_factor;

    //! Factor to convert pressure from this unit system to Pa
    double m_pressure_factor;

    //! Factor to convert quantity from this unit system to kmol
    double m_quantity_factor;
};

}

#endif
