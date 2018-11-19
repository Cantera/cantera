//! @file Units.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/Units.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/AnyMap.h"

namespace {
using namespace Cantera;

const std::map<std::string, Units> knownUnits{
    {"", Units(1.0)},
    {"1", Units(1.0)},

    // Mass [M]
    {"kg", Units(1.0, 1, 0, 0)},
    {"g", Units(1e-3, 1, 0, 0)},

    // Length [L]
    {"m", Units(1.0, 0, 1, 0)},
    {"micron", Units(1e-6, 0, 1, 0)},
    {"angstrom", Units(1e-10, 0, 1, 0)},
    {"Å", Units(1e-10, 0, 1, 0)},

    // Time [T]
    {"s", Units(1.0, 0, 0, 1)},
    {"min", Units(60, 0, 0, 1)},
    {"hr", Units(3600, 0, 0, 1)},

    // Temperature [K]
    {"K", Units(1.0, 0, 0, 0, 1)},
    {"C", Units(1.0, 0, 0, 0, 1)},

    // Current [A]
    {"A", Units(1.0, 0, 0, 0, 0, 1)},

    // Quantity [Q]
    {"mol", Units(1e-3, 0, 0, 0, 0, 0, 1)},
    {"gmol", Units(1e-3, 0, 0, 0, 0, 0, 1)},
    {"mole", Units(1e-3, 0, 0, 0, 0, 0, 1)},
    {"kmol", Units(1.0, 0, 0, 0, 0, 0, 1)},
    {"kgmol", Units(1.0, 0, 0, 0, 0, 0, 1)},
    {"molec", Units(1.0/Avogadro, 0, 0, 0, 0, 0, 1)},

    // Energy [M*L^2/T^2]
    {"J", Units(1.0, 1, 2, -2)},
    {"cal", Units(4.184, 1, 2, -2)},
    {"erg", Units(1e-7, 1, 2, -2)},
    {"eV", Units(ElectronCharge, 1, 2, -2)},

    // Pressure [M/L/T^2]
    {"Pa", Units(1.0, 1, -1, -2)},
    {"atm", Units(OneAtm, 1, -1, -2)},
    {"bar", Units(1.0e5, 1, -1, -2)},

    // Volume [L^3]
    {"m^3", Units(1.0, 0, 3, 0)},
    {"liter", Units(0.001, 0, 3, 0)},
    {"L", Units(0.001, 0, 3, 0)},
    {"l", Units(0.001, 0, 3, 0)},
    {"cc", Units(1.0e-6, 0, 3, 0)},

    // Other electrical units
    {"ohm", Units(1.0, 1, 2, -3, 0, -2)}, // kg*m^2/s^3/A^2
    {"V", Units(1.0, 1, 2, -3, 0, -1)}, // kg*m^2/s^3/A
    {"coulomb", Units(1.0, 0, 0, 1, 0, 1)}, // A*s
};

const std::map<std::string, double> prefixes{
    {"Y", 1e24},
    {"Z", 1e21},
    {"E", 1e18},
    {"P", 1e15},
    {"T", 1e12},
    {"G", 1e9},
    {"M", 1e6},
    {"k", 1e3},
    {"h", 1e2},
    {"d", 1e-1},
    {"c", 1e-2},
    {"m", 1e-3},
    {"u", 1e-6},
    {"n", 1e-9},
    {"p", 1e-12},
    {"f", 1e-15},
    {"a", 1e-18},
    {"z", 1e-21},
    {"y", 1e-24}
};
}

namespace Cantera
{

Units::Units(double factor, double mass, double length, double time,
             double temperature, double current, double quantity)
    : m_factor(factor)
    , m_mass_dim(mass)
    , m_length_dim(length)
    , m_time_dim(time)
    , m_temperature_dim(temperature)
    , m_current_dim(current)
    , m_quantity_dim(quantity)
{
}

Units::Units(const std::string& name)
    : m_factor(1.0)
    , m_mass_dim(0)
    , m_length_dim(0)
    , m_time_dim(0)
    , m_temperature_dim(0)
    , m_current_dim(0)
    , m_quantity_dim(0)
{
    size_t start = 0;
    while (true) {
        // Split into groups of the form 'unit^exponent'
        size_t stop = name.find_first_of("*/", start);
        size_t carat = name.find('^', start);
        if (carat > stop) {
            // No carat in this group
            carat = npos;
        }
        std::string unit = trimCopy(
            name.substr(start, std::min(carat, stop) - start));

        double exponent = 1.0;
        if (carat != npos) {
            exponent = fpValueCheck(name.substr(carat+1, stop-carat-1));
        }
        if (start != 0 && name[start-1] == '/') {
            // This unit is in the denominator
            exponent = -exponent;
        }

        if (knownUnits.find(unit) != knownUnits.end()) {
            // Incorporate the unit defined by the current group
            *this *= knownUnits.at(unit).pow(exponent);
        } else {
            // See if the unit looks like a prefix + base unit
            std::string prefix = unit.substr(0, 1);
            std::string suffix = unit.substr(1);
            if (prefixes.find(prefix) != prefixes.end() &&
                knownUnits.find(suffix) != knownUnits.end()) {
                Units u = knownUnits.at(suffix);
                u.scale(prefixes.at(prefix));
                *this *= u.pow(exponent);
            } else {
                throw CanteraError("Units::Units(string)",
                    "Unknown unit '{}' in unit string '{}'", unit, name);
            }
        }

        start = stop+1;
        if (stop == npos) {
            break;
        }
    }
}

bool Units::convertible(const Units& other) const
{
    return (m_mass_dim == other.m_mass_dim &&
            m_length_dim == other.m_length_dim &&
            m_time_dim == other.m_time_dim &&
            m_temperature_dim == other.m_temperature_dim &&
            m_current_dim == other.m_current_dim &&
            m_quantity_dim == other.m_quantity_dim);
}

Units& Units::operator*=(const Units& other)
{
    m_factor *= other.m_factor;
    m_mass_dim += other.m_mass_dim;
    m_length_dim += other.m_length_dim;
    m_time_dim += other.m_time_dim;
    m_temperature_dim += other.m_temperature_dim;
    m_current_dim += other.m_current_dim;
    m_quantity_dim += other.m_quantity_dim;
    return *this;
}

Units Units::pow(double exponent) const {
    return Units(std::pow(m_factor, exponent),
                 m_mass_dim * exponent,
                 m_length_dim * exponent,
                 m_time_dim * exponent,
                 m_temperature_dim * exponent,
                 m_current_dim * exponent,
                 m_quantity_dim * exponent);
}

std::string Units::str() const {
    return fmt::format("Units({} kg^{} * m^{} * s^{} * K^{} * A^{} * kmol^{})",
                       m_factor, m_mass_dim, m_length_dim, m_time_dim,
                       m_temperature_dim, m_current_dim, m_quantity_dim);
}


double UnitSystem::convert(double value, const std::string& src,
                           const std::string& dest) const
{
    return convert(value, Units(src), Units(dest));
}

double UnitSystem::convert(double value, const Units& src,
                           const Units& dest) const
{
    if (!src.convertible(dest)) {
        throw CanteraError("UnitSystem::convert",
            "Incompatible units:\n    {} and\n    {}", src.str(), dest.str());
    }
    return value * src.factor() / dest.factor();
}

}
