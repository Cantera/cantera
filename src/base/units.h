/**
 * @file units.h
 * Header for units conversion utilities, which are used to translate
 * user input from input files (See \ref inputfiles and
 * class \link Cantera::Unit Unit\endlink).
 *
 * This header is included only by file misc.cpp.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITS_H
#define CT_UNITS_H

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include <mutex>

namespace Cantera
{

//! Unit conversion utility
/*!
 * @ingroup inputfiles
 */
class Unit
{
public:
    //! Initialize the static Unit class.
    static Unit* units() {
        std::unique_lock<std::mutex> lock(units_mutex);
        if (!s_u) {
            s_u = new Unit;
        }
        return s_u;
    }

    //! Destroy the static Unit class
    /*!
     * Note this can't be done in a destructor.
     */
    static void deleteUnit() {
        std::unique_lock<std::mutex> lock(units_mutex);
        delete s_u;
        s_u = 0;
    }

    //! Empty Destructor
    virtual ~Unit() {}

    /**
     * Return the multiplier required to convert an activation
     * energy to SI units.
     * @param units_ activation energy units
     */
    doublereal actEnergyToSI(const std::string& units_) {
        if (m_act_u.find(units_) != m_act_u.end()) {
            return m_act_u[units_];
        } else {
            return toSI(units_);
        }
    }

    /**
     * Return the multiplier required to convert a dimensional quantity with
     * units specified by string 'units' to SI units. The list of recognized
     * units is stored as a stl map <string, doublereal>called m_u[] and m_act_u
     * for activity coefficients. These maps are initialized with likely values.
     *
     * @param units_ String containing the units description
     */
    doublereal toSI(const std::string& units_) {
        // if dimensionless, return 1.0
        if (units_ == "") {
            return 1.0;
        }

        doublereal f = 1.0, fctr;
        std::string u = units_, tok, tsub;
        std::string::size_type k;
        char action = '-';

        while (true) {
            // get token consisting of all characters up to the next
            // dash, slash, or the end of the string
            k = u.find_first_of("/-");
            if (k != std::string::npos) {
                tok = u.substr(0,k);
            } else {
                tok = u;
            }
            size_t tsize = tok.size();
            if (tsize == 0) {
                fctr = 1.0;
            } else if (tok[tsize - 1] == '2') {
                tsub = tok.substr(0,tsize-1);
                fctr = m_u[tsub];
                fctr *= fctr;
            } else if (tok[tsize - 1] == '3') {
                tsub = tok.substr(0,tsize-1);
                fctr = m_u[tsub];
                fctr *= fctr*fctr;
            } else if (tok[tsize - 1] == '4') {
                tsub = tok.substr(0,tsize-1);
                fctr = m_u[tsub];
                fctr *= fctr*fctr*fctr;
            } else if (tok[tsize - 1] == '5') {
                tsub = tok.substr(0,tsize-1);
                fctr = m_u[tsub];
                fctr *= fctr*fctr*fctr*fctr;
            } else if (tok[tsize - 1] == '6') {
                tsub = tok.substr(0,tsize-1);
                fctr = m_u[tsub];
                fctr *= fctr*fctr*fctr*fctr*fctr;
            } else {
                tsub = tok;
                fctr = m_u[tok];
            }

            // tok is not one of the entries in map m_u, then
            // m_u[tok] returns 0.0. Check for this.
            if (fctr == 0) {
                throw CanteraError("Unit::toSI", "unknown unit: '{}'", tsub);
            }
            if (action == '-') {
                f *= fctr;
            } else if (action == '/') {
                f /= fctr;
            }
            if (k == std::string::npos) {
                break;
            }
            action = u[k];
            u = u.substr(k+1,u.size());
        }
        return f;
    }

private:
    /// pointer to the single instance of Unit
    static Unit* s_u;

    //! Map between a string and a units double value
    /*!
     *  This map maps the dimension string to the units value adjustment. Example
     *   -  m_u["m"]    = 1.0;
     *   -  m_u["cm"]   = 0.01;
     */
    std::map<std::string, doublereal> m_u;

    //! Map between a string and a units double value for activation energy units
    /*!
     *  This map maps the dimension string to the units value adjustment. Example
     *   -    m_act_u["K"] =  GasConstant;
     */
    std::map<std::string, doublereal> m_act_u;

    //! Decl for static locker for Units singleton
    static std::mutex units_mutex;

    //! Units class constructor, containing the default mappings between
    //! strings and units.
    Unit() {
        // unity
        m_u["1"] = 1.0;

        // length
        m_u["m"] = 1.0;
        m_u["cm"] = 0.01;
        m_u["km"] = 1.0e3;
        m_u["mm"] = 1.0e-3;
        m_u["micron"] = 1.0e-6;
        m_u["nm"] = 1.0e-9;
        m_u["A"] = 1.0e-10;
        m_u["Angstrom"] = 1.0e-10;
        m_u["Angstroms"] = 1.0e-10;

        // energy
        m_u["J"] = 1.0;
        m_u["kJ"] = 1.0e3;
        m_u["cal"] = 4.184;
        m_u["kcal"] = 4184.0;
        m_u["eV"] = Faraday; //1.60217733e-19;

        // resistance
        m_u["ohm"] = 1.0;

        // quantity
        m_u["mol"] = 1.0e-3;
        m_u["gmol"] = 1.0e-3;
        m_u["mole"] = 1.0e-3;
        m_u["kmol"] = 1.0;
        m_u["kgmol"] = 1.0;
        m_u["molec"] = 1.0/Avogadro;

        // temperature
        m_u["K"] = 1.0;
        m_u["C"] = 1.0;
        m_u["Kelvin"] = 1.0;

        // mass
        m_u["gm"] = 1.0e-3;
        m_u["g"] = 1.0e-3;
        m_u["kg"] = 1.0;

        // pressure
        m_u["atm"] = 1.01325e5;
        m_u["bar"] = 1.0e5;
        m_u["Pa"] = 1.0;

        // time
        m_u["s"] = 1.0;
        m_u["min"] = 60.0;
        m_u["hr"] = 3600.0;
        m_u["ms"] = 0.001;

        // electric potential
        m_u["volt"] = 1.0;

        // charge
        m_u["coulomb"] = 1.0;

        /*
        // frequency  - Took frequency out to reevaluate it. Inverse cm is probably the wrong default unit
        m_u["hZ"] = 0.01/(lightSpeed);
        m_u["cm^-1"] = 1.0;
        m_u["m^-1"] = 0.1;
        m_u["cm-1"] = m_u["cm^-1"];
        m_u["m-1"] = m_u["m^-1"];
        m_u["wavenumbers"] = m_u["cm^-1"];
        */

        // viscosity
        m_u["Pa-s"] = 1;
        m_u["poise"] = 0.1;
        m_u["centipoise"] = 0.001;
        m_u["P"] = 0.1;
        m_u["cP"] = 0.001;

        // volume
        m_u["kL"] = 1.0;
        m_u["liter"] = 0.001;
        m_u["L"] = 0.001;
        m_u["l"] = 0.001;
        m_u["mL"] = 1.0e-6;
        m_u["ml"] = 1.0e-6;
        m_u["cc"] = 1.0e-6;

        m_act_u["eV"] = m_u["eV"]; // /m_u["molec"];
        m_act_u["K"] = GasConstant;
        m_act_u["Kelvin"] = GasConstant;
        m_act_u["Dimensionless"] = (GasConstant * 273.15);
    }
};
}

#endif
