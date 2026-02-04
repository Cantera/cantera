/**
 *  @file StoichSubstance.cpp
 * Definition file for the StoichSubstance class, which represents a fixed-composition
 * incompressible substance (see @ref thermoprops and
 * class @link Cantera::StoichSubstance StoichSubstance@endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"

namespace Cantera
{

// ----  Constructors -------

StoichSubstance::StoichSubstance(const string& infile, const string& id_)
{
    initThermoFile(infile, id_);
}

// ----- Mechanical Equation of State ------

double StoichSubstance::pressure() const
{
    return m_press;
}

void StoichSubstance::setPressure(double p)
{
    m_press = p;
}

double StoichSubstance::isothermalCompressibility() const
{
    return 0.0;
}

double StoichSubstance::thermalExpansionCoeff() const
{
    return 0.0;
}

// ---- Chemical Potentials and Activities ----

Units StoichSubstance::standardConcentrationUnits() const
{
    return Units(1.0);
}

void StoichSubstance::getActivityConcentrations(span<double> c) const
{
    checkArraySize("StoichSubstance::getActivityConcentrations", c.size(), 1);
    c[0] = 1.0;
}

double StoichSubstance::standardConcentration(size_t k) const
{
    return 1.0;
}

double StoichSubstance::logStandardConc(size_t k) const
{
    return 0.0;
}

// Properties of the Standard State of the Species in the Solution

void StoichSubstance::getStandardChemPotentials(span<double> mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= RT();
}

void StoichSubstance::getEnthalpy_RT(span<double> hrt) const
{
    getEnthalpy_RT_ref(hrt);
    double presCorrect = (m_press - m_p0) / molarDensity();
    hrt[0] += presCorrect / RT();
}

void StoichSubstance::getEntropy_R(span<double> sr) const
{
    getEntropy_R_ref(sr);
}

void StoichSubstance::getGibbs_RT(span<double> grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R;
}

void StoichSubstance::getCp_R(span<double> cpr) const
{
    checkArraySize("StoichSubstance::getCp_R", cpr.size(), 1);
    _updateThermo();
    cpr[0] = m_cp0_R;
}

void StoichSubstance::getIntEnergy_RT(span<double> urt) const
{
    checkArraySize("StoichSubstance::getIntEnergy_RT", urt.size(), 1);
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Thermodynamic Values for the Species Reference States ----

void StoichSubstance::getIntEnergy_RT_ref(span<double> urt) const
{
    checkArraySize("StoichSubstance::getIntEnergy_RT_ref", urt.size(), 1);
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Initialization and Internal functions

void StoichSubstance::initThermo()
{
    // Make sure there is one and only one species in this phase.
    if (m_kk != 1) {
        throw CanteraError("StoichSubstance::initThermo",
                           "stoichiometric substances may only contain one species.");
    }

    if (species(0)->input.hasKey("equation-of-state")) {
        auto& eos = species(0)->input["equation-of-state"].getMapWhere(
            "model", "constant-volume");
        if (eos.hasKey("density")) {
            assignDensity(eos.convert("density", "kg/m^3"));
        } else if (eos.hasKey("molar-density")) {
            assignDensity(meanMolecularWeight() *
                            eos.convert("molar-density", "kmol/m^3"));
        } else if (eos.hasKey("molar-volume")) {
            assignDensity(meanMolecularWeight() /
                            eos.convert("molar-volume", "m^3/kmol"));
        } else {
            throw InputFileError("StoichSubstance::initThermo", eos,
                "equation-of-state entry for species '{}' is missing 'density',"
                " 'molar-volume' or 'molar-density' specification",
                speciesName(0));
        }
    } else if (m_input.hasKey("density")) {
        assignDensity(m_input.convert("density", "kg/m^3"));
    }

    // Store the reference pressure in the variables for the class.
    m_p0 = refPressure();

    // Call the base class thermo initializer
    SingleSpeciesTP::initThermo();
}

void StoichSubstance::getSpeciesParameters(const string& name,
                                           AnyMap& speciesNode) const
{
    SingleSpeciesTP::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name, true);
    const auto S = species(k);
    auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
        "model", "constant-volume", true);
    // Output volume information in a form consistent with the input
    if (S->input.hasKey("equation-of-state")) {
        auto& eosIn = S->input["equation-of-state"];
        if (eosIn.hasKey("density")) {
            eosNode["density"].setQuantity(density(), "kg/m^3");
        } else if (eosIn.hasKey("molar-density")) {
            eosNode["molar-density"].setQuantity(density() / meanMolecularWeight(),
                                                 "kmol/m^3");
        } else {
            eosNode["molar-volume"].setQuantity(meanMolecularWeight() / density(),
                                                "m^3/kmol");
        }
    } else {
        eosNode["molar-volume"].setQuantity(meanMolecularWeight() / density(), "m^3/kmol");
    }
}

}
