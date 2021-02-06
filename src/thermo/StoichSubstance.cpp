/**
 *  @file StoichSubstance.cpp
 * Definition file for the StoichSubstance class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::StoichSubstance StoichSubstance\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/StoichSubstance.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/ctml.h"

namespace Cantera
{

// ----  Constructors -------

StoichSubstance::StoichSubstance(const std::string& infile, const std::string& id_)
{
    initThermoFile(infile, id_);
}

StoichSubstance::StoichSubstance(XML_Node& xmlphase, const std::string& id_)
{
    importPhase(xmlphase, this);
}

// ----- Mechanical Equation of State ------

doublereal StoichSubstance::pressure() const
{
    return m_press;
}

void StoichSubstance::setPressure(doublereal p)
{
    m_press = p;
}

doublereal StoichSubstance::isothermalCompressibility() const
{
    return 0.0;
}

doublereal StoichSubstance::thermalExpansionCoeff() const
{
    return 0.0;
}

// ---- Chemical Potentials and Activities ----

Units StoichSubstance::standardConcentrationUnits() const
{
    return Units(1.0);
}

void StoichSubstance::getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal StoichSubstance::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal StoichSubstance::logStandardConc(size_t k) const
{
    return 0.0;
}

// Properties of the Standard State of the Species in the Solution

void StoichSubstance::getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= RT();
}

void StoichSubstance::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
    doublereal presCorrect = (m_press - m_p0) / molarDensity();
    hrt[0] += presCorrect / RT();
}

void StoichSubstance::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
}

void StoichSubstance::getGibbs_RT(doublereal* grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R;
}

void StoichSubstance::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R;
}

void StoichSubstance::getIntEnergy_RT(doublereal* urt) const
{
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Thermodynamic Values for the Species Reference States ----

void StoichSubstance::getIntEnergy_RT_ref(doublereal* urt) const
{
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

void StoichSubstance::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("StoichSubstance::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    std::string model = tnode["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstance::initThermoXML",
                           "thermo model attribute must be StoichSubstance");
    }
    double dens = getFloat(tnode, "density", "toSI");
    assignDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

void StoichSubstance::setParameters(int n, doublereal* const c)
{
    assignDensity(c[0]);
}

void StoichSubstance::getParameters(int& n, doublereal* const c) const
{
    n = 1;
    c[0] = density();
}

void StoichSubstance::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstance::setParametersFromXML",
                           "thermo model attribute must be StoichSubstance");
    }
    assignDensity(getFloat(eosdata, "density", "toSI"));
}

}
