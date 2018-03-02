/**
 *  @file StoichSubstance.cpp
 * Definition file for the StoichSubstance class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::StoichSubstance StoichSubstance\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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

void StoichSubstance::getActivityConcentrations(double* c) const
{
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

void StoichSubstance::getStandardChemPotentials(double* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= RT();
}

void StoichSubstance::getEnthalpy_RT(double* hrt) const
{
    getEnthalpy_RT_ref(hrt);
    double presCorrect = (m_press - m_p0) / molarDensity();
    hrt[0] += presCorrect / RT();
}

void StoichSubstance::getEntropy_R(double* sr) const
{
    getEntropy_R_ref(sr);
}

void StoichSubstance::getGibbs_RT(double* grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R;
}

void StoichSubstance::getCp_R(double* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R;
}

void StoichSubstance::getIntEnergy_RT(double* urt) const
{
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Thermodynamic Values for the Species Reference States ----

void StoichSubstance::getIntEnergy_RT_ref(double* urt) const
{
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Initialization and Internal functions

void StoichSubstance::initThermo()
{
    // Make sure there is one and only one species in this phase.
    if (m_kk != 1) {
        throw CanteraError("initThermo",
                           "stoichiometric substances may only contain one species.");
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
    setDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

void StoichSubstance::setParameters(int n, double* const c)
{
    setDensity(c[0]);
}

void StoichSubstance::getParameters(int& n, double* const c) const
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
    setDensity(getFloat(eosdata, "density", "toSI"));
}

}
