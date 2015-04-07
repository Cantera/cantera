/**
 *  @file StoichSubstanceSSTP.cpp
 * Definition file for the StoichSubstanceSSTP class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::StoichSubstanceSSTP StoichSubstanceSSTP\endlink)
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 *
 * Copyright 2001 California Institute of Technology
 */
#include "cantera/base/ct_defs.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/StoichSubstanceSSTP.h"
#include "cantera/thermo/SpeciesThermo.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

/*
 * ----  Constructors -------
 */

StoichSubstanceSSTP::StoichSubstanceSSTP():
    SingleSpeciesTP()
{
}

StoichSubstanceSSTP::StoichSubstanceSSTP(const std::string& infile, std::string id_) :
    SingleSpeciesTP()
{
    XML_Node* root = get_XML_File(infile);
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                           "Couldn't find phase name in file:" + id_);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                           "thermo model attribute must be StoichSubstance");
    }
    importPhase(*xphase, this);
}

StoichSubstanceSSTP::StoichSubstanceSSTP(XML_Node& xmlphase, const std::string& id_) :
    SingleSpeciesTP()
{
    if (id_ != "") {
        std::string idxml = xmlphase["id"];
        if (id_ != idxml) {
            throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstanceSSTP::StoichSubstanceSSTP",
                           "thermo model attribute must be StoichSubstance");
    }
    importPhase(xmlphase, this);
}

StoichSubstanceSSTP::StoichSubstanceSSTP(const StoichSubstanceSSTP&  right) :
    SingleSpeciesTP()
{
    *this = operator=(right);
}

StoichSubstanceSSTP&
StoichSubstanceSSTP::operator=(const StoichSubstanceSSTP& right)
{
    if (&right != this) {
        SingleSpeciesTP::operator=(right);
    }
    return *this;
}

ThermoPhase* StoichSubstanceSSTP::duplMyselfAsThermoPhase() const
{
    return new StoichSubstanceSSTP(*this);
}

/*
 * ---- Utilities -----
 */

int StoichSubstanceSSTP::eosType() const
{
    return cStoichSubstance;
}

/*
 * ----- Mechanical Equation of State ------
 */

doublereal StoichSubstanceSSTP::pressure() const
{
    return m_press;
}

void StoichSubstanceSSTP::setPressure(doublereal p)
{
    m_press = p;
}

doublereal StoichSubstanceSSTP::isothermalCompressibility() const
{
    return 0.0;
}

doublereal StoichSubstanceSSTP::thermalExpansionCoeff() const
{
    return 0.0;
}

/*
 * ---- Chemical Potentials and Activities ----
 */

void StoichSubstanceSSTP::
getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal StoichSubstanceSSTP::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal StoichSubstanceSSTP::logStandardConc(size_t k) const
{
    return 0.0;
}

void StoichSubstanceSSTP::
getUnitsStandardConc(doublereal* uA, int k, int sizeUA) const
{
    for (int i = 0; i < 6; i++) {
        uA[i] = 0;
    }
}

/*
 * Properties of the Standard State of the Species in the Solution
 */

void StoichSubstanceSSTP::
getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= GasConstant * temperature();
}

void StoichSubstanceSSTP::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
    doublereal RT = GasConstant * temperature();
    doublereal presCorrect = (m_press - m_p0) /  molarDensity();
    hrt[0] += presCorrect / RT;
}

void StoichSubstanceSSTP::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
}

void StoichSubstanceSSTP::getGibbs_RT(doublereal* grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R[0];
}

void StoichSubstanceSSTP::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R[0];
}

void StoichSubstanceSSTP::getIntEnergy_RT(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    urt[0] = m_h0_RT[0] - m_p0 / molarDensity() / RT;
}

/*
 * ---- Thermodynamic Values for the Species Reference States ----
 */

void StoichSubstanceSSTP::getIntEnergy_RT_ref(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    urt[0] = m_h0_RT[0] - m_p0 / molarDensity() / RT;
}

/*
 * ---- Initialization and Internal functions
 */

void StoichSubstanceSSTP::initThermo()
{
    /*
     * Make sure there is one and only one species in this phase.
     */
    m_kk = nSpecies();
    if (m_kk != 1) {
        throw CanteraError("initThermo",
                           "stoichiometric substances may only contain one species.");
    }
    /*
     * Store the reference pressure in the variables for the class.
     */
    m_p0 = refPressure();

    /*
     * Resize temporary arrays.
     */
    int leng = 1;
    m_h0_RT.resize(leng);
    m_cp0_R.resize(leng);
    m_s0_R.resize(leng);
    /*
     * Call the base class thermo initializer
     */
    SingleSpeciesTP::initThermo();
}

void StoichSubstanceSSTP::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("StoichSubstanceSSTP::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& tnode = phaseNode.child("thermo");
    double dens = ctml::getFloatDefaultUnits(tnode, "density", "kg/m3");
    setDensity(dens);
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

void StoichSubstanceSSTP::setParameters(int n, doublereal* const c)
{
    warn_deprecated("StoichSubstanceSSTP::setParameters");
    doublereal rho = c[0];
    setDensity(rho);
}

void StoichSubstanceSSTP::getParameters(int& n, doublereal* const c) const
{
    warn_deprecated("StoichSubstanceSSTP::getParameters");
    doublereal rho = density();
    n = 1;
    c[0] = rho;
}

void StoichSubstanceSSTP::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "StoichSubstance" && model != "StoichSubstanceSSTP") {
        throw CanteraError("StoichSubstanceSSTP::setParametersFromXML",
                           "thermo model attribute must be StoichSubstance");
    }
    doublereal rho = ctml::getFloat(eosdata, "density", "toSI");
    setDensity(rho);
}

// ------ Methods of class electrodeElectron ------

electrodeElectron::electrodeElectron():
    StoichSubstanceSSTP()
{
}

electrodeElectron::electrodeElectron(const std::string& infile, std::string id_) :
    StoichSubstanceSSTP()
{
    XML_Node* root = get_XML_File(infile);
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("electrodeElectron::electrodeElectron",
                           "Couldn't find phase name in file:" + id_);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "electrodeElectron") {
        throw CanteraError("electrodeElectron::electrodeElectron",
                           "thermo model attribute must be electrodeElectron");
    }
    importPhase(*xphase, this);
}

electrodeElectron::electrodeElectron(XML_Node& xmlphase, const std::string& id_) :
    StoichSubstanceSSTP()
{
    if (id_ != "") {
        std::string idxml = xmlphase["id"];
        if (id_ != idxml) {
            throw CanteraError("electrodeElectron::electrodeElectron",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "electrodeElectron") {
        throw CanteraError("electrodeElectron::electrodeElectron",
                           "thermo model attribute must be electrodeElectron");
    }
    importPhase(xmlphase, this);
}

electrodeElectron::electrodeElectron(const electrodeElectron&  right) :
    StoichSubstanceSSTP()
{
    *this = operator=(right);
}

electrodeElectron&
electrodeElectron::operator=(const electrodeElectron& right)
{
    if (&right != this) {
        StoichSubstanceSSTP::operator=(right);
    }
    return *this;
}

electrodeElectron::~electrodeElectron()
{
}

void electrodeElectron::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "electrodeElectron") {
        throw CanteraError("electrodeElectron::setParametersFromXML",
                           "thermo model attribute must be electrodeElectron");
    }
}

void electrodeElectron::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    doublereal rho = 10.0;
    setDensity(rho);
    SingleSpeciesTP::initThermoXML(phaseNode, id_);
}

void electrodeElectron::setParameters(int n, doublereal* const c)
{
    warn_deprecated("electrodeElectron::setParameters");
    doublereal rho = 10.0;
    setDensity(rho);
}

}
