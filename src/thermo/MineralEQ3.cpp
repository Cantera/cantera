/**
 *  @file MineralEQ3.cpp
 *   Definition file for the MineralEQ3 class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::MineralEQ3 MineralEQ3\endlink)
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
#include "cantera/thermo/MineralEQ3.h"
#include "cantera/thermo/SpeciesThermo.h"

#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/MineralEQ3.h"

using namespace std;

namespace Cantera
{

/*
 * ----  Constructors -------
 */

MineralEQ3::MineralEQ3():
    StoichSubstanceSSTP()
{
}

MineralEQ3::MineralEQ3(const std::string& infile, std::string id_) :
    StoichSubstanceSSTP()
{
    XML_Node* root = get_XML_File(infile);
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("MineralEQ3::MineralEQ3",
                           "Couldn't find phase name in file:" + id_);
    }
    // Check the model name to ensure we have compatibility
    const XML_Node& th = xphase->child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "MineralEQ3") {
        throw CanteraError("MineralEQ3::MineralEQ3",
                           "thermo model attribute must be StoichSubstance");
    }
    importPhase(*xphase, this);
}

MineralEQ3::MineralEQ3(XML_Node& xmlphase, const std::string& id_) :
    StoichSubstanceSSTP()
{
    if (id_ != "") {
        std::string idxml = xmlphase["id"];
        if (id_ != idxml) {
            throw CanteraError("MineralEQ3::MineralEQ3",
                               "id's don't match");
        }
    }
    const XML_Node& th = xmlphase.child("thermo");
    std::string model = th["model"];
    if (model != "StoichSubstance" && model != "MineralEQ3") {
        throw CanteraError("MineralEQ3::MineralEQ3",
                           "thermo model attribute must be StoichSubstance");
    }
    importPhase(xmlphase, this);
}

MineralEQ3::MineralEQ3(const MineralEQ3&  right) :
    StoichSubstanceSSTP()
{
    *this = operator=(right);
}

MineralEQ3&
MineralEQ3::operator=(const MineralEQ3& right)
{
    if (&right == this) {
        return *this;
    }
    StoichSubstanceSSTP::operator=(right);
    m_Mu0_pr_tr = right.m_Mu0_pr_tr;
    m_Entrop_pr_tr = right.m_Entrop_pr_tr;
    m_deltaG_formation_pr_tr = right.m_deltaG_formation_pr_tr;
    m_deltaH_formation_pr_tr = right.m_deltaH_formation_pr_tr;
    m_V0_pr_tr               = right.m_V0_pr_tr;
    m_a                      = right.m_a;
    m_b                      = right.m_b;
    m_c                      = right.m_c;

    return *this;
}

ThermoPhase* MineralEQ3::duplMyselfAsThermoPhase() const
{
    return new MineralEQ3(*this);
}

/*
 * ---- Utilities -----
 */

int MineralEQ3::eosType() const
{
    return cStoichSubstance;
}

/*
 * ----- Mechanical Equation of State ------
 */

doublereal MineralEQ3::pressure() const
{
    return m_press;
}

void MineralEQ3::setPressure(doublereal p)
{
    m_press = p;
}

doublereal MineralEQ3::isothermalCompressibility() const
{
    return 0.0;
}

doublereal MineralEQ3::thermalExpansionCoeff() const
{
    return 0.0;
}

/*
 * ---- Chemical Potentials and Activities ----
 */

void MineralEQ3::getActivityConcentrations(doublereal* c) const
{
    c[0] = 1.0;
}

doublereal MineralEQ3::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal MineralEQ3::logStandardConc(size_t k) const
{
    return 0.0;
}

void MineralEQ3::getUnitsStandardConc(doublereal* uA, int k, int sizeUA) const
{
    for (int i = 0; i < 6; i++) {
        uA[i] = 0;
    }
}

/*
 * Properties of the Standard State of the Species in the Solution
 */

void MineralEQ3::getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= GasConstant * temperature();
}

void MineralEQ3::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
    doublereal RT = GasConstant * temperature();
    doublereal presCorrect = (m_press - m_p0) /  molarDensity();
    hrt[0] += presCorrect / RT;
}

void MineralEQ3::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
}

void MineralEQ3::getGibbs_RT(doublereal* grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R[0];
}

void MineralEQ3::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R[0];
}

void MineralEQ3::getIntEnergy_RT(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    urt[0] = m_h0_RT[0] - m_p0 / molarDensity() / RT;
}

/*
 * ---- Thermodynamic Values for the Species Reference States ----
 */

void MineralEQ3::getIntEnergy_RT_ref(doublereal* urt) const
{
    _updateThermo();
    doublereal RT = GasConstant * temperature();
    urt[0] = m_h0_RT[0] - m_p0 / molarDensity() / RT;
}

/*
 * ---- Initialization and Internal functions
 */

void MineralEQ3::initThermo()
{
    /*
     * Call the base class thermo initializer
     */
    StoichSubstanceSSTP::initThermo();
}

void MineralEQ3::setParameters(int n, doublereal* const c)
{
    doublereal rho = c[0];
    setDensity(rho);
}

void MineralEQ3::getParameters(int& n, doublereal* const c) const
{
    doublereal rho = density();
    n = 1;
    c[0] = rho;
}

void MineralEQ3::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("HMWSoln::initThermoXML",
                           "no thermo XML node");
    }

    std::vector<const XML_Node*> xspecies = speciesData();
    const XML_Node* xsp = xspecies[0];

    XML_Node* aStandardState = 0;
    if (xsp->hasChild("standardState")) {
        aStandardState = &xsp->child("standardState");
    } else {
        throw CanteraError("MineralEQ3::initThermoXML",
                           "no standard state mode");
    }
    doublereal volVal = 0.0;
    string smodel = (*aStandardState)["model"];
    if (smodel != "constantVolume") {
        throw CanteraError("MineralEQ3::initThermoXML",
                           "wrong standard state mode");
    }
    if (aStandardState->hasChild("V0_Pr_Tr")) {
        XML_Node& aV = aStandardState->child("V0_Pr_Tr");
        string Aunits = "";
        double Afactor = toSI("cm3/gmol");
        if (aV.hasAttrib("units")) {
            Aunits = aV.attrib("units");
            Afactor = toSI(Aunits);
        }
        volVal = ctml::getFloat(*aStandardState, "V0_Pr_Tr");
        m_V0_pr_tr= volVal;
        volVal *= Afactor;
        m_speciesSize[0] = volVal;
    } else {
        throw CanteraError("MineralEQ3::initThermoXML",
                           "wrong standard state mode");
    }
    doublereal rho = molecularWeight(0) / volVal;
    setDensity(rho);

    const XML_Node& sThermo = xsp->child("thermo");
    const XML_Node& MinEQ3node = sThermo.child("MinEQ3");


    m_deltaG_formation_pr_tr =
        ctml::getFloatDefaultUnits(MinEQ3node, "DG0_f_Pr_Tr", "cal/gmol", "actEnergy");
    m_deltaH_formation_pr_tr =
        ctml::getFloatDefaultUnits(MinEQ3node, "DH0_f_Pr_Tr", "cal/gmol", "actEnergy");
    m_Entrop_pr_tr = ctml::getFloatDefaultUnits(MinEQ3node, "S0_Pr_Tr", "cal/gmol/K");
    m_a = ctml::getFloatDefaultUnits(MinEQ3node, "a", "cal/gmol/K");
    m_b = ctml::getFloatDefaultUnits(MinEQ3node, "b", "cal/gmol/K2");
    m_c = ctml::getFloatDefaultUnits(MinEQ3node, "c", "cal-K/gmol");

    convertDGFormation();
}

void MineralEQ3::setParametersFromXML(const XML_Node& eosdata)
{
    std::string model = eosdata["model"];
    if (model != "MineralEQ3") {
        throw CanteraError("MineralEQ3::MineralEQ3",
                           "thermo model attribute must be MineralEQ3");
    }
}

doublereal MineralEQ3::LookupGe(const std::string& elemName)
{
    size_t iE = elementIndex(elemName);
    if (iE == npos) {
        throw CanteraError("PDSS_HKFT::LookupGe", "element " + elemName + " not found");
    }
    doublereal geValue = entropyElement298(iE);
    if (geValue == ENTROPY298_UNKNOWN) {
        throw CanteraError("PDSS_HKFT::LookupGe",
                           "element " + elemName + " does not have a supplied entropy298");
    }
    geValue *= (-298.15);
    return geValue;
}

void MineralEQ3::convertDGFormation()
{
    /*
     * Ok let's get the element compositions and conversion factors.
     */
    doublereal na;
    doublereal ge;
    string ename;

    doublereal totalSum = 0.0;
    for (size_t m = 0; m < nElements(); m++) {
        na = nAtoms(0, m);
        if (na > 0.0) {
            ename = elementName(m);
            ge = LookupGe(ename);
            totalSum += na * ge;
        }
    }
    // Ok, now do the calculation. Convert to joules kmol-1
    doublereal dg = m_deltaG_formation_pr_tr * 4.184 * 1.0E3;
    //! Store the result into an internal variable.
    m_Mu0_pr_tr = dg + totalSum;
}

}
