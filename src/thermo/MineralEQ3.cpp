/**
 *  @file MineralEQ3.cpp
 *   Definition file for the MineralEQ3 class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::MineralEQ3 MineralEQ3\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/thermo/MineralEQ3.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

// ----  Constructors -------

MineralEQ3::MineralEQ3(const std::string& infile, const std::string& id_)
{
    warn_deprecated("Class MineralEQ3", "To be removed after Cantera 2.4");
    initThermoFile(infile, id_);
}

MineralEQ3::MineralEQ3(XML_Node& xmlphase, const std::string& id_)
{
    warn_deprecated("Class MineralEQ3", "To be removed after Cantera 2.4");
    importPhase(xmlphase, this);
}

// ----- Mechanical Equation of State ------

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

// ---- Chemical Potentials and Activities ----

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

// Properties of the Standard State of the Species in the Solution

void MineralEQ3::getStandardChemPotentials(doublereal* mu0) const
{
    getGibbs_RT(mu0);
    mu0[0] *= RT();
}

void MineralEQ3::getEnthalpy_RT(doublereal* hrt) const
{
    getEnthalpy_RT_ref(hrt);
    doublereal presCorrect = (m_press - m_p0) / molarDensity();
    hrt[0] += presCorrect / RT();
}

void MineralEQ3::getEntropy_R(doublereal* sr) const
{
    getEntropy_R_ref(sr);
}

void MineralEQ3::getGibbs_RT(doublereal* grt) const
{
    getEnthalpy_RT(grt);
    grt[0] -= m_s0_R;
}

void MineralEQ3::getCp_R(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R;
}

void MineralEQ3::getIntEnergy_RT(doublereal* urt) const
{
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Thermodynamic Values for the Species Reference States ----

void MineralEQ3::getIntEnergy_RT_ref(doublereal* urt) const
{
    _updateThermo();
    urt[0] = m_h0_RT - m_p0 / molarDensity() / RT();
}

// ---- Initialization and Internal functions

void MineralEQ3::setParameters(int n, doublereal* const c)
{
    setDensity(c[0]);
}

void MineralEQ3::getParameters(int& n, doublereal* const c) const
{
    n = 1;
    c[0] = density();
}

void MineralEQ3::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    // Find the Thermo XML node
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("HMWSoln::initThermoXML",
                           "no thermo XML node");
    }

    const XML_Node* xsp = speciesData()[0];

    XML_Node* aStandardState = 0;
    if (xsp->hasChild("standardState")) {
        aStandardState = &xsp->child("standardState");
    } else {
        throw CanteraError("MineralEQ3::initThermoXML",
                           "no standard state mode");
    }
    doublereal volVal = 0.0;
    if (aStandardState->attrib("model") != "constantVolume") {
        throw CanteraError("MineralEQ3::initThermoXML",
                           "wrong standard state mode");
    }
    if (aStandardState->hasChild("V0_Pr_Tr")) {
        XML_Node& aV = aStandardState->child("V0_Pr_Tr");
        double Afactor = toSI("cm3/gmol");
        if (aV.hasAttrib("units")) {
            Afactor = toSI(aV.attrib("units"));
        }
        volVal = getFloat(*aStandardState, "V0_Pr_Tr");
        m_V0_pr_tr= volVal;
        volVal *= Afactor;
    } else {
        throw CanteraError("MineralEQ3::initThermoXML",
                           "wrong standard state mode");
    }
    setDensity(molecularWeight(0) / volVal);

    const XML_Node& MinEQ3node = xsp->child("thermo").child("MinEQ3");

    m_deltaG_formation_pr_tr =
        getFloat(MinEQ3node, "DG0_f_Pr_Tr", "actEnergy") / actEnergyToSI("cal/gmol");
    m_deltaH_formation_pr_tr =
        getFloat(MinEQ3node, "DH0_f_Pr_Tr", "actEnergy") / actEnergyToSI("cal/gmol");
    m_Entrop_pr_tr = getFloat(MinEQ3node, "S0_Pr_Tr", "toSI") / toSI("cal/gmol/K");
    m_a = getFloat(MinEQ3node, "a", "toSI") / toSI("cal/gmol/K");
    m_b = getFloat(MinEQ3node, "b", "toSI") / toSI("cal/gmol/K2");
    m_c = getFloat(MinEQ3node, "c", "toSI") / toSI("cal-K/gmol");

    convertDGFormation();
}

void MineralEQ3::setParametersFromXML(const XML_Node& eosdata)
{
    if (eosdata["model"] != "MineralEQ3") {
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
    // Ok let's get the element compositions and conversion factors.
    doublereal totalSum = 0.0;
    for (size_t m = 0; m < nElements(); m++) {
        double na = nAtoms(0, m);
        if (na > 0.0) {
            totalSum += na * LookupGe(elementName(m));
        }
    }
    // Ok, now do the calculation. Convert to joules kmol-1
    doublereal dg = m_deltaG_formation_pr_tr * toSI("cal/gmol");
    //! Store the result into an internal variable.
    m_Mu0_pr_tr = dg + totalSum;

    double Hcalc = m_Mu0_pr_tr + 298.15 * m_Entrop_pr_tr * toSI("cal/gmol");
    double DHjmol = m_deltaH_formation_pr_tr * toSI("kal/gmol");

    // If the discrepancy is greater than 100 cal gmol-1, print an error
    if (fabs(Hcalc -DHjmol) > 100 * toSI("cal/gmol")) {
        throw CanteraError("installMinEQ3asShomateThermoFromXML()",
                           "DHjmol is not consistent with G and S: {} vs {}",
                           Hcalc, DHjmol);
    }
}

}
