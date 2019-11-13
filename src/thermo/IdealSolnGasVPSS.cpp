/**
 *  @file IdealSolnGasVPSS.cpp
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::IdealSolnGasVPSS IdealSolnGasVPSS\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

IdealSolnGasVPSS::IdealSolnGasVPSS() :
    m_idealGas(-1),
    m_formGC(0)
{
}

IdealSolnGasVPSS::IdealSolnGasVPSS(const std::string& infile, std::string id_) :
    m_idealGas(0),
    m_formGC(0)
{
    XML_Node* root = get_XML_File(infile);
    if (id_ == "-") {
        id_ = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id_, root);
    if (!xphase) {
        throw CanteraError("IdealSolnGasVPSS::IdealSolnGasVPSS",
                           "Couldn't find phase named '{} in file '{}'",
                           id_, infile);
    }
    importPhase(*xphase, this);
}

void IdealSolnGasVPSS::setStandardConcentrationModel(const std::string& model)
{
    if (m_idealGas) {
        throw CanteraError("IdealSolnGasVPSS::setStandardConcentrationModel",
                           "Standard concentration model not applicable for ideal gas");
    }

    if (caseInsensitiveEquals(model, "unity")) {
        m_formGC = 0;
    } else if (caseInsensitiveEquals(model, "species-molar-volume")
               || caseInsensitiveEquals(model, "molar_volume")) {
        m_formGC = 1;
    } else if (caseInsensitiveEquals(model, "solvent-molar-volume")
               || caseInsensitiveEquals(model, "solvent_volume")) {
        m_formGC = 2;
    } else {
        throw CanteraError("IdealSolnGasVPSS::setStandardConcentrationModel",
                           "Unknown standard concentration model '{}'", model);
    }
}

// ------------Molar Thermodynamic Properties -------------------------

doublereal IdealSolnGasVPSS::enthalpy_mole() const
{
    updateStandardStateThermo();
    return RT() * mean_X(m_hss_RT);
}

doublereal IdealSolnGasVPSS::entropy_mole() const
{
    updateStandardStateThermo();
    return GasConstant * (mean_X(m_sss_R) - sum_xlogx());
}

doublereal IdealSolnGasVPSS::cp_mole() const
{
    updateStandardStateThermo();
    return GasConstant * mean_X(m_cpss_R);
}

doublereal IdealSolnGasVPSS::cv_mole() const
{
    return cp_mole() - GasConstant;
}

void IdealSolnGasVPSS::setPressure(doublereal p)
{
    m_Pcurrent = p;
    updateStandardStateThermo();
    calcDensity();
}

void IdealSolnGasVPSS::calcDensity()
{
    // Calculate the molarVolume of the solution (m**3 kmol-1)
    if (m_idealGas) {
        double dens = m_Pcurrent * meanMolecularWeight() / RT();
        Phase::assignDensity(dens);
    } else {
        const doublereal* const dtmp = moleFractdivMMW();
        const vector_fp& vss = getStandardVolumes();
        double dens = 1.0 / dot(vss.begin(), vss.end(), dtmp);

        // Set the density in the parent State object directly
        Phase::assignDensity(dens);
    }
}

doublereal IdealSolnGasVPSS::isothermalCompressibility() const
{
    if (m_idealGas) {
        return 1.0 / m_Pcurrent;
    } else {
        throw NotImplementedError("IdealSolnGasVPSS::isothermalCompressibility");
    }
}

Units IdealSolnGasVPSS::standardConcentrationUnits() const
{
    if (m_idealGas || m_formGC != 0) {
        return Units(1.0, 0, -3, 0, 0, 0, 1);
    } else {
        return Units(1.0);
    }
}

void IdealSolnGasVPSS::getActivityConcentrations(doublereal* c) const
{
    if (m_idealGas) {
        getConcentrations(c);
    } else {
        const vector_fp& vss = getStandardVolumes();
        switch (m_formGC) {
        case 0:
            for (size_t k = 0; k < m_kk; k++) {
                c[k] = moleFraction(k);
            }
            break;
        case 1:
            for (size_t k = 0; k < m_kk; k++) {
                c[k] = moleFraction(k) / vss[k];
            }
            break;
        case 2:
            for (size_t k = 0; k < m_kk; k++) {
                c[k] = moleFraction(k) / vss[0];
            }
            break;
        }
    }
}

doublereal IdealSolnGasVPSS::standardConcentration(size_t k) const
{
    if (m_idealGas) {
        return pressure() / RT();
    } else {
        const vector_fp& vss = getStandardVolumes();
        switch (m_formGC) {
        case 0:
            return 1.0;
        case 1:
            return 1.0 / vss[k];
        case 2:
            return 1.0/ vss[0];
        }
        return 0.0;
    }
}

void IdealSolnGasVPSS::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

// ---- Partial Molar Properties of the Solution -----------------

void IdealSolnGasVPSS::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += RT() * log(xx);
    }
}

void IdealSolnGasVPSS::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    scale(hbar, hbar+m_kk, hbar, RT());
}

void IdealSolnGasVPSS::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    scale(sbar, sbar+m_kk, sbar, GasConstant);
    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += GasConstant * (- log(xx));
    }
}

void IdealSolnGasVPSS::getPartialMolarIntEnergies(doublereal* ubar) const
{
    getIntEnergy_RT(ubar);
    scale(ubar, ubar+m_kk, ubar, RT());
}

void IdealSolnGasVPSS::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    scale(cpbar, cpbar+m_kk, cpbar, GasConstant);
}

void IdealSolnGasVPSS::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

void IdealSolnGasVPSS::setToEquilState(const doublereal* mu_RT)
{
    updateStandardStateThermo();

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    doublereal pres = 0.0;
    double m_p0 = refPressure();
    for (size_t k = 0; k < m_kk; k++) {
        double tmp = -m_g0_RT[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            double tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setState_PX(pres, &m_pp[0]);
}

bool IdealSolnGasVPSS::addSpecies(shared_ptr<Species> spec)
{
    bool added = VPStandardStateTP::addSpecies(spec);
    if (added) {
        m_pp.push_back(0.0);
    }
    return added;
}

void IdealSolnGasVPSS::initThermo()
{
    VPStandardStateTP::initThermo();
    if (m_input.hasKey("thermo")) {
        const string& model = m_input["thermo"].asString();
        if (model == "ideal-solution-VPSS") {
            setSolnMode();
        } else if (model == "ideal-gas-VPSS") {
            setGasMode();
        }
    }
    if (m_input.hasKey("standard-concentration-basis")) {
        setStandardConcentrationModel(m_input["standard-concentration-basis"].asString());
    }

    if (m_idealGas == -1) {
        throw CanteraError("IdealSolnGasVPSS::initThermo",
            "solution / gas mode not set");
    }
}

void IdealSolnGasVPSS::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thermoNode = phaseNode.child("thermo");
        std::string model = thermoNode["model"];
        if (model == "IdealGasVPSS") {
            setGasMode();
        } else if (model == "IdealSolnVPSS") {
            setSolnMode();
        } else {
            throw CanteraError("IdealSolnGasVPSS::initThermoXML",
                               "Unknown thermo model : " + model);
        }
    }

    // Form of the standard concentrations. Must have one of:
    //
    //     <standardConc model="unity" />
    //     <standardConc model="molar_volume" />
    //     <standardConc model="solvent_volume" />
    if (phaseNode.hasChild("standardConc")) {
        XML_Node& scNode = phaseNode.child("standardConc");
        setStandardConcentrationModel(scNode.attrib("model"));
    } else if (!m_idealGas) {
        throw CanteraError("IdealSolnGasVPSS::initThermoXML",
                           "Unspecified standardConc model");
    }

    VPStandardStateTP::initThermoXML(phaseNode, id_);
}

void IdealSolnGasVPSS::setParametersFromXML(const XML_Node& thermoNode)
{
    VPStandardStateTP::setParametersFromXML(thermoNode);
    std::string model = thermoNode["model"];
    if (model == "IdealGasVPSS") {
        setGasMode();
    } else if (model == "IdealSolnVPSS") {
        setSolnMode();
    } else {
        throw CanteraError("IdealSolnGasVPSS::setParametersFromXML",
                           "Unknown thermo model : " + model);
    }
}

}
