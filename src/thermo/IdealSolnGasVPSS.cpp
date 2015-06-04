/**
 *  @file IdealSolnGasVPSS.cpp
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::IdealSolnGasVPSS IdealSolnGasVPSS\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/IdealSolnGasVPSS.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/vec_functions.h"

using namespace std;

namespace Cantera
{

IdealSolnGasVPSS::IdealSolnGasVPSS() :
    m_idealGas(0),
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
        throw CanteraError("newPhase",
                           "Couldn't find phase named \"" + id_ + "\" in file, " + infile);
    }
    importPhase(*xphase, this);
}

IdealSolnGasVPSS::IdealSolnGasVPSS(const IdealSolnGasVPSS& b) :
    m_idealGas(0),
    m_formGC(0)
{
    *this = b;
}

IdealSolnGasVPSS& IdealSolnGasVPSS::operator=(const IdealSolnGasVPSS& b)
{
    if (&b != this) {
        /*
         * Mostly, this is a passthrough to the underlying
         * assignment operator for the ThermoPhae parent object.
         */
        VPStandardStateTP::operator=(b);
        /*
         * However, we have to handle data that we own.
         */
        m_idealGas = b.m_idealGas;
        m_formGC   = b.m_formGC;
    }
    return *this;
}

ThermoPhase* IdealSolnGasVPSS::duplMyselfAsThermoPhase() const
{
    return new IdealSolnGasVPSS(*this);
}

int IdealSolnGasVPSS::eosType() const
{
    if (m_idealGas) {
        return cIdealSolnGasVPSS;
    }
    return cIdealSolnGasVPSS_iscv;
}

/*
 * ------------Molar Thermodynamic Properties -------------------------
 */

doublereal IdealSolnGasVPSS::enthalpy_mole() const
{
    updateStandardStateThermo();
    return GasConstant * temperature() * mean_X(m_VPSS_ptr->enthalpy_RT());
}

doublereal IdealSolnGasVPSS::entropy_mole() const
{
    updateStandardStateThermo();
    return GasConstant * (mean_X(m_VPSS_ptr->entropy_R()) - sum_xlogx());

}

doublereal IdealSolnGasVPSS::cp_mole() const
{
    updateStandardStateThermo();
    return  GasConstant * mean_X(m_VPSS_ptr->cp_R());
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
    /*
     * Calculate the molarVolume of the solution (m**3 kmol-1)
     */
    if (m_idealGas) {
        double dens = (m_Pcurrent * meanMolecularWeight()
                       /(GasConstant * temperature()));
        Phase::setDensity(dens);
    } else {
        const doublereal* const dtmp = moleFractdivMMW();
        const vector_fp& vss = m_VPSS_ptr->getStandardVolumes();
        double dens = 1.0 / dot(vss.begin(), vss.end(), dtmp);
        /*
         * Set the density in the parent State object directly,
         * by calling the Phase::setDensity() function.
         */
        Phase::setDensity(dens);
    }
}

doublereal IdealSolnGasVPSS::isothermalCompressibility() const
{
    if (m_idealGas) {
        return 1.0 / m_Pcurrent;
    } else {
        throw CanteraError("IdealSolnGasVPSS::isothermalCompressibility() ",
                           "not implemented");
    }
    return 0.0;
}

void IdealSolnGasVPSS::getActivityConcentrations(doublereal* c) const
{
    if (m_idealGas) {
        getConcentrations(c);
    } else {
        const vector_fp& vss = m_VPSS_ptr->getStandardVolumes();
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
        double p = pressure();
        return p/(GasConstant * temperature());
    } else {
        const vector_fp& vss = m_VPSS_ptr->getStandardVolumes();
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

void IdealSolnGasVPSS::getUnitsStandardConc(double* uA, int, int sizeUA) const
{
    warn_deprecated("IdealSolnGasVPSS::getUnitsStandardConc",
                "To be removed after Cantera 2.2.");

    if (eosType() == cIdealSolnGasPhase0) {
        for (int i = 0; i < sizeUA; i++) {
            uA[i] = 0.0;
        }
    } else {
        for (int i = 0; i < sizeUA; i++) {
            if (i == 0) {
                uA[0] = 1.0;
            }
            if (i == 1) {
                uA[1] = -int(nDim());
            }
            if (i == 2) {
                uA[2] = 0.0;
            }
            if (i == 3) {
                uA[3] = 0.0;
            }
            if (i == 4) {
                uA[4] = 0.0;
            }
            if (i == 5) {
                uA[5] = 0.0;
            }
        }
    }
}

void IdealSolnGasVPSS::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

/*
 * ---- Partial Molar Properties of the Solution -----------------
 */

void IdealSolnGasVPSS::getChemPotentials_RT(doublereal* muRT) const
{
    getChemPotentials(muRT);
    doublereal invRT = 1.0 / _RT();
    for (size_t k = 0; k < m_kk; k++) {
        muRT[k] *= invRT;
    }
}

void IdealSolnGasVPSS::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
    doublereal rt = temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += rt*(log(xx));
    }
}

void IdealSolnGasVPSS::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    scale(hbar, hbar+m_kk, hbar, GasConstant * temperature());
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
    scale(ubar, ubar+m_kk, ubar, GasConstant * temperature());
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

void IdealSolnGasVPSS::initThermo()
{
    initLengths();
    VPStandardStateTP::initThermo();
}

void IdealSolnGasVPSS::setToEquilState(const doublereal* mu_RT)
{
    updateStandardStateThermo();
    const vector_fp& grt = m_VPSS_ptr->Gibbs_RT_ref();

    /*
     * Within the method, we protect against inf results if the
     * exponent is too high.
     *
     * If it is too low, we set
     * the partial pressure to zero. This capability is needed
     * by the elemental potential method.
     */
    doublereal pres = 0.0;
    double m_p0 = m_VPSS_ptr->refPressure();
    for (size_t k = 0; k < m_kk; k++) {
        double tmp = -grt[k] + mu_RT[k];
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

void IdealSolnGasVPSS::initLengths()
{
    m_pp.resize(m_kk, 0.0);
}

void IdealSolnGasVPSS::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    IdealSolnGasVPSS::initLengths();

    if (phaseNode.hasChild("thermo")) {
        XML_Node& thermoNode = phaseNode.child("thermo");
        std::string model = thermoNode["model"];
        if (model == "IdealGasVPSS") {
            m_idealGas = 1;
        } else if (model == "IdealSolnVPSS") {
            m_idealGas = 0;
        } else {
            throw CanteraError("IdealSolnGasVPSS::initThermoXML",
                               "Unknown thermo model : " + model);
        }
    }

    /*
     * Form of the standard concentrations. Must have one of:
     *
     *     <standardConc model="unity" />
     *     <standardConc model="molar_volume" />
     *     <standardConc model="solvent_volume" />
     */
    if (phaseNode.hasChild("standardConc")) {
        if (m_idealGas) {
            throw CanteraError("IdealSolnGasVPSS::initThermoXML",
                               "standardConc node for ideal gas");
        }
        XML_Node& scNode = phaseNode.child("standardConc");
        string formStringa = scNode.attrib("model");
        string formString = lowercase(formStringa);
        if (formString == "unity") {
            m_formGC = 0;
        } else if (formString == "molar_volume") {
            m_formGC = 1;
        } else if (formString == "solvent_volume") {
            m_formGC = 2;
        } else {
            throw CanteraError("initThermoXML",
                               "Unknown standardConc model: " + formStringa);
        }
    } else {
        if (!m_idealGas) {
            throw CanteraError("initThermoXML",
                               "Unspecified standardConc model");
        }
    }

    VPStandardStateTP::initThermoXML(phaseNode, id_);
}

void IdealSolnGasVPSS::setParametersFromXML(const XML_Node& thermoNode)
{
    VPStandardStateTP::setParametersFromXML(thermoNode);
    std::string model = thermoNode["model"];
    if (model == "IdealGasVPSS") {
        m_idealGas = 1;
    } else if (model == "IdealSolnVPSS") {
        m_idealGas = 0;
    } else {
        throw CanteraError("IdealSolnGasVPSS::initThermoXML",
                           "Unknown thermo model : " + model);
    }
}

}
