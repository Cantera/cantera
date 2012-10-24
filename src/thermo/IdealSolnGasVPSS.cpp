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
#include "cantera/thermo/VPSSMgr.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

/*
 * Default constructor
 */
IdealSolnGasVPSS::IdealSolnGasVPSS() :
    VPStandardStateTP(),
    m_idealGas(0),
    m_formGC(0)
{
}


IdealSolnGasVPSS::IdealSolnGasVPSS(const std::string& infile, std::string id) :
    VPStandardStateTP(),
    m_idealGas(0),
    m_formGC(0)
{
    XML_Node* root = get_XML_File(infile);
    if (id == "-") {
        id = "";
    }
    XML_Node* xphase = get_XML_NameID("phase", std::string("#")+id, root);
    if (!xphase) {
        throw CanteraError("newPhase",
                           "Couldn't find phase named \"" + id + "\" in file, " + infile);
    }
    importPhase(*xphase, this);
}

/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor.
 *
 *  The copy constructor just calls the assignment operator
 *  to do the heavy lifting.
 */
IdealSolnGasVPSS::IdealSolnGasVPSS(const IdealSolnGasVPSS& b) :
    VPStandardStateTP(),
    m_idealGas(0),
    m_formGC(0)
{
    *this = b;
}

/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
IdealSolnGasVPSS& IdealSolnGasVPSS::
operator=(const IdealSolnGasVPSS& b)
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

/*
 * ~IdealSolnGasVPSS():   (virtual)
 *
 */
IdealSolnGasVPSS::~IdealSolnGasVPSS()
{
}

/*
 * Duplication function.
 *  This calls the copy constructor for this object.
 */
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

/// Molar enthalpy. Units: J/kmol.
doublereal IdealSolnGasVPSS::enthalpy_mole() const
{
    updateStandardStateThermo();
    const vector_fp& enth_RT = m_VPSS_ptr->enthalpy_RT();
    return (GasConstant * temperature() *
            mean_X(DATA_PTR(enth_RT)));
}

/// Molar internal energy. Units: J/kmol.
doublereal IdealSolnGasVPSS::intEnergy_mole() const
{
    doublereal p0 = pressure();
    doublereal md = molarDensity();
    return (enthalpy_mole() - p0 / md);
}

/// Molar entropy. Units: J/kmol/K.
doublereal IdealSolnGasVPSS::entropy_mole() const
{
    updateStandardStateThermo();
    const vector_fp& entrop_R = m_VPSS_ptr->entropy_R();
    return GasConstant * (mean_X(DATA_PTR(entrop_R)) - sum_xlogx());

}

/// Molar Gibbs function. Units: J/kmol.
doublereal IdealSolnGasVPSS::gibbs_mole() const
{
    return enthalpy_mole() - temperature() * entropy_mole();
}

/// Molar heat capacity at constant pressure. Units: J/kmol/K.
doublereal IdealSolnGasVPSS::cp_mole() const
{
    updateStandardStateThermo();
    const vector_fp& cp_R = m_VPSS_ptr->cp_R();
    return  GasConstant * (mean_X(DATA_PTR(cp_R)));
}

/// Molar heat capacity at constant volume. Units: J/kmol/K.
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
        const vector_fp& vss = m_VPSS_ptr->standardVolumes();
        double invDens = dot(vss.begin(), vss.end(), dtmp);
        /*
         * Set the density in the parent State object directly,
         * by calling the Phase::setDensity() function.
         */
        double dens = 1.0/invDens;
        Phase::setDensity(dens);
    }
}

doublereal IdealSolnGasVPSS::isothermalCompressibility() const
{
    if (m_idealGas) {
        return -1.0 / m_Pcurrent;
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
        const vector_fp& vss = m_VPSS_ptr->standardVolumes();
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

/*
 * Returns the standard concentration \f$ C^0_k \f$, which is used to normalize
 * the generalized concentration.
 */
doublereal IdealSolnGasVPSS::standardConcentration(size_t k) const
{
    if (m_idealGas) {
        double p = pressure();
        return p/(GasConstant * temperature());
    } else {
        const vector_fp& vss = m_VPSS_ptr->standardVolumes();
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

/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal IdealSolnGasVPSS::logStandardConc(size_t k) const
{
    double c = standardConcentration(k);
    double lc = std::log(c);
    return lc;
}

/*
 *
 * getUnitsStandardConcentration()
 *
 * Returns the units of the standard and general concentrations
 * Note they have the same units, as their divisor is
 * defined to be equal to the activity of the kth species
 * in the solution, which is unitless.
 *
 * This routine is used in print out applications where the
 * units are needed. Usually, MKS units are assumed throughout
 * the program and in the XML input files.
 *
 *  uA[0] = kmol units - default  = 1
 *  uA[1] = m    units - default  = -nDim(), the number of spatial
 *                                dimensions in the Phase class.
 *  uA[2] = kg   units - default  = 0;
 *  uA[3] = Pa(pressure) units - default = 0;
 *  uA[4] = Temperature units - default = 0;
 *  uA[5] = time units - default = 0
 *
 *  For EOS types other than cIdealSolidSolnPhase1, the default
 *  kmol/m3 holds for standard concentration units. For
 *  cIdealSolidSolnPhase0 type, the standard concentration is
 *  unitless.
 */
void IdealSolnGasVPSS::getUnitsStandardConc(double* uA, int, int sizeUA) const
{
    int eos = eosType();
    if (eos == cIdealSolnGasPhase0) {
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


/*
 * Get the array of non-dimensional activity coefficients
 */
void IdealSolnGasVPSS::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

/*
 * ---- Partial Molar Properties of the Solution -----------------
 */

/*
 * Get the array of non-dimensional species chemical potentials
 * These are partial molar Gibbs free energies.
 * \f$ \mu_k / \hat R T \f$.
 * Units: unitless
 *
 * We close the loop on this function, here, calling
 * getChemPotentials() and then dividing by RT.
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
    doublereal xx;
    doublereal rt = temperature() * GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        mu[k] += rt*(log(xx));
    }
}


void IdealSolnGasVPSS::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    doublereal rt = GasConstant * temperature();
    scale(hbar, hbar+m_kk, hbar, rt);
}

void IdealSolnGasVPSS::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    doublereal r = GasConstant;
    scale(sbar, sbar+m_kk, sbar, r);
    for (size_t k = 0; k < m_kk; k++) {
        doublereal xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] += r * (- log(xx));
    }
}

void IdealSolnGasVPSS::getPartialMolarIntEnergies(doublereal* ubar) const
{
    getIntEnergy_RT(ubar);
    doublereal rt = GasConstant * temperature();
    scale(ubar, ubar+m_kk, ubar, rt);
}

void IdealSolnGasVPSS::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    doublereal r = GasConstant;
    scale(cpbar, cpbar+m_kk, cpbar, r);
}

void IdealSolnGasVPSS::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

/*
 * ----- Thermodynamic Values for the Species Reference States ----
 */




/*
 * Perform initializations after all species have been
 * added.
 */
void IdealSolnGasVPSS::initThermo()
{
    initLengths();
    VPStandardStateTP::initThermo();
}


void IdealSolnGasVPSS::setToEquilState(const doublereal* mu_RT)
{
    double tmp, tmp2;
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
        tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            tmp2 = tmp / 500.;
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

/*
 * Initialize the internal lengths.
 *       (this is not a virtual function)
 */
void IdealSolnGasVPSS::initLengths()
{
    m_kk = nSpecies();
    m_pp.resize(m_kk, 0.0);
}

/*
 *   Import and initialize a ThermoPhase object
 *
 * param phaseNode This object must be the phase node of a
 *             complete XML tree
 *             description of the phase, including all of the
 *             species data. In other words while "phase" must
 *             point to an XML phase object, it must have
 *             sibling nodes "speciesData" that describe
 *             the species in the phase.
 * param id   ID of the phase. If nonnull, a check is done
 *             to see if phaseNode is pointing to the phase
 *             with the correct id.
 *
 * This routine initializes the lengths in the current object and
 * then calls the parent routine.
 */
void IdealSolnGasVPSS::initThermoXML(XML_Node& phaseNode, const std::string& id)
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

    VPStandardStateTP::initThermoXML(phaseNode, id);
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


