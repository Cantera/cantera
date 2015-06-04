/**
 *  @file IdealSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see \ref
 *      thermoprops and \link Cantera::IdealSolidSolnPhase
 *      IdealSolidSolnPhase\endlink).
 */
/*
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 */

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/vec_functions.h"

using namespace std;

namespace Cantera
{

IdealSolidSolnPhase::IdealSolidSolnPhase(int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" IdealSolidSolnPhase Constructor",
                           " Illegal value of formGC");
    }
}

IdealSolidSolnPhase::IdealSolidSolnPhase(const std::string& inputFile,
        const std::string& id_, int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" IdealSolidSolnPhase Constructor",
                           " Illegal value of formGC");
    }
    initThermoFile(inputFile, id_);
}

IdealSolidSolnPhase::IdealSolidSolnPhase(XML_Node& root, const std::string& id_,
        int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" IdealSolidSolnPhase Constructor",
                           " Illegal value of formGC");
    }
    importPhase(*findXMLPhase(&root, id_), this);
}

IdealSolidSolnPhase::IdealSolidSolnPhase(const IdealSolidSolnPhase& b)
{
    *this = b;
}

IdealSolidSolnPhase& IdealSolidSolnPhase::operator=(const IdealSolidSolnPhase& b)
{
    if (this != &b) {
        ThermoPhase::operator=(b);

        m_formGC     = b.m_formGC;
        m_Pref       = b.m_Pref;
        m_Pcurrent   = b.m_Pcurrent;
        m_speciesMolarVolume = b.m_speciesMolarVolume;
        m_h0_RT      = b.m_h0_RT;
        m_cp0_R      = b.m_cp0_R;
        m_g0_RT      = b.m_g0_RT;
        m_s0_R       = b.m_s0_R;
        m_expg0_RT   = b.m_expg0_RT;
        m_pe         = b.m_pe;
        m_pp         = b.m_pp;
    }
    return *this;
}

ThermoPhase* IdealSolidSolnPhase::duplMyselfAsThermoPhase() const
{
    return new IdealSolidSolnPhase(*this);
}

int IdealSolidSolnPhase::eosType() const
{
    integer res;
    switch (m_formGC) {
    case 0:
        res = cIdealSolidSolnPhase0;
        break;
    case 1:
        res = cIdealSolidSolnPhase1;
        break;
    case 2:
        res = cIdealSolidSolnPhase2;
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
    return res;
}

/********************************************************************
 *            Molar Thermodynamic Properties of the Solution
 ********************************************************************/

doublereal IdealSolidSolnPhase::enthalpy_mole() const
{
    doublereal htp = GasConstant * temperature() * mean_X(enthalpy_RT_ref());
    return htp + (pressure() - m_Pref)/molarDensity();
}

doublereal IdealSolidSolnPhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx());
}

doublereal IdealSolidSolnPhase::gibbs_mole() const
{
    return GasConstant * temperature() * (mean_X(gibbs_RT_ref()) + sum_xlogx());
}

doublereal IdealSolidSolnPhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

/********************************************************************
 *                  Mechanical Equation of State
 ********************************************************************/

void IdealSolidSolnPhase::calcDensity()
{
    /*
     * Calculate the molarVolume of the solution (m**3 kmol-1)
     */
    const doublereal* const dtmp = moleFractdivMMW();
    double invDens = dot(m_speciesMolarVolume.begin(),
                         m_speciesMolarVolume.end(), dtmp);
    /*
     * Set the density in the parent State object directly,
     * by calling the Phase::setDensity() function.
     */
    Phase::setDensity(1.0/invDens);
}

void IdealSolidSolnPhase::setDensity(const doublereal rho)
{
    /*
     * Unless the input density is exactly equal to the density
     * calculated and stored in the State object, we throw an
     * exception. This is because the density is NOT an
     * independent variable.
     */
    if (rho != density()) {
        throw CanteraError("IdealSolidSolnPhase::setDensity",
                           "Density is not an independent variable");
    }
}

void IdealSolidSolnPhase::setPressure(doublereal p)
{
    m_Pcurrent = p;
    calcDensity();
}

void IdealSolidSolnPhase::setMolarDensity(const doublereal n)
{
    throw CanteraError("IdealSolidSolnPhase::setMolarDensity",
                       "Density is not an independent variable");
}

void IdealSolidSolnPhase::setMoleFractions(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    calcDensity();
}

void IdealSolidSolnPhase::setMoleFractions_NoNorm(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    calcDensity();
}

void IdealSolidSolnPhase::setMassFractions(const doublereal* const y)
{
    Phase::setMassFractions(y);
    calcDensity();
}

void IdealSolidSolnPhase::setMassFractions_NoNorm(const doublereal* const y)
{
    Phase::setMassFractions_NoNorm(y);
    calcDensity();
}

void IdealSolidSolnPhase::setConcentrations(const doublereal* const c)
{
    Phase::setConcentrations(c);
    calcDensity();
}

/********************************************************************
 *        Chemical Potentials and Activities
 ********************************************************************/

void IdealSolidSolnPhase::getActivityConcentrations(doublereal* c) const
{
    const doublereal* const dtmp = moleFractdivMMW();
    const double mmw = meanMolecularWeight();
    switch (m_formGC) {
    case 0:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = dtmp[k] * mmw;
        }
        break;
    case 1:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = dtmp[k] * mmw / m_speciesMolarVolume[k];
        }
        break;
    case 2:
        double atmp = mmw / m_speciesMolarVolume[m_kk-1];
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = dtmp[k] * atmp;
        }
        break;
    }
}

doublereal IdealSolidSolnPhase::standardConcentration(size_t k) const
{
    switch (m_formGC) {
    case 0:
        return 1.0;
    case 1:
        return 1.0 / m_speciesMolarVolume[k];
    case 2:
        return 1.0/m_speciesMolarVolume[m_kk-1];
    }
    return 0.0;
}
doublereal IdealSolidSolnPhase::referenceConcentration(int k) const
{
    switch (m_formGC) {
    case 0:
        return 1.0;
    case 1:
        return 1.0 / m_speciesMolarVolume[k];
    case 2:
        return 1.0 / m_speciesMolarVolume[m_kk-1];
    }
    return 0.0;
}

doublereal IdealSolidSolnPhase::logStandardConc(size_t k) const
{
    _updateThermo();
    double res;
    switch (m_formGC) {
    case 0:
        res = 0.0;
        break;
    case 1:
        res = log(1.0/m_speciesMolarVolume[k]);
        break;
    case 2:
        res =  log(1.0/m_speciesMolarVolume[m_kk-1]);
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
    return res;
}

void IdealSolidSolnPhase::getUnitsStandardConc(double* uA, int, int sizeUA) const
{
    warn_deprecated("IdealSolidSolnPhase::getUnitsStandardConc",
                    "To be removed after Cantera 2.2.");
    int eos = eosType();
    if (eos == cIdealSolidSolnPhase0) {
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

void IdealSolidSolnPhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

void IdealSolidSolnPhase::getChemPotentials(doublereal* mu) const
{
    doublereal delta_p = m_Pcurrent - m_Pref;
    doublereal RT = temperature() * GasConstant;
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT * (g_RT[k] + log(xx))
                + delta_p * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getChemPotentials_RT(doublereal* mu) const
{
    doublereal delta_pdRT = (m_Pcurrent - m_Pref) / (temperature() * GasConstant);
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = (g_RT[k] + log(xx))
                + delta_pdRT * m_speciesMolarVolume[k];
    }
}

/********************************************************************
 *                    Partial Molar Properties
 ********************************************************************/

void IdealSolidSolnPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    scale(_h.begin(), _h.end(), hbar, GasConstant * temperature());
}

void IdealSolidSolnPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] = GasConstant * (_s[k] - log(xx));
    }
}

void IdealSolidSolnPhase::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void IdealSolidSolnPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

/*****************************************************************
 * Properties of the Standard State of the Species in the Solution
 *****************************************************************/

void IdealSolidSolnPhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal RT = _RT();
    const doublereal* const gk = DATA_PTR(gibbsrt);
    doublereal delta_p = (m_Pcurrent - m_Pref);
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = RT * gk[k] + delta_p * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal RT = _RT();
    const doublereal* const gk = DATA_PTR(gibbsrt);
    doublereal delta_prt = (m_Pcurrent - m_Pref)/ RT;
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = gk[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal delta_prt = ((m_Pcurrent - m_Pref) /
                            (GasConstant * temperature()));
    for (size_t k = 0; k < m_kk; k++) {
        hrt[k] = _h[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEntropy_R(doublereal* sr) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
}

void IdealSolidSolnPhase::getIntEnergy_RT(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal prefrt = m_Pref / (GasConstant * temperature());
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - prefrt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getCp_R(doublereal* cpr) const
{
    const vector_fp& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

void IdealSolidSolnPhase::getStandardVolumes(doublereal* vol) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), vol);
}

/*********************************************************************
 *     Thermodynamic Values for the Species Reference States
 *********************************************************************/

void IdealSolidSolnPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        hrt[k] = m_h0_RT[k];
    }
}

void IdealSolidSolnPhase::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        grt[k] = m_g0_RT[k];
    }
}

void IdealSolidSolnPhase::getGibbs_ref(doublereal* g) const
{
    _updateThermo();
    double tmp = GasConstant * temperature();
    for (size_t k = 0; k != m_kk; k++) {
        g[k] = tmp * m_g0_RT[k];
    }
}

void IdealSolidSolnPhase::getIntEnergy_RT_ref(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal prefrt = m_Pref / (GasConstant * temperature());
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - prefrt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEntropy_R_ref(doublereal* er) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        er[k] = m_s0_R[k];
    }
}

void IdealSolidSolnPhase::getCp_R_ref(doublereal* cpr) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        cpr[k] = m_cp0_R[k];
    }
}

const vector_fp& IdealSolidSolnPhase::enthalpy_RT_ref() const
{
    _updateThermo();
    return m_h0_RT;
}

const vector_fp& IdealSolidSolnPhase::entropy_R_ref() const
{
    _updateThermo();
    return m_s0_R;
}

/*********************************************************************
 *    Utility Functions
 *********************************************************************/

void IdealSolidSolnPhase::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0) {
        if (phaseNode.id() != id_) {
            throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Check on the thermo field. Must have:
     * <thermo model="IdealSolidSolution" />
     */
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thNode = phaseNode.child("thermo");
        string mString = thNode.attrib("model");
        if (lowercase(mString) != "idealsolidsolution") {
            throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                               "Unknown thermo model: " + mString);
        }
    } else {
        throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                           "Unspecified thermo model");
    }

    /*
     * Form of the standard concentrations. Must have one of:
     *
     *     <standardConc model="unity" />
     *     <standardConc model="molar_volume" />
     *     <standardConc model="solvent_volume" />
     */
    if (phaseNode.hasChild("standardConc")) {
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
            throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                               "Unknown standardConc model: " + formStringa);
        }
    } else {
        throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                           "Unspecified standardConc model");
    }

    /*
     * Initialize all of the lengths now that we know how many species
     * there are in the phase.
     */
    initLengths();
    /*
     * Now go get the molar volumes
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &phaseNode.root());

    for (size_t k = 0; k < m_kk; k++) {
        XML_Node* s =  speciesDB->findByAttr("name", speciesName(k));
        XML_Node* ss = s->findByName("standardState");
        m_speciesMolarVolume[k] = getFloat(*ss, "molarVolume", "toSI");
    }

    /*
     * Call the base initThermo, which handles setting the initial
     * state.
     */
    ThermoPhase::initThermoXML(phaseNode, id_);
}

void IdealSolidSolnPhase::initLengths()
{
    /*
     * Obtain the reference pressure by calling the ThermoPhase
     * function refPressure, which in turn calls the
     * species thermo reference pressure function of the
     * same name.
     */
    m_Pref = refPressure();

    m_h0_RT.resize(m_kk);
    m_g0_RT.resize(m_kk);
    m_expg0_RT.resize(m_kk);
    m_cp0_R.resize(m_kk);
    m_s0_R.resize(m_kk);
    m_pe.resize(m_kk, 0.0);
    m_pp.resize(m_kk);
    m_speciesMolarVolume.resize(m_kk);
}

void IdealSolidSolnPhase::setToEquilState(const doublereal* lambda_RT)
{
    const vector_fp& grt = gibbs_RT_ref();

    // set the pressure and composition to be consistent with
    // the temperature,
    doublereal pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = -grt[k];
        for (size_t m = 0; m < nElements(); m++) {
            m_pp[k] += nAtoms(k,m)*lambda_RT[m];
        }
        m_pp[k] = m_Pref * exp(m_pp[k]);
        pres += m_pp[k];
    }
    setState_PX(pres, &m_pp[0]);
}

double IdealSolidSolnPhase::speciesMolarVolume(int k) const
{
    return  m_speciesMolarVolume[k];
}

void IdealSolidSolnPhase::getSpeciesMolarVolumes(doublereal* smv) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), smv);
}

void IdealSolidSolnPhase::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        /*
         * Update the thermodynamic functions of the reference state.
         */
        m_spthermo->update(tnow, DATA_PTR(m_cp0_R), DATA_PTR(m_h0_RT),
                           DATA_PTR(m_s0_R));
        m_tlast = tnow;
        doublereal rrt = 1.0 / (GasConstant * tnow);
        for (size_t k = 0; k < m_kk; k++) {
            double deltaE = rrt * m_pe[k];
            m_h0_RT[k] += deltaE;
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

}  // end namespace Cantera
