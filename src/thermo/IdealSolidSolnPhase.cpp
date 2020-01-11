/**
 *  @file IdealSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see \ref
 *      thermoprops and \link Cantera::IdealSolidSolnPhase
 *      IdealSolidSolnPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{

IdealSolidSolnPhase::IdealSolidSolnPhase(int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError("IdealSolidSolnPhase::IdealSolidSolnPhase",
                           "Illegal value of formGC");
    }
}

IdealSolidSolnPhase::IdealSolidSolnPhase(const std::string& inputFile,
        const std::string& id_, int formGC) :
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError("IdealSolidSolnPhase::IdealSolidSolnPhase",
                           "Illegal value of formGC");
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
        throw CanteraError("IdealSolidSolnPhase::IdealSolidSolnPhase",
                           "Illegal value of formGC");
    }
    importPhase(root, this);
}

// Molar Thermodynamic Properties of the Solution

doublereal IdealSolidSolnPhase::enthalpy_mole() const
{
    doublereal htp = RT() * mean_X(enthalpy_RT_ref());
    return htp + (pressure() - m_Pref)/molarDensity();
}

doublereal IdealSolidSolnPhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx());
}

doublereal IdealSolidSolnPhase::gibbs_mole() const
{
    return RT() * (mean_X(gibbs_RT_ref()) + sum_xlogx());
}

doublereal IdealSolidSolnPhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

// Mechanical Equation of State

void IdealSolidSolnPhase::calcDensity()
{
    // Calculate the molarVolume of the solution (m**3 kmol-1)
    const doublereal* const dtmp = moleFractdivMMW();
    double invDens = dot(m_speciesMolarVolume.begin(),
                         m_speciesMolarVolume.end(), dtmp);

    // Set the density in the parent State object directly, by calling the
    // Phase::assignDensity() function.
    Phase::assignDensity(1.0/invDens);
}

void IdealSolidSolnPhase::setDensity(const doublereal rho)
{
    // Unless the input density is exactly equal to the density calculated and
    // stored in the State object, we throw an exception. This is because the
    // density is NOT an independent variable.
    warn_deprecated("IdealSolidSolnPhase::setDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setDensity instead");
    if (std::abs(rho/density() - 1.0) > 1e-15) {
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
    warn_deprecated("IdealSolidSolnPhase::setMolarDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setMolarDensity instead");
    throw CanteraError("IdealSolidSolnPhase::setMolarDensity",
                       "Density is not an independent variable");
}

void IdealSolidSolnPhase::compositionChanged()
{
    Phase::compositionChanged();
    calcDensity();
}

// Chemical Potentials and Activities

Units IdealSolidSolnPhase::standardConcentrationUnits() const
{
    if (m_formGC == 0) {
        return Units(1.0); // dimensionless
    } else {
        // kmol/m^3 for bulk phases
        return Units(1.0, 0, -static_cast<double>(nDim()), 0, 0, 0, 1);
    }
}

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

void IdealSolidSolnPhase::getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

void IdealSolidSolnPhase::getChemPotentials(doublereal* mu) const
{
    doublereal delta_p = m_Pcurrent - m_Pref;
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT() * (g_RT[k] + log(xx))
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

// Partial Molar Properties

void IdealSolidSolnPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    scale(_h.begin(), _h.end(), hbar, RT());
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

// Properties of the Standard State of the Species in the Solution

void IdealSolidSolnPhase::getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal delta_p = (m_Pcurrent - m_Pref);
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = RT() * gibbsrt[k] + delta_p * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal delta_prt = (m_Pcurrent - m_Pref)/ RT();
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = gibbsrt[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal delta_prt = (m_Pcurrent - m_Pref) / RT();
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
    doublereal prefrt = m_Pref / RT();
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

// Thermodynamic Values for the Species Reference States

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
    double tmp = RT();
    for (size_t k = 0; k != m_kk; k++) {
        g[k] = tmp * m_g0_RT[k];
    }
}

void IdealSolidSolnPhase::getIntEnergy_RT_ref(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal prefrt = m_Pref / RT();
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

// Utility Functions

bool IdealSolidSolnPhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            // Obtain the reference pressure by calling the ThermoPhase function
            // refPressure, which in turn calls the species thermo reference
            // pressure function of the same name.
            m_Pref = refPressure();
        }

        m_h0_RT.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_expg0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_s0_R.push_back(0.0);
        m_pe.push_back(0.0);;
        m_pp.push_back(0.0);
        if (spec->input.hasKey("equation-of-state")) {
            auto& eos = spec->input["equation-of-state"].getMapWhere("model", "constant-volume");
            double mv;
            if (eos.hasKey("density")) {
                mv = molecularWeight(m_kk-1) / eos.convert("density", "kg/m^3");
            } else if (eos.hasKey("molar-density")) {
                mv = 1.0 / eos.convert("molar-density", "kmol/m^3");
            } else if (eos.hasKey("molar-volume")) {
                mv = eos.convert("molar-volume", "m^3/kmol");
            } else {
                throw CanteraError("IdealSolidSolnPhase::addSpecies",
                    "equation-of-state entry for species '{}' is missing "
                    "'density', 'molar-volume', or 'molar-density' "
                    "specification", spec->name);
            }
            m_speciesMolarVolume.push_back(mv);
        } else if (spec->extra.hasKey("molar_volume")) {
            // From XML
            m_speciesMolarVolume.push_back(spec->extra["molar_volume"].asDouble());
        } else {
            throw CanteraError("IdealSolidSolnPhase::addSpecies",
                "Molar volume not specified for species '{}'", spec->name);
        }
        calcDensity();
    }
    return added;
}

void IdealSolidSolnPhase::initThermo()
{
   if (m_input.hasKey("standard-concentration-basis")) {
        setStandardConcentrationModel(m_input["standard-concentration-basis"].asString());
    }
    ThermoPhase::initThermo();
}

void IdealSolidSolnPhase::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0 && phaseNode.id() != id_) {
        throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    // Check on the thermo field. Must have:
    // <thermo model="IdealSolidSolution" />
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thNode = phaseNode.child("thermo");
        if (!caseInsensitiveEquals(thNode["model"], "idealsolidsolution")) {
            throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                               "Unknown thermo model: " + thNode["model"]);
        }
    } else {
        throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                           "Unspecified thermo model");
    }

    // Form of the standard concentrations. Must have one of:
    //
    //     <standardConc model="unity" />
    //     <standardConc model="molar_volume" />
    //     <standardConc model="solvent_volume" />
    if (phaseNode.hasChild("standardConc")) {
        setStandardConcentrationModel(phaseNode.child("standardConc")["model"]);
    } else {
        throw CanteraError("IdealSolidSolnPhase::initThermoXML",
                           "Unspecified standardConc model");
    }

    // Call the base initThermo, which handles setting the initial state.
    ThermoPhase::initThermoXML(phaseNode, id_);
}

void IdealSolidSolnPhase::setToEquilState(const doublereal* lambda_RT)
{
    const vector_fp& grt = gibbs_RT_ref();

    // set the pressure and composition to be consistent with the temperature
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

void IdealSolidSolnPhase::setStandardConcentrationModel(const std::string& model)
{
    if (caseInsensitiveEquals(model, "unity")) {
        m_formGC = 0;
    } else if (caseInsensitiveEquals(model, "species-molar-volume")
               || caseInsensitiveEquals(model, "molar_volume")) {
        m_formGC = 1;
    } else if (caseInsensitiveEquals(model, "solvent-molar-volume")
               || caseInsensitiveEquals(model, "solvent_volume")) {
        m_formGC = 2;
    } else {
        throw CanteraError("IdealSolidSolnPhase::setStandardConcentrationModel",
                           "Unknown standard concentration model '{}'", model);
    }
}

double IdealSolidSolnPhase::speciesMolarVolume(int k) const
{
    return m_speciesMolarVolume[k];
}

void IdealSolidSolnPhase::getSpeciesMolarVolumes(doublereal* smv) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), smv);
}

void IdealSolidSolnPhase::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {

        // Update the thermodynamic functions of the reference state.
        m_spthermo.update(tnow, m_cp0_R.data(), m_h0_RT.data(), m_s0_R.data());
        m_tlast = tnow;
        doublereal rrt = 1.0 / RT();
        for (size_t k = 0; k < m_kk; k++) {
            double deltaE = rrt * m_pe[k];
            m_h0_RT[k] += deltaE;
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

} // end namespace Cantera
