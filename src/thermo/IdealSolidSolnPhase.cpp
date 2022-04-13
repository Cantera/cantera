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
#include "cantera/thermo/Species.h"
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
    // @todo: After Cantera 2.6, this constructor can be deleted and the default
    // construction option can be provided by adding "" as the default argument
    // to the constructor from input file name and phase id.
    if (formGC == -1) {
        formGC = 0;
    } else {
        warn_deprecated("IdealSolidSolnPhase(int formGC)",
            "The formGC constructor argument is deprecated and will be removed"
            " after Cantera 2.6. Use the setStandardConcentrationModel"
            " method instead.");
    }
    m_formGC = formGC;
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
    if (formGC == -1) {
        formGC = 0;
    } else {
        warn_deprecated("IdealSolidSolnPhase(string inputFile, string id_, int formGC)",
            "The formGC constructor argument is deprecated and will be removed"
            " after Cantera 2.6. Use the setStandardConcentrationModel"
            " method instead.");
    }
    m_formGC = formGC;
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
    if (formGC == -1) {
        formGC = 0;
    } else {
        warn_deprecated("IdealSolidSolnPhase(XML_Node root, string id_, int formGC)",
            "The formGC constructor argument is deprecated and will be removed"
            " after Cantera 2.6. Use the setStandardConcentrationModel"
            " method instead.");
    }
    m_formGC = formGC;
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

void IdealSolidSolnPhase::setPressure(doublereal p)
{
    m_Pcurrent = p;
    calcDensity();
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
    double delta_p = m_Pcurrent - m_Pref;
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT() * (g_RT[k] + log(xx))
                + delta_p * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getChemPotentials_RT(doublereal* mu) const
{
    double delta_pdRT = (m_Pcurrent - m_Pref) / (temperature() * GasConstant);
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
    double delta_p = m_Pcurrent - m_Pref;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] = _h[k]*RT() + delta_p * m_speciesMolarVolume[k];
    }
    // scale(_h.begin(), _h.end(), hbar, RT());
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
    double delta_p = (m_Pcurrent - m_Pref);
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
        } else if (spec->input.hasKey("molar_volume")) {
            // @Deprecated - remove this case for Cantera 3.0 with removal of the XML format
            m_speciesMolarVolume.push_back(spec->input["molar_volume"].asDouble());
        } else {
            throw CanteraError("IdealSolidSolnPhase::addSpecies",
                "Molar volume not specified for species '{}'", spec->name);
        }
        if (ready()) {
            calcDensity();
        }
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

void IdealSolidSolnPhase::getParameters(AnyMap& phaseNode) const
{
    ThermoPhase::getParameters(phaseNode);
    if (m_formGC == 1) {
        phaseNode["standard-concentration-basis"] = "species-molar-volume";
    } else if (m_formGC == 2) {
        phaseNode["standard-concentration-basis"] = "solvent-molar-volume";
    }
}

void IdealSolidSolnPhase::getSpeciesParameters(const std::string &name,
                                               AnyMap& speciesNode) const
{
    ThermoPhase::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name);
    const auto S = species(k);
    auto& eosNode = speciesNode["equation-of-state"].getMapWhere(
        "model", "constant-volume", true);
    // Output volume information in a form consistent with the input
    if (S->input.hasKey("equation-of-state")) {
        auto& eosIn = S->input["equation-of-state"];
        if (eosIn.hasKey("density")) {
            eosNode["density"].setQuantity(
                molecularWeight(k) / m_speciesMolarVolume[k], "kg/m^3");
        } else if (eosIn.hasKey("molar-density")) {
            eosNode["molar-density"].setQuantity(1.0 / m_speciesMolarVolume[k],
                                                 "kmol/m^3");
        } else {
            eosNode["molar-volume"].setQuantity(m_speciesMolarVolume[k],
                                                "m^3/kmol");
        }
    } else {
        eosNode["molar-volume"].setQuantity(m_speciesMolarVolume[k],
                                            "m^3/kmol");
    }
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

void IdealSolidSolnPhase::setToEquilState(const doublereal* mu_RT)
{
    const vector_fp& grt = gibbs_RT_ref();
    
    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    doublereal pres = 0.0;
    double m_p0 = refPressure();
    for (size_t k = 0; k < m_kk; k++) {
        double tmp = -grt[k] + mu_RT[k];
        if (tmp < -600.) {
            m_pp[k] = 0.0;
        } else if (tmp > 500.0) {
            // Protect against inf results if the exponent is too high
            double tmp2 = tmp / 500.;
            tmp2 *= tmp2;
            m_pp[k] = m_p0 * exp(500.) * tmp2;
        } else {
            m_pp[k] = m_p0 * exp(tmp);
        }
        pres += m_pp[k];
    }
    // set state
    setState_PX(pres, m_pp.data());
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
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

} // end namespace Cantera
