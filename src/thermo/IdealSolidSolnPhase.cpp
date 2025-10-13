/**
 *  @file IdealSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see @ref
 *      thermoprops and @link Cantera::IdealSolidSolnPhase
 *      IdealSolidSolnPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

IdealSolidSolnPhase::IdealSolidSolnPhase(const string& inputFile, const string& id_)
{
    initThermoFile(inputFile, id_);
}

// Molar Thermodynamic Properties of the Solution

double IdealSolidSolnPhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx());
}

double IdealSolidSolnPhase::gibbs_mole() const
{
    double Pv = (pressure() - m_Pref)/molarDensity();
    return RT() * (mean_X(gibbs_RT_ref()) + sum_xlogx()) + Pv;
}

double IdealSolidSolnPhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

// Mechanical Equation of State

void IdealSolidSolnPhase::calcDensity()
{
    // Calculate the molarVolume of the solution (m**3 kmol-1)
    double v_mol = mean_X(m_speciesMolarVolume);

    // Set the density in the parent object directly, by calling the
    // Phase::assignDensity() function.
    Phase::assignDensity(meanMolecularWeight()/v_mol);
}

void IdealSolidSolnPhase::setPressure(double p)
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

void IdealSolidSolnPhase::getActivityConcentrations(double* c) const
{
    getMoleFractions(c);
    switch (m_formGC) {
    case 0:
        break;
    case 1:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] /= m_speciesMolarVolume[k];
        }
        break;
    case 2:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] /= m_speciesMolarVolume[m_kk-1];
        }
        break;
    }
}

double IdealSolidSolnPhase::standardConcentration(size_t k) const
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

void IdealSolidSolnPhase::getActivityCoefficients(double* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

void IdealSolidSolnPhase::getChemPotentials(double* mu) const
{
    double delta_p = m_Pcurrent - m_Pref;
    const vector<double>& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT() * (g_RT[k] + log(xx))
                + delta_p * m_speciesMolarVolume[k];
    }
}

// Partial Molar Properties

void IdealSolidSolnPhase::getPartialMolarEnthalpies(double* hbar) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    double delta_p = m_Pcurrent - m_Pref;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] = _h[k]*RT() + delta_p * m_speciesMolarVolume[k];
    }
    // scale(_h.begin(), _h.end(), hbar, RT());
}

void IdealSolidSolnPhase::getPartialMolarEntropies(double* sbar) const
{
    const vector<double>& _s = entropy_R_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] = GasConstant * (_s[k] - log(xx));
    }
}

void IdealSolidSolnPhase::getPartialMolarCp(double* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void IdealSolidSolnPhase::getPartialMolarVolumes(double* vbar) const
{
    getStandardVolumes(vbar);
}

// Properties of the Standard State of the Species in the Solution

void IdealSolidSolnPhase::getPureGibbs(double* gpure) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    double delta_p = (m_Pcurrent - m_Pref);
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = RT() * gibbsrt[k] + delta_p * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getGibbs_RT(double* grt) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    double delta_prt = (m_Pcurrent - m_Pref)/ RT();
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = gibbsrt[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEnthalpy_RT(double* hrt) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    double delta_prt = (m_Pcurrent - m_Pref) / RT();
    for (size_t k = 0; k < m_kk; k++) {
        hrt[k] = _h[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEntropy_R(double* sr) const
{
    const vector<double>& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
}

void IdealSolidSolnPhase::getIntEnergy_RT(double* urt) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    double prefrt = m_Pref / RT();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - prefrt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getCp_R(double* cpr) const
{
    const vector<double>& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

void IdealSolidSolnPhase::getStandardVolumes(double* vol) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), vol);
}

// Thermodynamic Values for the Species Reference States

void IdealSolidSolnPhase::getEnthalpy_RT_ref(double* hrt) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        hrt[k] = m_h0_RT[k];
    }
}

void IdealSolidSolnPhase::getGibbs_RT_ref(double* grt) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        grt[k] = m_g0_RT[k];
    }
}

void IdealSolidSolnPhase::getGibbs_ref(double* g) const
{
    _updateThermo();
    double tmp = RT();
    for (size_t k = 0; k != m_kk; k++) {
        g[k] = tmp * m_g0_RT[k];
    }
}

void IdealSolidSolnPhase::getIntEnergy_RT_ref(double* urt) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    double prefrt = m_Pref / RT();
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - prefrt * m_speciesMolarVolume[k];
    }
}

void IdealSolidSolnPhase::getEntropy_R_ref(double* er) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        er[k] = m_s0_R[k];
    }
}

void IdealSolidSolnPhase::getCp_R_ref(double* cpr) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        cpr[k] = m_cp0_R[k];
    }
}

const vector<double>& IdealSolidSolnPhase::enthalpy_RT_ref() const
{
    _updateThermo();
    return m_h0_RT;
}

const vector<double>& IdealSolidSolnPhase::entropy_R_ref() const
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

void IdealSolidSolnPhase::getSpeciesParameters(const string &name,
                                               AnyMap& speciesNode) const
{
    ThermoPhase::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name, true);
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

void IdealSolidSolnPhase::setToEquilState(const double* mu_RT)
{
    const vector<double>& grt = gibbs_RT_ref();

    // Within the method, we protect against inf results if the exponent is too
    // high.
    //
    // If it is too low, we set the partial pressure to zero. This capability is
    // needed by the elemental potential method.
    double pres = 0.0;
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
    setMoleFractions(m_pp.data());
    setPressure(pres);
}

void IdealSolidSolnPhase::setStandardConcentrationModel(const string& model)
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

void IdealSolidSolnPhase::getSpeciesMolarVolumes(double* smv) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), smv);
}

void IdealSolidSolnPhase::_updateThermo() const
{
    double tnow = temperature();
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
