/**
 *  @file LatticePhase.cpp
 *  Definitions for a simple thermodynamics model of a bulk phase
 *  derived from ThermoPhase,
 *  assuming a lattice of solid atoms
 *  (see @ref thermoprops and class @link Cantera::LatticePhase LatticePhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/LatticePhase.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/Species.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

namespace Cantera
{

LatticePhase::LatticePhase(const string& inputFile, const string& id_)
{
    warn_deprecated("class LatticePhase", "To be removed after Cantera 3.2. Can be "
        "replaced by use of IdealSolidSolnPhase with the site density used to set the "
        "molar density of each constituent species.");
    initThermoFile(inputFile, id_);
}

double LatticePhase::enthalpy_mole() const
{
    return RT() * mean_X(enthalpy_RT_ref()) +
            (pressure() - m_Pref)/molarDensity();
}

double LatticePhase::entropy_mole() const
{
    return GasConstant * (mean_X(entropy_R_ref()) - sum_xlogx());
}

double LatticePhase::gibbs_mole() const {
    return enthalpy_mole() - temperature() * entropy_mole();
}

double LatticePhase::cp_mole() const
{
    return GasConstant * mean_X(cp_R_ref());
}

double LatticePhase::cv_mole() const
{
    return cp_mole();
}

double LatticePhase::calcDensity()
{
    assignDensity(std::max(meanMolecularWeight() * m_site_density, SmallNumber));
    return meanMolecularWeight() * m_site_density;
}

void LatticePhase::setPressure(double p)
{
    m_Pcurrent = p;
    calcDensity();
}

void LatticePhase::compositionChanged()
{
    Phase::compositionChanged();
    calcDensity();
}

Units LatticePhase::standardConcentrationUnits() const
{
    return Units(1.0);
}

void LatticePhase::getActivityConcentrations(double* c) const
{
    getMoleFractions(c);
}

void LatticePhase::getActivityCoefficients(double* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}

double LatticePhase::standardConcentration(size_t k) const
{
    return 1.0;
}

double LatticePhase::logStandardConc(size_t k) const
{
    return 0.0;
}

void LatticePhase::getChemPotentials(double* mu) const
{
    double delta_p = m_Pcurrent - m_Pref;
    const vector<double>& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT() * (g_RT[k] + log(xx))
                + delta_p * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getPartialMolarEnthalpies(double* hbar) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    scale(_h.begin(), _h.end(), hbar, RT());
}

void LatticePhase::getPartialMolarEntropies(double* sbar) const
{
    const vector<double>& _s = entropy_R_ref();
    for (size_t k = 0; k < m_kk; k++) {
        double xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] = GasConstant * (_s[k] - log(xx));
    }
}

void LatticePhase::getPartialMolarCp(double* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

void LatticePhase::getPartialMolarVolumes(double* vbar) const
{
    getStandardVolumes(vbar);
}

void LatticePhase::getStandardChemPotentials(double* mu0) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    scale(gibbsrt.begin(), gibbsrt.end(), mu0, RT());
}

void LatticePhase::getPureGibbs(double* gpure) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    double delta_p = (m_Pcurrent - m_Pref);
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = RT() * gibbsrt[k] + delta_p * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getEnthalpy_RT(double* hrt) const
{
    const vector<double>& _h = enthalpy_RT_ref();
    double delta_prt = (m_Pcurrent - m_Pref) / RT();
    for (size_t k = 0; k < m_kk; k++) {
        hrt[k] = _h[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getEntropy_R(double* sr) const
{
    const vector<double>& _s = entropy_R_ref();
    std::copy(_s.begin(), _s.end(), sr);
}

void LatticePhase::getGibbs_RT(double* grt) const
{
    const vector<double>& gibbsrt = gibbs_RT_ref();
    double delta_prt = (m_Pcurrent - m_Pref) / RT();
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = gibbsrt[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

void LatticePhase::getGibbs_ref(double* g) const
{
    getGibbs_RT_ref(g);
    for (size_t k = 0; k < m_kk; k++) {
        g[k] *= RT();
    }
}

void LatticePhase::getCp_R(double* cpr) const
{
    const vector<double>& _cpr = cp_R_ref();
    std::copy(_cpr.begin(), _cpr.end(), cpr);
}

void LatticePhase::getStandardVolumes(double* vbar) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), vbar);
}

const vector<double>& LatticePhase::enthalpy_RT_ref() const
{
    _updateThermo();
    return m_h0_RT;
}

const vector<double>& LatticePhase::gibbs_RT_ref() const
{
    _updateThermo();
    return m_g0_RT;
}

void LatticePhase::getGibbs_RT_ref(double* grt) const
{
    _updateThermo();
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = m_g0_RT[k];
    }
}

const vector<double>& LatticePhase::entropy_R_ref() const
{
    _updateThermo();
    return m_s0_R;
}

const vector<double>& LatticePhase::cp_R_ref() const
{
    _updateThermo();
    return m_cp0_R;
}

bool LatticePhase::addSpecies(shared_ptr<Species> spec)
{
    bool added = ThermoPhase::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            m_Pref = refPressure();
        }
        m_h0_RT.push_back(0.0);
        m_g0_RT.push_back(0.0);
        m_cp0_R.push_back(0.0);
        m_s0_R.push_back(0.0);
        double mv = 1.0 / m_site_density;
        if (spec->input.hasKey("equation-of-state")) {
            auto& eos = spec->input["equation-of-state"].getMapWhere(
                "model", "constant-volume");
            if (eos.hasKey("density")) {
                mv = molecularWeight(m_kk-1) / eos.convert("density", "kg/m^3");
            } else if (eos.hasKey("molar-density")) {
                mv = 1.0 / eos.convert("molar-density", "kmol/m^3");
            } else if (eos.hasKey("molar-volume")) {
                mv = eos.convert("molar-volume", "m^3/kmol");
            }
        }
        m_speciesMolarVolume.push_back(mv);
    }
    return added;
}

void LatticePhase::setSiteDensity(double sitedens)
{
    m_site_density = sitedens;
    for (size_t k = 0; k < m_kk; k++) {
        if (species(k)->input.hasKey("equation-of-state")) {
            auto& eos = species(k)->input["equation-of-state"].getMapWhere(
                "model", "constant-volume");
            if (eos.hasKey("molar-volume") || eos.hasKey("density")
                || eos.hasKey("molar-density")) {
                continue;
            }
        }
        m_speciesMolarVolume[k] = 1.0 / m_site_density;
    }
}

void LatticePhase::_updateThermo() const
{
    double tnow = temperature();
    if (m_tlast != tnow) {
        m_spthermo.update(tnow, &m_cp0_R[0], &m_h0_RT[0], &m_s0_R[0]);
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

void LatticePhase::initThermo()
{
    if (m_input.hasKey("site-density")) {
        setSiteDensity(m_input.convert("site-density", "kmol/m^3"));
    }
}

void LatticePhase::getParameters(AnyMap& phaseNode) const
{
    ThermoPhase::getParameters(phaseNode);
    phaseNode["site-density"].setQuantity(m_site_density, "kmol/m^3");
}

void LatticePhase::getSpeciesParameters(const string& name, AnyMap& speciesNode) const
{
    ThermoPhase::getSpeciesParameters(name, speciesNode);
    size_t k = speciesIndex(name);
    // Output volume information in a form consistent with the input
    const auto S = species(k);
    if (S->input.hasKey("equation-of-state")) {
        auto& eosIn = S->input["equation-of-state"].getMapWhere(
            "model", "constant-volume");
        auto& eosOut = speciesNode["equation-of-state"].getMapWhere(
            "model", "constant-volume", true);

        if (eosIn.hasKey("density")) {
            eosOut["model"] = "constant-volume";
            eosOut["density"].setQuantity(
                molecularWeight(k) / m_speciesMolarVolume[k], "kg/m^3");
        } else if (eosIn.hasKey("molar-density")) {
            eosOut["model"] = "constant-volume";
            eosOut["molar-density"].setQuantity(1.0 / m_speciesMolarVolume[k],
                                                "kmol/m^3");
        } else if (eosIn.hasKey("molar-volume")) {
            eosOut["model"] = "constant-volume";
            eosOut["molar-volume"].setQuantity(m_speciesMolarVolume[k],
                                               "m^3/kmol");
        }
    }
    // Otherwise, species volume is determined by the phase-level site density
}

}
