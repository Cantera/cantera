/**
 *  @file SingleSpeciesTP.cpp
 *  Definitions for the SingleSpeciesTP class, which is a filter class for ThermoPhase,
 *  that eases the construction of single species phases
 *  ( see @ref thermoprops and class @link Cantera::SingleSpeciesTP SingleSpeciesTP@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/SingleSpeciesTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

namespace Cantera
{

// ------------ Molar Thermodynamic Properties --------------------

double SingleSpeciesTP::enthalpy_mole() const
{
    double hbar;
    getPartialMolarEnthalpies(&hbar);
    return hbar;
}

double SingleSpeciesTP::intEnergy_mole() const
{
    double ubar;
    getPartialMolarIntEnergies(&ubar);
    return ubar;
}

double SingleSpeciesTP::entropy_mole() const
{
    double sbar;
    getPartialMolarEntropies(&sbar);
    return sbar;
}

double SingleSpeciesTP::gibbs_mole() const
{
    double gbar;

    // Get the chemical potential of the first species. This is the same as the
    // partial molar Gibbs free energy.
    getChemPotentials(&gbar);
    return gbar;
}

double SingleSpeciesTP::cp_mole() const
{
    double cpbar;

    // Really should have a partial molar heat capacity function in ThermoPhase.
    // However, the standard state heat capacity will do fine here for now.
    getCp_R(&cpbar);
    cpbar *= GasConstant;
    return cpbar;
}

double SingleSpeciesTP::cv_mole() const
{
    // For single species, we go directory to the general Cp - Cv relation
    //
    //     Cp = Cv + alpha**2 * V * T / beta
    //
    // where
    //     alpha = volume thermal expansion coefficient
    //     beta  = isothermal compressibility
    double cvbar = cp_mole();
    double alpha = thermalExpansionCoeff();
    double beta = isothermalCompressibility();
    double V = molecularWeight(0)/density();
    double T = temperature();
    if (beta != 0.0) {
        cvbar -= alpha * alpha * V * T / beta;
    }
    return cvbar;
}

// ----------- Partial Molar Properties of the Solution -----------------

void SingleSpeciesTP::getChemPotentials(double* mu) const
{
    getStandardChemPotentials(mu);
}

void SingleSpeciesTP::getPartialMolarEnthalpies(double* hbar) const
{
    getEnthalpy_RT(hbar);
    hbar[0] *= RT();
}

void SingleSpeciesTP::getPartialMolarIntEnergies(double* ubar) const
{
    getIntEnergy_RT(ubar);
    ubar[0] *= RT();
}

void SingleSpeciesTP::getPartialMolarEntropies(double* sbar) const
{
    getEntropy_R(sbar);
    sbar[0] *= GasConstant;
}

void SingleSpeciesTP::getPartialMolarCp(double* cpbar) const
{
    getCp_R(cpbar);
    cpbar[0] *= GasConstant;
}

void SingleSpeciesTP::getPartialMolarVolumes(double* vbar) const
{
    vbar[0] = molecularWeight(0) / density();
}

// Properties of the Standard State of the Species in the Solution

void SingleSpeciesTP::getPureGibbs(double* gpure) const
{
    getGibbs_RT(gpure);
    gpure[0] *= RT();
}

void SingleSpeciesTP::getStandardVolumes(double* vbar) const
{
    vbar[0] = molecularWeight(0) / density();
}

// ---- Thermodynamic Values for the Species Reference States -------

void SingleSpeciesTP::getEnthalpy_RT_ref(double* hrt) const
{
    _updateThermo();
    hrt[0] = m_h0_RT;
}

void SingleSpeciesTP::getGibbs_RT_ref(double* grt) const
{
    _updateThermo();
    grt[0] = m_h0_RT - m_s0_R;
}

void SingleSpeciesTP::getGibbs_ref(double* g) const
{
    getGibbs_RT_ref(g);
    g[0] *= RT();
}

void SingleSpeciesTP::getEntropy_R_ref(double* er) const
{
    _updateThermo();
    er[0] = m_s0_R;
}

void SingleSpeciesTP::getCp_R_ref(double* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R;
}

bool SingleSpeciesTP::addSpecies(shared_ptr<Species> spec)
{
    if (m_kk != 0) {
        throw CanteraError("SingleSpeciesTP::addSpecies",
            "Stoichiometric substances may only contain one species.");
    }
    return ThermoPhase::addSpecies(spec);
}

void SingleSpeciesTP::_updateThermo() const
{
    double tnow = temperature();
    if (m_tlast != tnow) {
        m_spthermo.update(tnow, &m_cp0_R, &m_h0_RT, &m_s0_R);
        m_tlast = tnow;
    }
}

}
