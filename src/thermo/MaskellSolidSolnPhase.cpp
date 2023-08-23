/**
 *  @file MaskellSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see @ref
 *      thermoprops and @link Cantera::MaskellSolidSolnPhase
 *      MaskellSolidSolnPhase@endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

#include <cassert>

namespace Cantera
{

MaskellSolidSolnPhase::MaskellSolidSolnPhase()
{
    warn_deprecated("class MaskellSolidSolnPhase", "To be removed after Cantera 3.0");
}

void MaskellSolidSolnPhase::getActivityConcentrations(double* c) const
{
    getActivityCoefficients(c);
    for (size_t sp = 0; sp < m_kk; ++sp) {
        c[sp] *= moleFraction(sp);
    }
}

// Molar Thermodynamic Properties of the Solution

double MaskellSolidSolnPhase::enthalpy_mole() const
{
    const double h0 = RT() * mean_X(m_h0_RT);
    const double r = moleFraction(product_species_index);
    const double fmval = fm(r);
    return h0 + r * fmval * h_mixing;
}

double xlogx(double x)
{
    return x * std::log(x);
}

double MaskellSolidSolnPhase::entropy_mole() const
{
    const double s0 = GasConstant * mean_X(m_s0_R);
    const double r = moleFraction(product_species_index);
    const double fmval = fm(r);
    const double rfm = r * fmval;
    return s0 + GasConstant * (xlogx(1-rfm) - xlogx(rfm) - xlogx(1-r-rfm) - xlogx((1-fmval)*r) - xlogx(1-r) - xlogx(r));
}

// Mechanical Equation of State

void MaskellSolidSolnPhase::calcDensity()
{
    const vector<double>& vbar = getStandardVolumes();

    vector<double> moleFracs(m_kk);
    Phase::getMoleFractions(&moleFracs[0]);
    double vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFracs[i];
    }
    Phase::assignDensity(meanMolecularWeight() / vtotal);
}

void MaskellSolidSolnPhase::setPressure(double p)
{
    m_Pcurrent = p;
}

// Chemical Potentials and Activities

void MaskellSolidSolnPhase::getActivityCoefficients(double* ac) const
{
    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);
    if (!cached.validate(temperature(), pressure(), stateMFNumber())) {
        cached.value.resize(2);

        const double r = moleFraction(product_species_index);
        const double pval = p(r);
        const double rfm = r * fm(r);
        const double A = (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval)) /
                             (std::pow(1 - r - rfm, 1 + pval) * (1 - r));
        const double B = pval * h_mixing / RT();
        cached.value[product_species_index] = A * std::exp(B);
        cached.value[reactant_species_index] = 1 / (A * r * (1-r) ) * std::exp(-B);
    }
    std::copy(cached.value.begin(), cached.value.end(), ac);
}

void MaskellSolidSolnPhase::getChemPotentials(double* mu) const
{
    const double r = moleFraction(product_species_index);
    const double pval = p(r);
    const double rfm = r * fm(r);
    const double DgbarDr = pval * h_mixing +
                               RT() *
                               std::log( (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval) * r) /
                               (std::pow(1 - r - rfm, 1 + pval) * (1 - r)) );
    mu[product_species_index] = RT() * m_g0_RT[product_species_index] + DgbarDr;
    mu[reactant_species_index] = RT() * m_g0_RT[reactant_species_index] - DgbarDr;
}

void MaskellSolidSolnPhase::getChemPotentials_RT(double* mu) const
{
    warn_deprecated("MaskellSolidSolnPhase::getChemPotentials_RT",
                    "To be removed after Cantera 3.0. Use getChemPotentials instead.");
    getChemPotentials(mu);
    for (size_t sp=0; sp < m_kk; ++sp) {
        mu[sp] *= 1.0 / RT();
    }
}

// Partial Molar Properties

void MaskellSolidSolnPhase::getPartialMolarEnthalpies(double* hbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarEnthalpies");
}

void MaskellSolidSolnPhase::getPartialMolarEntropies(double* sbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarEntropies");
}

void MaskellSolidSolnPhase::getPartialMolarCp(double* cpbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarCp");
}

void MaskellSolidSolnPhase::getPartialMolarVolumes(double* vbar) const
{
    getStandardVolumes(vbar);
}

void MaskellSolidSolnPhase::getPureGibbs(double* gpure) const
{
    for (size_t sp=0; sp < m_kk; ++sp) {
        gpure[sp] = RT() * m_g0_RT[sp];
    }
}

void MaskellSolidSolnPhase::getStandardChemPotentials(double* mu) const
{
    // What is the difference between this and getPureGibbs? IdealSolidSolnPhase
    // gives the same for both
    getPureGibbs(mu);
}

// Utility Functions

void MaskellSolidSolnPhase::initThermo()
{
    if (!m_input.empty()) {
        set_h_mix(m_input.convert("excess-enthalpy", "J/kmol"));
        setProductSpecies(m_input["product-species"].asString());
    }
    VPStandardStateTP::initThermo();
}

void MaskellSolidSolnPhase::getParameters(AnyMap& phaseNode) const
{
    VPStandardStateTP::getParameters(phaseNode);
    phaseNode["excess-enthalpy"].setQuantity(h_mixing, "J/kmol");
    phaseNode["product-species"] = speciesName(product_species_index);
}

void MaskellSolidSolnPhase::setProductSpecies(const string& name)
{
    product_species_index = static_cast<int>(speciesIndex(name));
    if (product_species_index == -1) {
        throw CanteraError("MaskellSolidSolnPhase::setProductSpecies",
                           "Species '{}' not found", name);
    }
    reactant_species_index = (product_species_index == 0) ? 1 : 0;
}

double MaskellSolidSolnPhase::s() const
{
    return 1 + std::exp(h_mixing / RT());
}

double MaskellSolidSolnPhase::fm(const double r) const
{
    return (1 - std::sqrt(1 - 4*r*(1-r)/s())) / (2*r);
}

double MaskellSolidSolnPhase::p(const double r) const
{
    const double sval = s();
    return (1 - 2*r) / std::sqrt(sval*sval - 4 * sval * r + 4 * sval * r * r);
}

} // end namespace Cantera
