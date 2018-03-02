/**
 *  @file MaskellSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see \ref
 *      thermoprops and \link Cantera::MaskellSolidSolnPhase
 *      MaskellSolidSolnPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/xml.h"

#include <cassert>

namespace Cantera
{

MaskellSolidSolnPhase::MaskellSolidSolnPhase() :
    m_Pcurrent(OneAtm),
    m_h0_RT(2),
    m_cp0_R(2),
    m_g0_RT(2),
    m_s0_R(2),
    h_mixing(0.0),
    product_species_index(-1),
    reactant_species_index(-1)
{
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
    _updateThermo();
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
    _updateThermo();
    const double s0 = GasConstant * mean_X(m_s0_R);
    const double r = moleFraction(product_species_index);
    const double fmval = fm(r);
    const double rfm = r * fmval;
    return s0 + GasConstant * (xlogx(1-rfm) - xlogx(rfm) - xlogx(1-r-rfm) - xlogx((1-fmval)*r) - xlogx(1-r) - xlogx(r));
}

// Mechanical Equation of State

void MaskellSolidSolnPhase::setDensity(const double rho)
{
    // Unless the input density is exactly equal to the density calculated and
    // stored in the State object, we throw an exception. This is because the
    // density is NOT an independent variable.
    double dens = density();
    if (rho != dens) {
        throw CanteraError("MaskellSolidSolnPhase::setDensity",
                           "Density is not an independent variable");
    }
}

void MaskellSolidSolnPhase::calcDensity()
{
    const vector_fp& vbar = getStandardVolumes();

    vector_fp moleFracs(m_kk);
    Phase::getMoleFractions(&moleFracs[0]);
    double vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFracs[i];
    }
    Phase::setDensity(meanMolecularWeight() / vtotal);
}

void MaskellSolidSolnPhase::setPressure(double p)
{
    m_Pcurrent = p;
}

void MaskellSolidSolnPhase::setMolarDensity(const double n)
{
    throw CanteraError("MaskellSolidSolnPhase::setMolarDensity",
                       "Density is not an independent variable");
}

// Chemical Potentials and Activities

void MaskellSolidSolnPhase::getActivityCoefficients(double* ac) const
{
    _updateThermo();
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
    _updateThermo();
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
    getChemPotentials(mu);
    for (size_t sp=0; sp < m_kk; ++sp) {
        mu[sp] *= 1.0 / RT();
    }
}

// Partial Molar Properties

void MaskellSolidSolnPhase::getPartialMolarEnthalpies(double* hbar) const
{
    throw CanteraError("MaskellSolidSolnPhase::getPartialMolarEnthalpies()", "Not yet implemented.");
}

void MaskellSolidSolnPhase::getPartialMolarEntropies(double* sbar) const
{
    throw CanteraError("MaskellSolidSolnPhase::getPartialMolarEntropies()", "Not yet implemented.");
}

void MaskellSolidSolnPhase::getPartialMolarCp(double* cpbar) const
{
    throw CanteraError("MaskellSolidSolnPhase::getPartialMolarCp()", "Not yet implemented.");
}

void MaskellSolidSolnPhase::getPartialMolarVolumes(double* vbar) const
{
    getStandardVolumes(vbar);
}

void MaskellSolidSolnPhase::getPureGibbs(double* gpure) const
{
    _updateThermo();
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

void MaskellSolidSolnPhase::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    if (id_.size() > 0 && phaseNode.id() != id_) {
        throw CanteraError("MaskellSolidSolnPhase::initThermoXML",
                           "phasenode and Id are incompatible");
    }

    // Check on the thermo field. Must have:
    // <thermo model="MaskellSolidSolution" />
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thNode = phaseNode.child("thermo");
        if (!caseInsensitiveEquals(thNode["model"], "maskellsolidsolnphase")) {
            throw CanteraError("MaskellSolidSolnPhase::initThermoXML",
                               "Unknown thermo model: " + thNode["model"]);
        }

        // Parse the enthalpy of mixing constant
        if (thNode.hasChild("h_mix")) {
            set_h_mix(fpValue(thNode.child("h_mix").value()));
        } else {
            throw CanteraError("MaskellSolidSolnPhase::initThermoXML",
                               "Mixing enthalpy parameter not specified.");
        }

        if (thNode.hasChild("product_species")) {
            setProductSpecies(thNode.child("product_species").value());
        } else {
            setProductSpecies(speciesName(0)); // default
        }
    } else {
        throw CanteraError("MaskellSolidSolnPhase::initThermoXML",
                           "Unspecified thermo model");
    }

    // Confirm that the phase only contains 2 species
    if (m_kk != 2) {
        throw CanteraError("MaskellSolidSolnPhase::initThermoXML",
                "MaskellSolidSolution model requires exactly 2 species.");
    }

    // Call the base initThermo, which handles setting the initial state.
    VPStandardStateTP::initThermoXML(phaseNode, id_);
}

void MaskellSolidSolnPhase::setProductSpecies(const std::string& name)
{
    product_species_index = static_cast<int>(speciesIndex(name));
    if (product_species_index == -1) {
        throw CanteraError("MaskellSolidSolnPhase::setProductSpecies",
                           "Species '{}' not found", name);
    }
    reactant_species_index = (product_species_index == 0) ? 1 : 0;
}

void MaskellSolidSolnPhase::_updateThermo() const
{
    assert(m_kk == 2);
    static const int cacheId = m_cache.getId();
    CachedScalar cached = m_cache.getScalar(cacheId);

    // Update the thermodynamic functions of the reference state.
    double tnow = temperature();
    if (!cached.validate(tnow)) {
        m_spthermo.update(tnow, m_cp0_R.data(), m_h0_RT.data(), m_s0_R.data());
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
    }
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
