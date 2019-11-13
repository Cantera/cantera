/**
 *  @file MaskellSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see \ref
 *      thermoprops and \link Cantera::MaskellSolidSolnPhase
 *      MaskellSolidSolnPhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/MaskellSolidSolnPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/xml.h"

#include <cassert>

namespace Cantera
{

MaskellSolidSolnPhase::MaskellSolidSolnPhase() :
    m_Pcurrent(OneAtm),
    h_mixing(0.0),
    product_species_index(-1),
    reactant_species_index(-1)
{
}

void MaskellSolidSolnPhase::getActivityConcentrations(doublereal* c) const
{
    getActivityCoefficients(c);
    for (size_t sp = 0; sp < m_kk; ++sp) {
        c[sp] *= moleFraction(sp);
    }
}

// Molar Thermodynamic Properties of the Solution

doublereal MaskellSolidSolnPhase::enthalpy_mole() const
{
    const doublereal h0 = RT() * mean_X(m_h0_RT);
    const doublereal r = moleFraction(product_species_index);
    const doublereal fmval = fm(r);
    return h0 + r * fmval * h_mixing;
}

doublereal xlogx(doublereal x)
{
    return x * std::log(x);
}

doublereal MaskellSolidSolnPhase::entropy_mole() const
{
    const doublereal s0 = GasConstant * mean_X(m_s0_R);
    const doublereal r = moleFraction(product_species_index);
    const doublereal fmval = fm(r);
    const doublereal rfm = r * fmval;
    return s0 + GasConstant * (xlogx(1-rfm) - xlogx(rfm) - xlogx(1-r-rfm) - xlogx((1-fmval)*r) - xlogx(1-r) - xlogx(r));
}

// Mechanical Equation of State

void MaskellSolidSolnPhase::setDensity(const doublereal rho)
{
    // Unless the input density is exactly equal to the density calculated and
    // stored in the State object, we throw an exception. This is because the
    // density is NOT an independent variable.
    warn_deprecated("MaskellSolidSolnPhase::setDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setDensity instead");
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
    doublereal vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFracs[i];
    }
    Phase::assignDensity(meanMolecularWeight() / vtotal);
}

void MaskellSolidSolnPhase::setPressure(doublereal p)
{
    m_Pcurrent = p;
}

void MaskellSolidSolnPhase::setMolarDensity(const doublereal n)
{
    warn_deprecated("MaskellSolidSolnPhase::setMolarDensity",
        "Overloaded function to be removed after Cantera 2.5. "
        "Error will be thrown by Phase::setMolarDensity instead");
    throw CanteraError("MaskellSolidSolnPhase::setMolarDensity",
                       "Density is not an independent variable");
}

// Chemical Potentials and Activities

void MaskellSolidSolnPhase::getActivityCoefficients(doublereal* ac) const
{
    static const int cacheId = m_cache.getId();
    CachedArray cached = m_cache.getArray(cacheId);
    if (!cached.validate(temperature(), pressure(), stateMFNumber())) {
        cached.value.resize(2);

        const doublereal r = moleFraction(product_species_index);
        const doublereal pval = p(r);
        const doublereal rfm = r * fm(r);
        const doublereal A = (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval)) /
                             (std::pow(1 - r - rfm, 1 + pval) * (1 - r));
        const doublereal B = pval * h_mixing / RT();
        cached.value[product_species_index] = A * std::exp(B);
        cached.value[reactant_species_index] = 1 / (A * r * (1-r) ) * std::exp(-B);
    }
    std::copy(cached.value.begin(), cached.value.end(), ac);
}

void MaskellSolidSolnPhase::getChemPotentials(doublereal* mu) const
{
    const doublereal r = moleFraction(product_species_index);
    const doublereal pval = p(r);
    const doublereal rfm = r * fm(r);
    const doublereal DgbarDr = pval * h_mixing +
                               RT() *
                               std::log( (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval) * r) /
                               (std::pow(1 - r - rfm, 1 + pval) * (1 - r)) );
    mu[product_species_index] = RT() * m_g0_RT[product_species_index] + DgbarDr;
    mu[reactant_species_index] = RT() * m_g0_RT[reactant_species_index] - DgbarDr;
}

void MaskellSolidSolnPhase::getChemPotentials_RT(doublereal* mu) const
{
    getChemPotentials(mu);
    for (size_t sp=0; sp < m_kk; ++sp) {
        mu[sp] *= 1.0 / RT();
    }
}

// Partial Molar Properties

void MaskellSolidSolnPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarEnthalpies");
}

void MaskellSolidSolnPhase::getPartialMolarEntropies(doublereal* sbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarEntropies");
}

void MaskellSolidSolnPhase::getPartialMolarCp(doublereal* cpbar) const
{
    throw NotImplementedError("MaskellSolidSolnPhase::getPartialMolarCp");
}

void MaskellSolidSolnPhase::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

void MaskellSolidSolnPhase::getPureGibbs(doublereal* gpure) const
{
    for (size_t sp=0; sp < m_kk; ++sp) {
        gpure[sp] = RT() * m_g0_RT[sp];
    }
}

void MaskellSolidSolnPhase::getStandardChemPotentials(doublereal* mu) const
{
    // What is the difference between this and getPureGibbs? IdealSolidSolnPhase
    // gives the same for both
    getPureGibbs(mu);
}

// Utility Functions

void MaskellSolidSolnPhase::initThermo()
{
    if (m_input.hasKey("excess-enthalpy")) {
        set_h_mix(m_input.convert("excess-enthalpy", "J/kmol"));
    }
    if (m_input.hasKey("product-species")) {
        setProductSpecies(m_input["product-species"].asString());
    }
    VPStandardStateTP::initThermo();
}


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

doublereal MaskellSolidSolnPhase::s() const
{
    return 1 + std::exp(h_mixing / RT());
}

doublereal MaskellSolidSolnPhase::fm(const doublereal r) const
{
    return (1 - std::sqrt(1 - 4*r*(1-r)/s())) / (2*r);
}

doublereal MaskellSolidSolnPhase::p(const doublereal r) const
{
    const doublereal sval = s();
    return (1 - 2*r) / std::sqrt(sval*sval - 4 * sval * r + 4 * sval * r * r);
}

} // end namespace Cantera
