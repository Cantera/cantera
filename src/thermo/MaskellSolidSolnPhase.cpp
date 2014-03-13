/**
 *  @file MaskellSolidSolnPhase.cpp Implementation file for an ideal solid
 *      solution model with incompressible thermodynamics (see \ref
 *      thermoprops and \link Cantera::MaskellSolidSolnPhase
 *      MaskellSolidSolnPhase\endlink).
 */
/*
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 */

#include "cantera/thermo/MaskellSolidSolnPhase.h"

#include "cantera/base/stringUtils.h"

#include <cassert>
#include <iostream>

namespace Cantera
{
//=====================================================================================================
MaskellSolidSolnPhase::MaskellSolidSolnPhase() :
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_updateThermo_valid(false),
    m_activity_coeffs(2),
    m_activity_coeffs_valid(false),
    m_h0_RT(2),
    m_cp0_R(2),
    m_g0_RT(2),
    m_s0_R(2),
    h_mixing(0.0),
    product_species_index(0),
    reactant_species_index(1)
{
}
//=====================================================================================================
MaskellSolidSolnPhase::MaskellSolidSolnPhase(const MaskellSolidSolnPhase& b) :
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_updateThermo_valid(false),
    m_activity_coeffs(2),
    m_activity_coeffs_valid(false),
    m_h0_RT(2),
    m_cp0_R(2),
    m_g0_RT(2),
    m_s0_R(2),
    h_mixing(0.0),
    product_species_index(0),
    reactant_species_index(1)
{
    *this = b;
}
//=====================================================================================================
MaskellSolidSolnPhase& MaskellSolidSolnPhase::
operator=(const MaskellSolidSolnPhase& b)
{
    if (this != &b) {
        VPStandardStateTP::operator=(b);
    }
    return *this;
}
//=====================================================================================================
ThermoPhase* MaskellSolidSolnPhase::duplMyselfAsThermoPhase() const
{
    return new MaskellSolidSolnPhase(*this);
}
//=====================================================================================================
void MaskellSolidSolnPhase::
getActivityConcentrations(doublereal* c) const
{
    getActivityCoefficients(c);
    for(unsigned sp=0; sp < m_kk; ++sp)
    {
      c[sp] *= moleFraction(sp);
    }
}

/********************************************************************
 *            Molar Thermodynamic Properties of the Solution
 ********************************************************************/
//=====================================================================================================
doublereal MaskellSolidSolnPhase::
enthalpy_mole() const
{
    _updateThermo();
    const doublereal h0 = GasConstant * temperature() * mean_X(&m_h0_RT[0]);
    const doublereal r = moleFraction(product_species_index);
    const doublereal fmval = fm(r);
    return h0 + r * fmval * h_mixing;
}
//=====================================================================================================
doublereal xlogx(doublereal x)
{
    return x * std::log(x);
}

doublereal MaskellSolidSolnPhase::entropy_mole() const
{
    _updateThermo();
    const doublereal s0 = GasConstant * mean_X(&m_s0_R[0]);
    const doublereal r = moleFraction(product_species_index);
    const doublereal fmval = fm(r);
    const doublereal rfm = r * fmval;
    return s0 + GasConstant * (xlogx(1-rfm) - xlogx(rfm) - xlogx(1-r-rfm) - xlogx((1-fmval)*r) - xlogx(1-r) - xlogx(r));
}

/********************************************************************
 *                  Mechanical Equation of State
 ********************************************************************/

void MaskellSolidSolnPhase::
setDensity(const doublereal rho)
{
    /*
     * Unless the input density is exactly equal to the density
     * calculated and stored in the State object, we throw an
     * exception. This is because the density is NOT an
     * independent variable.
     */
    double dens = density();
    if (rho != dens) {
        throw CanteraError("MaskellSolidSolnPhase::setDensity",
                           "Density is not an independent variable");
    }
}

void MaskellSolidSolnPhase::
calcDensity()
{
    const vector_fp & vbar = getStandardVolumes();

    vector_fp moleFracs(m_kk);
    Phase::getMoleFractions(&moleFracs[0]);
    doublereal vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFracs[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;
    Phase::setDensity(dd);
}

void MaskellSolidSolnPhase::setPressure(doublereal p)
{
    m_Pcurrent = p;
}

void MaskellSolidSolnPhase::setMolarDensity(const doublereal n)
{
    throw CanteraError("MaskellSolidSolnPhase::setMolarDensity",
                       "Density is not an independent variable");
}

/********************************************************************
 *        Chemical Potentials and Activities
 ********************************************************************/

void MaskellSolidSolnPhase::
getActivityCoefficients(doublereal* ac) const
{
    _updateThermo();
    if( !m_activity_coeffs_valid ) {
      m_activity_coeffs_valid = true;
      const doublereal r = moleFraction(product_species_index);
      const doublereal pval = p(r);
      const doublereal fmval = fm(r);
      const doublereal rfm = r * fmval;
      const doublereal RT = GasConstant * temperature();
      const doublereal A = (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval))  /
                           (std::pow(1 - r - rfm, 1 + pval) * (1 - r));
      const doublereal B = pval * h_mixing / RT;
      m_activity_coeffs[product_species_index] = A * std::exp(B);
      m_activity_coeffs[reactant_species_index] = 1 / (A * r * (1-r) ) * std::exp(-B);
    }
    std::copy(m_activity_coeffs.begin(), m_activity_coeffs.end(), ac);
}

void MaskellSolidSolnPhase::
getChemPotentials(doublereal* mu) const
{
    _updateThermo();
    const doublereal r = moleFraction(product_species_index);
    const doublereal pval = p(r);
    const doublereal fmval = fm(r);
    const doublereal rfm = r * fmval;
    const doublereal RT = GasConstant * temperature();
    const doublereal DgbarDr = pval * h_mixing +
                               GasConstant * temperature() *
                               std::log( (std::pow(1 - rfm, pval) * std::pow(rfm, pval) * std::pow(r - rfm, 1 - pval) * r)  /
                               (std::pow(1 - r - rfm, 1 + pval) * (1 - r)) );
    mu[product_species_index] = RT * m_g0_RT[product_species_index] + DgbarDr;
    mu[reactant_species_index] = RT * m_g0_RT[reactant_species_index] - DgbarDr;
}

void MaskellSolidSolnPhase::
getChemPotentials_RT(doublereal* mu) const
{
  const doublereal invRT = 1.0 / (GasConstant * temperature());
  getChemPotentials(mu);
  for(unsigned sp=0; sp < m_kk; ++sp)
  {
    mu[sp] *= invRT;
  }
}

/********************************************************************
 *                    Partial Molar Properties
 ********************************************************************/

void MaskellSolidSolnPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
  throw CanteraError("MaskellSolidSolnPhase::getPartialMolarEnthalpies()", "Not yet implemented.");
}

void MaskellSolidSolnPhase::
getPartialMolarEntropies(doublereal* sbar) const
{
  throw CanteraError("MaskellSolidSolnPhase::getPartialMolarEntropies()", "Not yet implemented.");
}

void MaskellSolidSolnPhase::
getPartialMolarCp(doublereal* cpbar) const
{
  throw CanteraError("MaskellSolidSolnPhase::getPartialMolarCp()", "Not yet implemented.");
}

void MaskellSolidSolnPhase::
getPartialMolarVolumes(doublereal* vbar) const
{
  getStandardVolumes(vbar);
}

void MaskellSolidSolnPhase::
getPureGibbs(doublereal* gpure) const
{
    _updateThermo();
    const doublereal RT = GasConstant * temperature();
    for(unsigned sp=0; sp < m_kk; ++sp)
    {
        gpure[sp] = RT * m_g0_RT[sp];
    }
}

void MaskellSolidSolnPhase::
getStandardChemPotentials(doublereal* mu) const
{
  // What is the difference between this and getPureGibbs? IdealSolidSolnPhase gives the same for both
  getPureGibbs(mu);
}

/*********************************************************************
 *    Utility Functions
 *********************************************************************/
void MaskellSolidSolnPhase::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
    std::string subname = "MaskellSolidSolnPhase::initThermoXML";
    if (id_.size() > 0) {
        std::string idp = phaseNode.id();
        if (idp != id_) {
            throw CanteraError(subname.c_str(),
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Check on the thermo field. Must have:
     * <thermo model="MaskellSolidSolution" />
     */
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thNode = phaseNode.child("thermo");
        std::string mStringa = thNode.attrib("model");
        std::string mString = lowercase(mStringa);
        if (mString != "maskellsolidsolnphase") {
            throw CanteraError(subname.c_str(),
                               "Unknown thermo model: " + mStringa);
        }


        /*
         * Parse the enthalpy of mixing constant
         */
        if (thNode.hasChild("h_mix")) {
            XML_Node& scNode = thNode.child("h_mix");
            set_h_mix(fpValue(scNode.value()));
        } else {
            throw CanteraError(subname.c_str(),
                               "Mixing enthalpy parameter not specified.");
        }

        if (thNode.hasChild("product_species")) {
            XML_Node& scNode = thNode.child("product_species");
            std::string product_species_name = scNode.value();
            product_species_index = speciesIndex(product_species_name);
            if( product_species_index == static_cast<int>(npos) )
            {
              throw CanteraError(subname.c_str(),
                                 "Species " + product_species_name + " not found.");
            }
            if( product_species_index == 0 )
            {
              reactant_species_index = 1;
            } else {
              reactant_species_index = 0;
            }
        }

    } else {
        throw CanteraError(subname.c_str(),
                           "Unspecified thermo model");
    }


    // Confirm that the phase only contains 2 species
    if( m_kk != 2 )
    {
      throw CanteraError( subname.c_str(), "MaskellSolidSolution model requires exactly 2 species.");
    }

    /*
     * Call the base initThermo, which handles setting the initial
     * state.
     */
    VPStandardStateTP::initThermoXML(phaseNode, id_);
}

void MaskellSolidSolnPhase::_updateThermo() const
{
    assert(m_kk == 2);
    if (!m_updateThermo_valid) {
        m_updateThermo_valid = true;
        /*
         * Update the thermodynamic functions of the reference state.
         */
        doublereal tnow = temperature();
        m_spthermo->update(tnow, DATA_PTR(m_cp0_R), DATA_PTR(m_h0_RT),
                           DATA_PTR(m_s0_R));
        m_tlast = tnow;
        for (size_t k = 0; k < m_kk; k++) {
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}

doublereal MaskellSolidSolnPhase::s() const
{
  return 1 + std::exp(h_mixing / (GasConstant * temperature()));
}

doublereal MaskellSolidSolnPhase::fm(const doublereal r) const
{
  const doublereal sval = s();
  return (1 - std::sqrt(1 - 4*r*(1-r)/sval)) / (2*r);
}

doublereal MaskellSolidSolnPhase::p(const doublereal r) const
{
  const doublereal sval = s();
  return (1 - 2*r) / std::sqrt(sval*sval - 4 * sval * r + 4 * sval * r * r);
}

void MaskellSolidSolnPhase::invalidateCachedDataOnStateChange(StateVariable changed_var)
{
  m_updateThermo_valid = false;
  m_activity_coeffs_valid = false;
  VPStandardStateTP::invalidateCachedDataOnStateChange(changed_var);
}

}  // end namespace Cantera
