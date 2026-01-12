/**
 *  @file GibbsExcessVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess Gibbs free energy formulations
 *  (see @ref thermoprops and class @link Cantera::GibbsExcessVPSSTP GibbsExcessVPSSTP@endlink).
 *
 * Header file for a derived class of ThermoPhase that handles variable pressure
 * standard state methods for calculating thermodynamic properties that are
 * further based upon expressions for the excess Gibbs free energy expressed as
 * a function of the mole fractions.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/GibbsExcessVPSSTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"
#include <numeric>

namespace Cantera
{

void GibbsExcessVPSSTP::compositionChanged()
{
    Phase::compositionChanged();
    getMoleFractions(moleFractions_);
}

// ------------ Mechanical Properties ------------------------------

void GibbsExcessVPSSTP::calcDensity()
{
    const vector<double>& vbar = getStandardVolumes();
    double vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFractions_[i];
    }
    double dd = meanMolecularWeight() / vtotal;
    Phase::assignDensity(dd);
}

// - Activities, Standard States, Activity Concentrations -----------

Units GibbsExcessVPSSTP::standardConcentrationUnits() const
{
    return Units(1.0); // dimensionless
}

void GibbsExcessVPSSTP::getActivityConcentrations(double* c) const
{
    getActivities(c);
}

double GibbsExcessVPSSTP::standardConcentration(size_t k) const
{
    return 1.0;
}

double GibbsExcessVPSSTP::logStandardConc(size_t k) const
{
    return 0.0;
}

void GibbsExcessVPSSTP::getActivities(double* ac) const
{
    getActivityCoefficients(ac);
    getMoleFractions(moleFractions_);
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] *= moleFractions_[k];
    }
}

void GibbsExcessVPSSTP::getActivityCoefficients(double* const ac) const
{
    getLnActivityCoefficients(ac);
    for (size_t k = 0; k < m_kk; k++) {
        if (ac[k] > 700.) {
            ac[k] = exp(700.0);
        } else if (ac[k] < -700.) {
            ac[k] = exp(-700.0);
        } else {
            ac[k] = exp(ac[k]);
        }
    }
}

// ------------ Partial Molar Properties of the Solution ------------

void GibbsExcessVPSSTP::getPartialMolarVolumes(double* vbar) const
{
    // Get the standard state values in m^3 kmol-1
    getStandardVolumes(vbar);
}

bool GibbsExcessVPSSTP::addSpecies(shared_ptr<Species> spec)
{
    bool added = VPStandardStateTP::addSpecies(spec);
    if (added) {
        if (m_kk == 1) {
            moleFractions_.push_back(1.0);
        } else {
            moleFractions_.push_back(0.0);
        }
        lnActCoeff_Scaled_.push_back(0.0);
        dlnActCoeffdT_Scaled_.push_back(0.0);
        d2lnActCoeffdT2_Scaled_.push_back(0.0);
        dlnActCoeffdlnX_diag_.push_back(0.0);
        dlnActCoeffdlnN_diag_.push_back(0.0);
        dlnActCoeffdlnN_.resize(m_kk, m_kk);
    }
    return added;
}

} // end of namespace Cantera
