/**
 *  @file GibbsExcessVPSSTP.cpp
 *   Definitions for intermediate ThermoPhase object for phases which
 *   employ excess gibbs free energy formulations
 *  (see \ref thermoprops
 * and class \link Cantera::GibbsExcessVPSSTP GibbsExcessVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon expressions
 * for the excess gibbs free energy expressed as a function of
 * the mole fractions.
 */
/*
 * Copyright (2009) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/thermo/GibbsExcessVPSSTP.h"
#include "cantera/base/stringUtils.h"

#include <iomanip>

using namespace std;

namespace Cantera
{

GibbsExcessVPSSTP::GibbsExcessVPSSTP() :
    VPStandardStateTP(),
    moleFractions_(0),
    lnActCoeff_Scaled_(0),
    dlnActCoeffdT_Scaled_(0),
    d2lnActCoeffdT2_Scaled_(0),
    dlnActCoeffdlnN_diag_(0),
    dlnActCoeffdlnX_diag_(0),
    dlnActCoeffdlnN_(0,0),
    m_pp(0)
{
}

GibbsExcessVPSSTP::GibbsExcessVPSSTP(const GibbsExcessVPSSTP& b) :
    VPStandardStateTP(),
    moleFractions_(0),
    lnActCoeff_Scaled_(0),
    dlnActCoeffdT_Scaled_(0),
    d2lnActCoeffdT2_Scaled_(0),
    dlnActCoeffdlnN_diag_(0),
    dlnActCoeffdlnX_diag_(0),
    dlnActCoeffdlnN_(0,0),
    m_pp(0)
{
    GibbsExcessVPSSTP::operator=(b);
}

GibbsExcessVPSSTP& GibbsExcessVPSSTP::
operator=(const GibbsExcessVPSSTP& b)
{
    if (&b == this) {
        return *this;
    }

    VPStandardStateTP::operator=(b);

    moleFractions_       = b.moleFractions_;
    lnActCoeff_Scaled_   = b.lnActCoeff_Scaled_;
    dlnActCoeffdT_Scaled_   = b.dlnActCoeffdT_Scaled_;
    d2lnActCoeffdT2_Scaled_   = b.d2lnActCoeffdT2_Scaled_;
    dlnActCoeffdlnX_diag_ = b.dlnActCoeffdlnX_diag_;
    dlnActCoeffdlnN_diag_ = b.dlnActCoeffdlnN_diag_;
    dlnActCoeffdlnN_  = b.dlnActCoeffdlnN_;
    m_pp                 = b.m_pp;

    return *this;
}

ThermoPhase*
GibbsExcessVPSSTP::duplMyselfAsThermoPhase() const
{
    return new GibbsExcessVPSSTP(*this);
}

void GibbsExcessVPSSTP::setMassFractions(const doublereal* const y)
{
    Phase::setMassFractions(y);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void GibbsExcessVPSSTP::setMassFractions_NoNorm(const doublereal* const y)
{
    Phase::setMassFractions_NoNorm(y);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void GibbsExcessVPSSTP::setMoleFractions(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void GibbsExcessVPSSTP::setMoleFractions_NoNorm(const doublereal* const x)
{
    Phase::setMoleFractions_NoNorm(x);
    getMoleFractions(DATA_PTR(moleFractions_));
}

void GibbsExcessVPSSTP::setConcentrations(const doublereal* const c)
{
    Phase::setConcentrations(c);
    getMoleFractions(DATA_PTR(moleFractions_));
}

int GibbsExcessVPSSTP::eosType() const
{
    return 0;
}

/*
 * ------------ Mechanical Properties ------------------------------
 */

void GibbsExcessVPSSTP::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
}

void GibbsExcessVPSSTP::calcDensity()
{
    vector_fp vbar = getPartialMolarVolumes();
    //    double *vbar = &m_pp[0];
    //    getPartialMolarVolumes(&vbar[0]);

    doublereal vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * moleFractions_[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;
    Phase::setDensity(dd);
}

void GibbsExcessVPSSTP::setState_TP(doublereal t, doublereal p)
{
    Phase::setTemperature(t);
    /*
     * Store the current pressure
     */
    m_Pcurrent = p;
    /*
     * update the standard state thermo
     * -> This involves calling the water function and setting the pressure
     */
    updateStandardStateThermo();

    /*
     * Calculate the partial molar volumes, and then the density of the fluid
     */
    calcDensity();
}

/*
 * - Activities, Standard States, Activity Concentrations -----------
 */

void GibbsExcessVPSSTP::getActivityConcentrations(doublereal* c) const
{
    getActivities(c);
}

doublereal GibbsExcessVPSSTP::standardConcentration(size_t k) const
{
    return 1.0;
}

doublereal GibbsExcessVPSSTP::logStandardConc(size_t k) const
{
    return 0.0;
}

void GibbsExcessVPSSTP::getActivities(doublereal* ac) const
{
    getActivityCoefficients(ac);
    getMoleFractions(DATA_PTR(moleFractions_));
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] *= moleFractions_[k];
    }
}

void GibbsExcessVPSSTP::getActivityCoefficients(doublereal* const ac) const
{

    getLnActivityCoefficients(ac);

    // Protect against roundoff when taking exponentials
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

void GibbsExcessVPSSTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
    double ve = Faraday * electricPotential();
    for (size_t k = 0; k < m_kk; k++) {
        mu[k] += ve*charge(k);
    }
}

/*
 * ------------ Partial Molar Properties of the Solution ------------
 */

void GibbsExcessVPSSTP::getPartialMolarVolumes(doublereal* vbar) const
{
    /*
     * Get the standard state values in m^3 kmol-1
     */
    getStandardVolumes(vbar);
}

const vector_fp& GibbsExcessVPSSTP::getPartialMolarVolumes() const
{
    return getStandardVolumes();
}

doublereal GibbsExcessVPSSTP::err(const std::string& msg) const
{
    throw CanteraError("GibbsExcessVPSSTP","Base class method "
                       +msg+" called. Equation of state type: "+int2str(eosType()));
    return 0;
}

double GibbsExcessVPSSTP::checkMFSum(const doublereal* const x) const
{
    doublereal norm = accumulate(x, x + m_kk, 0.0);
    if (fabs(norm - 1.0) > 1.0E-9) {
        throw CanteraError("GibbsExcessVPSSTP::checkMFSum",
                           "(MF sum - 1) exceeded tolerance of 1.0E-9:" + fp2str(norm));
    }
    return norm;
}

void GibbsExcessVPSSTP::getUnitsStandardConc(double* uA, int k, int sizeUA) const
{
    for (int i = 0; i < sizeUA; i++) {
        if (i == 0) {
            uA[0] = 0.0;
        }
        if (i == 1) {
            uA[1] = 0.0;
        }
        if (i == 2) {
            uA[2] = 0.0;
        }
        if (i == 3) {
            uA[3] = 0.0;
        }
        if (i == 4) {
            uA[4] = 0.0;
        }
        if (i == 5) {
            uA[5] = 0.0;
        }
    }
}

void GibbsExcessVPSSTP::initThermo()
{
    initLengths();
    VPStandardStateTP::initThermo();
    getMoleFractions(DATA_PTR(moleFractions_));
}

void  GibbsExcessVPSSTP::initLengths()
{
    m_kk = nSpecies();
    moleFractions_.resize(m_kk);
    lnActCoeff_Scaled_.resize(m_kk);
    dlnActCoeffdT_Scaled_.resize(m_kk);
    d2lnActCoeffdT2_Scaled_.resize(m_kk);
    dlnActCoeffdlnX_diag_.resize(m_kk);
    dlnActCoeffdlnN_diag_.resize(m_kk);
    dlnActCoeffdlnN_.resize(m_kk, m_kk);
    m_pp.resize(m_kk);
}

}
