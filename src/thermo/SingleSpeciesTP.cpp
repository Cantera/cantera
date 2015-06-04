/**
 *  @file SingleSpeciesTP.cpp
 *  Definitions for the SingleSpeciesTP class, which is a filter class for ThermoPhase,
 *  that eases the construction of single species phases
 *  ( see \ref thermoprops and class \link Cantera::SingleSpeciesTP SingleSpeciesTP\endlink).
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/thermo/SingleSpeciesTP.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{
SingleSpeciesTP::SingleSpeciesTP() :
    m_press(OneAtm),
    m_p0(OneAtm)
{
}

SingleSpeciesTP::SingleSpeciesTP(const SingleSpeciesTP& right):
    m_press(OneAtm),
    m_p0(OneAtm)
{
    *this = right;
}

SingleSpeciesTP& SingleSpeciesTP::operator=(const SingleSpeciesTP& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_press      = right.m_press;
        m_p0         = right.m_p0;
        m_h0_RT      = right.m_h0_RT;
        m_cp0_R      = right.m_cp0_R;
        m_s0_R       = right.m_s0_R;
    }
    return *this;
}

ThermoPhase* SingleSpeciesTP::duplMyselfAsThermoPhase() const
{
    return new SingleSpeciesTP(*this);
}

int SingleSpeciesTP::eosType() const
{
    throw NotImplementedError("SingleSpeciesTP::eosType");
}

/*
 * ------------ Molar Thermodynamic Properties --------------------
 */

doublereal SingleSpeciesTP::enthalpy_mole() const
{
    double hbar;
    getPartialMolarEnthalpies(&hbar);
    return hbar;
}

doublereal SingleSpeciesTP::intEnergy_mole() const
{
    double ubar;
    getPartialMolarIntEnergies(&ubar);
    return ubar;
}

doublereal SingleSpeciesTP::entropy_mole() const
{
    double sbar;
    getPartialMolarEntropies(&sbar);
    return sbar;
}

doublereal SingleSpeciesTP::gibbs_mole() const
{
    double gbar;
    /*
     * Get the chemical potential of the first species.
     * This is the same as the partial molar Gibbs
     * free energy.
     */
    getChemPotentials(&gbar);
    return gbar;
}

doublereal SingleSpeciesTP::cp_mole() const
{
    double cpbar;
    /*
     * Really should have a partial molar heat capacity
     * function in ThermoPhase. However, the standard
     * state heat capacity will do fine here for now.
     */
    getCp_R(&cpbar);
    cpbar *= GasConstant;
    return cpbar;
}

doublereal SingleSpeciesTP::cv_mole() const
{
    /*
     *  For single species, we go directory to the general Cp - Cv relation
     *
     *  Cp = Cv + alpha**2 * V * T / beta
     *
     * where
     *     alpha = volume thermal expansion coefficient
     *     beta  = isothermal compressibility
     */
    doublereal cvbar = cp_mole();
    doublereal alpha = thermalExpansionCoeff();
    doublereal beta = isothermalCompressibility();
    doublereal V = molecularWeight(0)/density();
    doublereal T = temperature();
    if (beta != 0.0) {
        cvbar -= alpha * alpha * V * T / beta;
    }
    return cvbar;
}

/*
 * ----------- Partial Molar Properties of the Solution -----------------
 */

void SingleSpeciesTP::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
}

void SingleSpeciesTP::getChemPotentials_RT(doublereal* murt) const
{
    getStandardChemPotentials(murt);
    murt[0] /= GasConstant * temperature();
}

void SingleSpeciesTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
}

void SingleSpeciesTP::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    hbar[0] *= GasConstant * temperature();
}

void SingleSpeciesTP::getPartialMolarIntEnergies(doublereal* ubar) const
{
    getIntEnergy_RT(ubar);
    ubar[0] *= GasConstant * temperature();
}

void SingleSpeciesTP::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    sbar[0] *= GasConstant;
}

void SingleSpeciesTP::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    cpbar[0] *= GasConstant;
}

void SingleSpeciesTP::getPartialMolarVolumes(doublereal* vbar) const
{
    vbar[0] = molecularWeight(0) / density();
}

/*
 * Properties of the Standard State of the Species in the Solution
 */

void SingleSpeciesTP::getPureGibbs(doublereal* gpure) const
{
    getGibbs_RT(gpure);
    gpure[0] *= GasConstant * temperature();
}

void SingleSpeciesTP::getStandardVolumes(doublereal* vbar) const
{
    vbar[0] = molecularWeight(0) / density();
}

/*
 * ---- Thermodynamic Values for the Species Reference States -------
 */

void SingleSpeciesTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    _updateThermo();
    hrt[0] = m_h0_RT[0];
}

void SingleSpeciesTP::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    grt[0] = m_h0_RT[0] - m_s0_R[0];
}

void SingleSpeciesTP::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    g[0] *= GasConstant * temperature();
}

void SingleSpeciesTP::getEntropy_R_ref(doublereal* er) const
{
    _updateThermo();
    er[0] = m_s0_R[0];
}

void SingleSpeciesTP::getCp_R_ref(doublereal* cpr) const
{
    _updateThermo();
    cpr[0] = m_cp0_R[0];
}

/*
 * ------------------ Setting the State ------------------------
 */

void SingleSpeciesTP::setState_HP(doublereal h, doublereal p,
                                  doublereal tol)
{
    doublereal dt;
    setPressure(p);
    for (int n = 0; n < 50; n++) {
        dt = clip((h - enthalpy_mass())/cp_mass(), -100.0, 100.0);
        setState_TP(temperature() + dt, p);
        if (fabs(dt) < tol) {
            return;
        }
    }
    throw CanteraError("setState_HP","no convergence. dt = " + fp2str(dt));
}

void SingleSpeciesTP::setState_UV(doublereal u, doublereal v,
                                  doublereal tol)
{
    doublereal dt;
    if (v == 0.0) {
        setDensity(1.0E100);
    } else {
        setDensity(1.0/v);
    }
    for (int n = 0; n < 50; n++) {
        dt = clip((u - intEnergy_mass())/cv_mass(), -100.0, 100.0);
        setTemperature(temperature() + dt);
        if (fabs(dt) < tol) {
            return;
        }
    }
    throw CanteraError("setState_UV",
                       "no convergence. dt = " + fp2str(dt)+"\n"
                       +"u = "+fp2str(u)+" v = "+fp2str(v)+"\n");
}

void SingleSpeciesTP::setState_SP(doublereal s, doublereal p,
                                  doublereal tol)
{
    doublereal dt;
    setPressure(p);
    for (int n = 0; n < 50; n++) {
        dt = clip((s - entropy_mass())*temperature()/cp_mass(), -100.0, 100.0);
        setState_TP(temperature() + dt, p);
        if (fabs(dt) < tol) {
            return;
        }
    }
    throw CanteraError("setState_SP","no convergence. dt = " + fp2str(dt));
}

void SingleSpeciesTP::setState_SV(doublereal s, doublereal v,
                                  doublereal tol)
{
    doublereal dt;
    if (v == 0.0) {
        setDensity(1.0E100);
    } else {
        setDensity(1.0/v);
    }
    for (int n = 0; n < 50; n++) {
        dt = clip((s - entropy_mass())*temperature()/cv_mass(), -100.0, 100.0);
        setTemperature(temperature() + dt);
        if (fabs(dt) < tol) {
            return;
        }
    }
    throw CanteraError("setState_SV","no convergence. dt = " + fp2str(dt));
}

void SingleSpeciesTP::initThermo()
{
    /*
     * Make sure there is one and only one species in this phase.
     */
    if (nSpecies() != 1) {
        throw CanteraError("initThermo",
                           "stoichiometric substances may only contain one species.");
    }

    /*
     * Resize temporary arrays.
     */
    m_h0_RT.resize(1);
    m_cp0_R.resize(1);
    m_s0_R.resize(1);

    /*
     *  Make sure the species mole fraction is equal to 1.0;
     */
    double x = 1.0;
    ThermoPhase::setMoleFractions(&x);
    /*
     * Call the base class initThermo object.
     */
    ThermoPhase::initThermo();
}

void SingleSpeciesTP::_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        m_spthermo->update(tnow, DATA_PTR(m_cp0_R), DATA_PTR(m_h0_RT),
                           DATA_PTR(m_s0_R));
        m_tlast = tnow;
    }
}

}
