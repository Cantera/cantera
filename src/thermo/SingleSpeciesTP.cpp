/**
 *  @file SingleSpeciesTP.cpp
 *  Definitions for the %SingleSpeciesTP class, which is a filter class for %ThermoPhase,
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

using namespace std;

namespace Cantera
{
SingleSpeciesTP::SingleSpeciesTP() :
    ThermoPhase(),
    m_press(OneAtm),
    m_p0(OneAtm),
    m_tlast(-1.0)
{
}

SingleSpeciesTP::SingleSpeciesTP(const SingleSpeciesTP& right):
    ThermoPhase(),
    m_press(OneAtm),
    m_p0(OneAtm),
    m_tlast(-1.0)
{
    *this = operator=(right);
}

SingleSpeciesTP& SingleSpeciesTP::operator=(const SingleSpeciesTP& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_press      = right.m_press;
        m_p0         = right.m_p0;
        m_tlast      = right.m_tlast;
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
    err("eosType");
    return -1;
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
    //getPartialMolarCp(&cpbar);
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
    doublereal molecW = molecularWeight(0);
    doublereal V = molecW/density();
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
    double rt = GasConstant * temperature();
    murt[0] /= rt;
}

void SingleSpeciesTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
}

void SingleSpeciesTP::
getPartialMolarEnthalpies(doublereal* hbar) const
{
    double _rt = GasConstant * temperature();
    getEnthalpy_RT(hbar);
    hbar[0] *= _rt;
}

void SingleSpeciesTP::
getPartialMolarIntEnergies(doublereal* ubar) const
{
    double _rt = GasConstant * temperature();
    getIntEnergy_RT(ubar);
    ubar[0] *= _rt;
}

void SingleSpeciesTP::
getPartialMolarEntropies(doublereal* sbar) const
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
    double mw = molecularWeight(0);
    double dens = density();
    vbar[0] = mw / dens;
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
    double mw = molecularWeight(0);
    double dens = density();
    vbar[0] = mw / dens;
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

void SingleSpeciesTP::setState_TPX(doublereal t, doublereal p,
                                   const doublereal* x)
{
    setTemperature(t);
    setPressure(p);
}

void SingleSpeciesTP::setState_TPX(doublereal t, doublereal p,
                                   compositionMap& x)
{
    setTemperature(t);
    setPressure(p);
}

void SingleSpeciesTP::setState_TPX(doublereal t, doublereal p,
                                   const std::string& x)
{
    setTemperature(t);
    setPressure(p);
}

void SingleSpeciesTP::setState_TPY(doublereal t, doublereal p,
                                   const doublereal* y)
{
    setTemperature(t);
    setPressure(p);
}

void SingleSpeciesTP::setState_TPY(doublereal t, doublereal p,
                                   compositionMap& y)
{
    setTemperature(t);
    setPressure(p);
}

void SingleSpeciesTP::setState_TPY(doublereal t, doublereal p,
                                   const std::string& y)
{
    setTemperature(t);
    setPressure(p);
}

void SingleSpeciesTP::setState_PX(doublereal p, doublereal* x)
{
    if (x[0] != 1.0) {
        err("setStatePX -> x[0] not 1.0");
    }
    setPressure(p);
}

void SingleSpeciesTP::setState_PY(doublereal p, doublereal* y)
{
    if (y[0] != 1.0) {
        err("setStatePY -> x[0] not 1.0");
    }
    setPressure(p);
}

void SingleSpeciesTP::setState_HP(doublereal h, doublereal p,
                                  doublereal tol)
{
    doublereal dt;
    setPressure(p);
    for (int n = 0; n < 50; n++) {
        dt = (h - enthalpy_mass())/cp_mass();
        if (dt > 100.0) {
            dt = 100.0;
        } else if (dt < -100.0) {
            dt = -100.0;
        }
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
        dt = (u - intEnergy_mass())/cv_mass();
        if (dt > 100.0) {
            dt = 100.0;
        } else if (dt < -100.0) {
            dt = -100.0;
        }
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
        dt = (s - entropy_mass())*temperature()/cp_mass();
        if (dt > 100.0) {
            dt = 100.0;
        } else if (dt < -100.0) {
            dt = -100.0;
        }
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
        dt = (s - entropy_mass())*temperature()/cv_mass();
        if (dt > 100.0) {
            dt = 100.0;
        } else if (dt < -100.0) {
            dt = -100.0;
        }
        setTemperature(temperature() + dt);
        if (fabs(dt) < tol) {
            return;
        }
    }
    throw CanteraError("setState_SV","no convergence. dt = " + fp2str(dt));
}

doublereal SingleSpeciesTP::err(const std::string& msg) const
{
    throw CanteraError("SingleSpeciesTP","Base class method "
                       +msg+" called. Equation of state type: "
                       +int2str(eosType()));
    return 0;
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
    int leng = 1;
    m_h0_RT.resize(leng);
    m_cp0_R.resize(leng);
    m_s0_R.resize(leng);

    /*
     *  Make sure the species mole fraction is equal to 1.0;
     */
    double x = 1.0;
    setMoleFractions(&x);
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
