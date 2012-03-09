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

/*
 * --------------  Constructors ------------------------------------
 *
 */

// Base empty constructor.
/*
 *   Base constructor -> does nothing but called the inherited
 *   class constructor
 */
SingleSpeciesTP::SingleSpeciesTP() :
    ThermoPhase(),
    m_tmin(0.0),
    m_tmax(0.0),
    m_press(OneAtm),
    m_p0(OneAtm),
    m_tlast(-1.0)
{
}


//! Copy constructor
/*!
 * @param right Object to be copied
 */
SingleSpeciesTP::SingleSpeciesTP(const SingleSpeciesTP& right):
    ThermoPhase(),
    m_tmin(0.0),
    m_tmax(0.0),
    m_press(OneAtm),
    m_p0(OneAtm),
    m_tlast(-1.0)
{
    *this = operator=(right);
}

//! Assignment operator
/*!
 * @param right Object to be copied
 */
SingleSpeciesTP& SingleSpeciesTP::operator=(const SingleSpeciesTP& right)
{
    if (&right != this) {
        ThermoPhase::operator=(right);
        m_tmin       = right.m_tmin;
        m_tmax       = right.m_tmax;
        m_press      = right.m_press;
        m_p0         = right.m_p0;
        m_tlast      = right.m_tlast;
        m_h0_RT      = right.m_h0_RT;
        m_cp0_R      = right.m_cp0_R;
        m_s0_R       = right.m_s0_R;
    }
    return *this;
}

/*
 *  destructor -> does nothing but implicitly calls the inherited
 *                class destructors.
 */
SingleSpeciesTP::~SingleSpeciesTP()
{
}

//! Duplication function
/*!
 * This virtual function is used to create a duplicate of the
 * current phase. It's used to duplicate the phase when given
 * a ThermoPhase pointer to the phase.
 *
 * @return It returns a ThermoPhase pointer.
 */
ThermoPhase* SingleSpeciesTP::duplMyselfAsThermoPhase() const
{
    SingleSpeciesTP* stp = new SingleSpeciesTP(*this);
    return (ThermoPhase*) stp;
}

/**
 *
 * ------------------- Utilities ----------------------------------
 *
 */

/**
 * eosType():
 *      Creates an error because this is not a fully formed
 *      class
 */
int SingleSpeciesTP::eosType() const
{
    err("eosType");
    return -1;
}

/**
 * ------------ Molar Thermodynamic Properties --------------------
 *
 *
 *   For this single species template, the molar properties of
 *   the mixture are identified with the partial molar properties
 *   of species number 0. The partial molar property routines
 *   are called to evaluate these functions.
 */

/**
 * enthalpy_mole():
 *
 *  Molar enthalpy. Units: J/kmol.
 */
doublereal SingleSpeciesTP::enthalpy_mole() const
{
    double hbar;
    getPartialMolarEnthalpies(&hbar);
    return hbar;
}

/**
 * enthalpy_mole():
 *
 *  Molar internal energy. Units: J/kmol.
 */
doublereal SingleSpeciesTP::intEnergy_mole() const
{
    double ubar;
    getPartialMolarIntEnergies(&ubar);
    return ubar;
}

/**
 * entropy_mole():
 *
 *  Molar entropy of the mixture. Units: J/kmol/K.
 */
doublereal SingleSpeciesTP::entropy_mole() const
{
    double sbar;
    getPartialMolarEntropies(&sbar);
    return sbar;
}

/**
 * gibbs_mole():
 *
 *  Molar Gibbs free energy of the mixture. Units: J/kmol/K.
 */
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

/**
 * cp_mole():
 *
 *  Molar heat capacity at constant pressure of the mixture.
 *  Units: J/kmol/K.
 */
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

/*
 * cv_mole():
 *
 *  Molar heat capacity at constant volume of the mixture.
 *  Units: J/kmol/K.
 *
 *  For single species, we go directory to the
 *  general Cp - Cv relation
 *
 *  Cp = Cv + alpha**2 * V * T / beta
 *
 * where
 *     alpha = volume thermal expansion coefficient
 *     beta  = isothermal compressibility
 */
doublereal SingleSpeciesTP::cv_mole() const
{
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
 * ----------- Chemical Potentials and Activities ----------------------
 */

/*
 * ----------- Partial Molar Properties of the Solution -----------------
 *
 *  These are calculated by reference to the standard state properties
 *  of the zeroeth species.
 */


// Get the array of chemical potentials at unit activity
/*
 * These are the standard state chemical potentials.  \f$ \mu^0_k \f$.
 *
 *  @param mu   On return, Contains the chemical potential of the single species
 *              and the phase. Units are J / kmol . Length = 1
 */
void SingleSpeciesTP::getChemPotentials(doublereal* mu) const
{
    getStandardChemPotentials(mu);
}


//  Get the array of non-dimensional species chemical potentials
// These are partial molar Gibbs free energies.
/*
 *  These are the standard state dimensionless chemical potentials.
 *  \f$ \mu_k / \hat R T \f$.
 *
 * Units: unitless
 *
 *  @param murt   On return, Contains the chemical potential / RT of the single species
 *                and the phase. Units are unitless. Length = 1
 */
void SingleSpeciesTP::getChemPotentials_RT(doublereal* murt) const
{
    getStandardChemPotentials(murt);
    double rt = GasConstant * temperature();
    murt[0] /= rt;
}

// Get the species electrochemical potentials. Units: J/kmol.
/*
 * This method adds a term \f$ Fz_k \phi_k \f$ to
 * each chemical potential.
 *
 * This is resolved here. A single  species phase
 * is not allowed to have anything other than a zero charge.
 *
 *  @param murt   On return, Contains the chemical potential / RT of the single species
 *                and the phase. Units are unitless. Length = 1
 */
void SingleSpeciesTP::getElectrochemPotentials(doublereal* mu) const
{
    getChemPotentials(mu);
}

// Get the species partial molar enthalpies. Units: J/kmol.
/*
 * These are the phase enthalpies.  \f$ h_k \f$.
 *
 *  @param hbar On return, Contains the enthalpy of the single species
 *              and the phase. Units are J / kmol . Length = 1
 */
void SingleSpeciesTP::
getPartialMolarEnthalpies(doublereal* hbar) const
{
    double _rt = GasConstant * temperature();
    getEnthalpy_RT(hbar);
    hbar[0] *= _rt;
}

// Get the species partial molar internal energies. Units: J/kmol.
/*
 * These are the phase internal energies.  \f$ u_k \f$.
 *
 * This  member function is resolved here. A single species phase obtains its
 * thermo from the standard state function.
 *
 *  @param ubar On return, Contains the internal energy of the single species
 *              and the phase. Units are J / kmol . Length = 1
 */
void SingleSpeciesTP::
getPartialMolarIntEnergies(doublereal* ubar) const
{
    double _rt = GasConstant * temperature();
    getIntEnergy_RT(ubar);
    ubar[0] *= _rt;
}

// Get the species partial molar entropy. Units: J/kmol K.
/*
 * This is the phase entropy.  \f$ s(T,P) = s_o(T,P) \f$.
 *
 * This member function is resolved here. A single species phase obtains its
 * thermo from the standard state function.
 *
 *  @param sbar On return, Contains the entropy of the single species
 *              and the phase. Units are J / kmol / K . Length = 1
 */
void SingleSpeciesTP::
getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    sbar[0] *= GasConstant;
}

// Get the species partial molar Heat Capacities. Units: J/ kmol K.
/*
 * This is the phase heat capacity.  \f$ Cp(T,P) = Cp_o(T,P) \f$.
 *
 * This member function is resolved here. A single species phase obtains its
 * thermo from the standard state function.
 *
 *  @param cpbar On return, Contains the heat capacity of the single species
 *              and the phase. Units are J / kmol / K . Length = 1
 */
void SingleSpeciesTP::getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    cpbar[0] *= GasConstant;
}

// Get the species partial molar volumes. Units: m^3/kmol.
/*
 * This is the phase molar volume.  \f$ V(T,P) = V_o(T,P) \f$.
 *
 * This member function is resolved here. A single species phase obtains its
 * thermo from the standard state function.
 *
 *  @param vbar On return, Contains the molar volume of the single species
 *              and the phase. Units are m^3 / kmol. Length = 1
 */
void SingleSpeciesTP::getPartialMolarVolumes(doublereal* vbar) const
{
    double mw = molecularWeight(0);
    double dens = density();
    vbar[0] = mw / dens;
}

/*
 * ----- Properties of the Standard State of the Species in the Solution
 *  -----
 */

/*
 * Get the dimensional Gibbs functions for the standard
 * state of the species at the current T and P.
 */
void SingleSpeciesTP::getPureGibbs(doublereal* gpure) const
{
    getGibbs_RT(gpure);
    gpure[0] *= GasConstant * temperature();
}


// Get the molar volumes of each species in their standard
// states at the current  <I>T</I> and <I>P</I> of the solution.
/*
 *   units = m^3 / kmol
 *
 * We resolve this function at this level, by assigning
 * the molecular weight divided by the phase density
 *
 * @param vbar On output this contains the standard volume of the species
 *             and phase (m^3/kmol). Vector of length 1
 */
void SingleSpeciesTP::getStandardVolumes(doublereal* vbar) const
{
    double mw = molecularWeight(0);
    double dens = density();
    vbar[0] = mw / dens;
}

/*
 * ---- Thermodynamic Values for the Species Reference States -------
 */


/**
 *  Returns the vector of nondimensional
 *  enthalpies of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *
 *
 */
void SingleSpeciesTP::getEnthalpy_RT_ref(doublereal* hrt) const
{
    _updateThermo();
    hrt[0] = m_h0_RT[0];
}


/**
 *  Returns the vector of nondimensional
 *  enthalpies of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 */
void SingleSpeciesTP::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    grt[0] = m_h0_RT[0] - m_s0_R[0];
}

/**
 *  Returns the vector of the
 *  gibbs function of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  units = J/kmol
 */
void SingleSpeciesTP::getGibbs_ref(doublereal* g) const
{
    getGibbs_RT_ref(g);
    g[0] *= GasConstant * temperature();
}

/**
 *  Returns the vector of nondimensional
 *  entropies of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 */
void SingleSpeciesTP::getEntropy_R_ref(doublereal* er) const
{
    _updateThermo();
    er[0] = m_s0_R[0];
}

/**
 * Get the nondimensional Gibbs functions for the standard
 * state of the species at the current T and reference pressure
 * for the species.
 */
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

/*
 *  This private function throws a cantera exception. It's used when
 * this class doesn't have an answer for the question given to it,
 *  because the derived class isn't overriding a function.
 */
doublereal SingleSpeciesTP::err(std::string msg) const
{
    throw CanteraError("SingleSpeciesTP","Base class method "
                       +msg+" called. Equation of state type: "
                       +int2str(eosType()));
    return 0;
}

/*
 * @internal Initialize. This method is provided to allow
 * subclasses to perform any initialization required after all
 * species have been added. For example, it might be used to
 * resize internal work arrays that must have an entry for
 * each species.  The base class implementation does nothing,
 * and subclasses that do not require initialization do not
 * need to overload this method.  When importing a CTML phase
 * description, this method is called just prior to returning
 * from function importPhase.
 *
 * Inheriting objects should call this function
 *
 * @see importCTML.cpp
 */
void SingleSpeciesTP::initThermo()
{

    /*
     * Make sure there is one and only one species in this phase.
     */
    m_kk = nSpecies();
    if (m_kk != 1) {
        throw CanteraError("initThermo",
                           "stoichiometric substances may only contain one species.");
    }
    doublereal tmin = m_spthermo->minTemp();
    doublereal tmax = m_spthermo->maxTemp();
    if (tmin > 0.0) {
        m_tmin = tmin;
    }
    if (tmax > 0.0) {
        m_tmax = tmax;
    }

    /*
     * Store the reference pressure in the variables for the class.
     */
    m_p0 = refPressure();

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

/*
 * _updateThermo():
 *
 *        This crucial internal routine calls the species thermo
 *        update program to calculate new species Cp0, H0, and
 *        S0 whenever the temperature has changed.
 */
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




