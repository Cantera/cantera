/**
 *  @file IdealSolidSolnPhase.cpp
 *      Implementation file for an ideal solid solution model
 *      with incompressible thermodynamics (see \ref thermoprops and
 *      \link Cantera::IdealSolidSolnPhase IdealSolidSolnPhase\endlink).
 */
/*
 * Copyright 2006 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000, with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 */

#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/base/stringUtils.h"

#include <iostream>
#include <fstream>

using namespace std;

namespace Cantera
{

/*
 * Constructor for IdealSolidSolnPhase class:
 *  The default form for the generalized concentrations is 0
 * i.e., unity.
 */
IdealSolidSolnPhase::IdealSolidSolnPhase(int formGC) :
    ThermoPhase(),
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_tlast(0.0)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" IdealSolidSolnPhase Constructor",
                           " Illegal value of formGC");
    }
}

IdealSolidSolnPhase::IdealSolidSolnPhase(const std::string& inputFile,
                                         const std::string& id,
                                         int formGC) :
    ThermoPhase(),
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_tlast(0.0)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" IdealSolidSolnPhase Constructor",
                           " Illegal value of formGC");
    }
    initThermoFile(inputFile, id);
}
//====================================================================================================================
IdealSolidSolnPhase::IdealSolidSolnPhase(XML_Node& root, const std::string& id,
        int formGC) :
    ThermoPhase(),
    m_formGC(formGC),
    m_Pref(OneAtm),
    m_Pcurrent(OneAtm),
    m_tlast(0.0)
{
    if (formGC < 0 || formGC > 2) {
        throw CanteraError(" IdealSolidSolnPhase Constructor",
                           " Illegal value of formGC");
    }
    importPhase(*findXMLPhase(&root, id), this);
}
//====================================================================================================================
IdealSolidSolnPhase::IdealSolidSolnPhase(const IdealSolidSolnPhase& b)
{
    *this = b;
}
//====================================================================================================================

IdealSolidSolnPhase& IdealSolidSolnPhase::
operator=(const IdealSolidSolnPhase& b)
{
    if (this != &b) {
        ThermoPhase::operator=(b);

        m_formGC     = b.m_formGC;
        m_Pref       = b.m_Pref;
        m_Pcurrent   = b.m_Pcurrent;
        m_speciesMolarVolume = b.m_speciesMolarVolume;
        m_tlast      = b.m_tlast;
        m_h0_RT      = b.m_h0_RT;
        m_cp0_R      = b.m_cp0_R;
        m_g0_RT      = b.m_g0_RT;
        m_s0_R       = b.m_s0_R;
        m_expg0_RT   = b.m_expg0_RT;
        m_pe         = b.m_pe;
        m_pp         = b.m_pp;
    }
    return *this;
}

/*
 * Base Class Duplication Function
 *  -> given a pointer to ThermoPhase, this function can
 *     duplicate the object. (note has to be a separate function
 *     not the copy constructor, because it has to be
 *     a virtual function)
 */
ThermoPhase* IdealSolidSolnPhase::duplMyselfAsThermoPhase() const
{
    return new IdealSolidSolnPhase(*this);
}
//====================================================================================================================
/**
 * Equation of state flag. Returns the value cIdealGas, defined
 * in mix_defs.h.
 */
int IdealSolidSolnPhase::eosType() const
{
    integer res;
    switch (m_formGC) {
    case 0:
        res = cIdealSolidSolnPhase0;
        break;
    case 1:
        res = cIdealSolidSolnPhase1;
        break;
    case 2:
        res = cIdealSolidSolnPhase2;
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
    return res;
}

/********************************************************************
 *            Molar Thermodynamic Properties of the Solution
 ********************************************************************/
/**
 * Molar enthalpy of the solution. Units: J/kmol.
 * For an ideal, constant partial molar volume solution mixture with
 * pure species phases which exhibit zero volume expansivity and
 * zero isothermal compressibility:
 * \f[
 * \hat h(T,P) = \sum_k X_k \hat h^0_k(T) + (P - P_{ref}) (\sum_k X_k \hat V^0_k)
 * \f]
 * The reference-state pure-species enthalpies at the reference pressure Pref
 * \f$ \hat h^0_k(T) \f$, are computed by the species thermodynamic
 * property manager. They are polynomial functions of temperature.
 * @see SpeciesThermo
 */
doublereal IdealSolidSolnPhase::
enthalpy_mole() const
{
    const double* eptr = &(enthalpy_RT_ref()[0]);
    doublereal htp = (GasConstant * temperature() * mean_X(eptr));
    return (htp + (pressure() - m_Pref)/molarDensity());
}

/**
 * Molar internal energy of the solution. J/kmol.
 * For an ideal, constant partial molar volume solution mixture with
 * pure species phases which exhibit zero volume expansivity and
 * zero isothermal compressibility:
 * \f[
 * \hat u(T) = \hat h(T,P) - p \hat V =  \sum_k X_k \hat h^0_k(T)
 *              - P_{ref} (\sum_k X_k \hat V^0_k)
 * \f]
 * and is a function only of temperature.
 * The reference-state pure-species enthalpies
 * \f$ \hat h^0_k(T) \f$ are computed by the species thermodynamic
 * property manager.
 * @see SpeciesThermo
 */
doublereal IdealSolidSolnPhase::intEnergy_mole() const
{
    const double* eptr = DATA_PTR(enthalpy_RT_ref().begin());
    doublereal htp = (GasConstant * temperature() *
                      mean_X(eptr));
    return (htp - m_Pref / molarDensity());
}

/**
 * Molar entropy of the solution. Units: J/kmol/K.
 * For an ideal, constant partial molar volume solution mixture with
 * pure species phases which exhibit zero volume expansivity:
 * \f[
 * \hat s(T, P, X_k) = \sum_k X_k \hat s^0_k(T)
 *                   - \hat R  \sum_k X_k log(X_k)
 * \f]
 * The reference-state pure-species entropies
 * \f$ \hat s^0_k(T,p_{ref}) \f$ are computed by the species thermodynamic
 * property manager. The pure species entropies are independent of
 * temperature since the volume expansivities are equal to zero.
 * @see SpeciesThermo
 */
doublereal IdealSolidSolnPhase::entropy_mole() const
{
    const double* dptr = DATA_PTR(entropy_R_ref());
    return GasConstant * (mean_X(dptr) - sum_xlogx());
}

/**
 * Molar gibbs free energy of the solution. Units: J/kmol.
 * For an ideal, constant partial molar volume solution mixture with
 * pure species phases which exhibit zero volume expansivity:
 * \f[
 * \hat g(T, P) = \sum_k X_k \hat g^0_k(T,P) + \hat R T \sum_k X_k log(X_k)
 * \f]
 * The reference-state pure-species gibbs free energies
 * \f$ \hat g^0_k(T) \f$ are computed by the species thermodynamic
 * property manager, while the standard state gibbs free energies
 * \f$ \hat g^0_k(T,P) \f$ are computed by the member function, gibbs_RT().
 * @see SpeciesThermo
 */
doublereal IdealSolidSolnPhase::gibbs_mole() const
{
    const double* dptr = DATA_PTR(gibbs_RT_ref());
    doublereal g = mean_X(dptr);
    return (GasConstant * temperature() * (g + sum_xlogx()));
}

/**
 * Molar heat capacity at constant pressure of the solution.
 * Units: J/kmol/K.
 * For an ideal, constant partial molar volume solution mixture with
 * pure species phases which exhibit zero volume expansivity:
 * \f[
 * \hat c_p(T,P) = \sum_k X_k \hat c^0_{p,k}(T) .
 * \f]
 * The heat capacity is independent of pressure.
 * The reference-state pure-species heat capacities
 * \f$ \hat c^0_{p,k}(T) \f$ are computed by the species thermodynamic
 * property manager.
 * @see SpeciesThermo
 */
doublereal IdealSolidSolnPhase::cp_mole() const
{
    const double* dptr = DATA_PTR(cp_R_ref());
    return GasConstant * mean_X(dptr);
}

/********************************************************************
 *                  Mechanical Equation of State
 ********************************************************************/
/**
 *
 * Calculate the density of the mixture using the partial
 * molar volumes and mole fractions as input
 *
 * The formula for this is
 *
 * \f[
 * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
 * \f]
 *
 * where \f$ X_k \f$ are the mole fractions, \f$W_k\f$ are
 * the molecular weights, and \f$V_k\f$ are the pure species
 * molar volumes.
 *
 * Note, the basis behind this formula is that in an ideal
 * solution the partial molar volumes are equal to the pure
 * species molar volumes. We have additionally specified that
 * in this class that the pure species molar volumes are
 * independent of temperature and pressure.
 */
void IdealSolidSolnPhase::calcDensity()
{
    /*
     * Calculate the molarVolume of the solution (m**3 kmol-1)
     */
    const doublereal* const dtmp = moleFractdivMMW();
    double invDens = dot(m_speciesMolarVolume.begin(),
                         m_speciesMolarVolume.end(), dtmp);
    /*
     * Set the density in the parent State object directly,
     * by calling the Phase::setDensity() function.
     */
    double dens = 1.0/invDens;
    Phase::setDensity(dens);
}

/**
 * Overwritten setDensity() function is necessary because the
 * density is not an independent variable.
 *
 * This function will now throw an error condition
 *
 * @internal May have to adjust the strategy here to make
 * the eos for these materials slightly compressible, in order
 * to create a condition where the density is a function of
 * the pressure.
 *
 * This function will now throw an error condition.
 *
 *  NOTE: This is a virtual function that overwrites the State.h
 *        class
 */
void IdealSolidSolnPhase::
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
        throw CanteraError("IdealSolidSolnPhase::setDensity",
                           "Density is not an independent variable");
    }
}

/*
 * setPressure(double)               (virtual from ThermoPhase)
 *
 * Set the pressure at constant temperature. Units: Pa.
 * This method sets a constant within the object.
 * The mass density is not a function of pressure.
 * Note: This function overrides the setPressure() function
 *       in the ThermoPhase object.
 *       We calculate the density and store it in the
 *       State object, because this density is supposed to
 *       be current after setting the pressure, and is now
 *       a dependent variable.
 */
void IdealSolidSolnPhase::setPressure(doublereal p)
{
    m_Pcurrent = p;
    calcDensity();
}

/*
 * setMolarDensity()                    (virtual from State)
 * Overwritten setMolarDensity() function is necessary because the
 * density is not an independent variable.
 *
 * This function will now throw an error condition.
 *
 *  NOTE: This is a virtual function that overrides the State.h
 *        class
 */
void IdealSolidSolnPhase::setMolarDensity(const doublereal n)
{
    throw CanteraError("IdealSolidSolnPhase::setMolarDensity",
                       "Density is not an independent variable");
}

/*
 * setMoleFractions()                     (virtual from State)
 *
 * Sets the mole fractions and adjusts the internal density.
 */
void IdealSolidSolnPhase::setMoleFractions(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    calcDensity();
}

/**
 * setMoleFractions_NoNorm()               (virtual from State)
 *
 * Sets the mole fractions and adjusts the internal density.
 */
void IdealSolidSolnPhase::setMoleFractions_NoNorm(const doublereal* const x)
{
    Phase::setMoleFractions(x);
    calcDensity();
}

/*
 * setMassFractions()                      (virtual from State)
 *
 * Sets the mass fractions and adjusts the internal density.
 */
void IdealSolidSolnPhase::setMassFractions(const doublereal* const y)
{
    Phase::setMassFractions(y);
    calcDensity();
}

/*
 * setMassFractions_NoNorm()               (virtual from State)
 *
 * Sets the mass fractions and adjusts the internal density.
 */
void IdealSolidSolnPhase::setMassFractions_NoNorm(const doublereal* const y)
{
    Phase::setMassFractions_NoNorm(y);
    calcDensity();
}

/*
 * setConcentrations                        (virtual from State)
 *
 * Sets the concentrations and adjusts the internal density
 */
void IdealSolidSolnPhase::setConcentrations(const doublereal* const c)
{
    Phase::setConcentrations(c);
    calcDensity();
}

/********************************************************************
 *        Chemical Potentials and Activities
 ********************************************************************/

/********************************************************************
 *
 * getActivitConcentrations():
 *
 * This method returns the array of generalized
 * concentrations. The generalized concentrations are used
 * in the evaluation of the rates of progress for reactions
 * involving species in this phase. The generalized
 * concentration divided by the standard concentration is also
 * equal to the activity of species.
 *
 * For this implementation the activity is defined to be the
 * mole fraction of the species. The generalized concentration
 * is defined to be equal to the mole fraction divided by
 * the partial molar volume. The generalized concentrations
 * for species in this phase therefore have units of
 * kmol m<SUP>-3</SUP>. Rate constants must reflect this fact.
 *
 * On a general note, the following must be true.
 * For an ideal solution, the generalized concentration  must consist
 * of the mole fraction multiplied by a constant. The constant may be
 * fairly arbitrarily chosen, with differences adsorbed into the
 * reaction rate expression. 1/V_N, 1/V_k, or 1 are equally good,
 * as long as the standard concentration is adjusted accordingly.
 * However, it must be a constant (and not the concentration, btw,
 * which is a function of the mole fractions) in order for the
 * ideal solution properties to hold at the same time having the
 * standard concentration to be independent of the mole fractions.
 *
 * In this implementation the form of the generalized concentrations
 * depend upon the member attribute, m_formGC:
 *
 *                          <TABLE>
 *  <TR><TD> m_formGC </TD><TD> GeneralizedConc </TD><TD> StandardConc </TD></TR>
 *  <TR><TD> 0        </TD><TD> X_k             </TD><TD> 1.0          </TD></TR>
 *  <TR><TD> 1        </TD><TD> X_k / V_k       </TD><TD> 1.0 / V_k    </TD></TR>
 *  <TR><TD> 2        </TD><TD> X_k / V_N       </TD><TD> 1.0 / V_N    </TD></TR>
 *                         </TABLE>
 *
 * HKM Note: We have absorbed the pressure dependence of the pure species
 *        state into the thermodynamics functions. Therefore the
 *        standard state on which the activities are based depend
 *        on both temperature and pressure. If we hadn't, it would have
 *        appeared in this function in a very awkward exp[] format.
 *
 * @param c[] Pointer to array of doubles of length m_kk, which on exit
 *          will contain the generalized concentrations.
 */
void IdealSolidSolnPhase::
getActivityConcentrations(doublereal* c) const
{
    const doublereal* const dtmp = moleFractdivMMW();
    const double mmw = meanMolecularWeight();
    switch (m_formGC) {
    case 0:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = dtmp[k] * mmw;
        }
        break;
    case 1:
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = dtmp[k] * mmw / m_speciesMolarVolume[k];
        }
        break;
    case 2:
        double atmp = mmw / m_speciesMolarVolume[m_kk-1];
        for (size_t k = 0; k < m_kk; k++) {
            c[k] = dtmp[k] * atmp;
        }
        break;
    }
}

/*********************************************************************
 *
 * standardConcentration()
 *
 * The standard concentration \f$ C^0_k \f$ used to normalize
 * the generalized concentration.
 * In many cases, this quantity
 * will be the same for all species in a phase.
 * However, for this case, we will return a distinct concentration
 * for each species. This is the inverse of the species molar
 * volume. Units are m<SUP>3</SUP> kmol<SUP>-1</SUP>.
 *
 *
 * @param k Species number: this is a require parameter,
 * a change from the ThermoPhase base class, where it was
 * an optional parameter.
 */
doublereal IdealSolidSolnPhase::
standardConcentration(size_t k) const
{
    switch (m_formGC) {
    case 0:
        return 1.0;
    case 1:
        return 1.0 / m_speciesMolarVolume[k];
    case 2:
        return 1.0/m_speciesMolarVolume[m_kk-1];
    }
    return 0.0;
}
doublereal IdealSolidSolnPhase::
referenceConcentration(int k) const
{
    switch (m_formGC) {
    case 0:
        return 1.0;
    case 1:
        return 1.0 / m_speciesMolarVolume[k];
    case 2:
        return 1.0 / m_speciesMolarVolume[m_kk-1];
    }
    return 0.0;
}

/*********************************************************************
 *
 * logStandardConc()
 *
 * Returns the log of the standard concentration
 *
 * @param k Species number: this is a require parameter,
 * a change from the ThermoPhase base class, where it was
 * an optional parameter.
 */
doublereal IdealSolidSolnPhase::
logStandardConc(size_t k) const
{
    _updateThermo();
    double res;
    switch (m_formGC) {
    case 0:
        res = 0.0;
        break;
    case 1:
        res = log(1.0/m_speciesMolarVolume[k]);
        break;
    case 2:
        res =  log(1.0/m_speciesMolarVolume[m_kk-1]);
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
    return res;
}

/***********************************************************************
 *
 * getUnitsStandardConcentration()
 *
 * Returns the units of the standard and general concentrations
 * Note they have the same units, as their divisor is
 * defined to be equal to the activity of the kth species
 * in the solution, which is unitless.
 *
 * This routine is used in print out applications where the
 * units are needed. Usually, MKS units are assumed throughout
 * the program and in the XML input files.
 *
 *  uA[0] = kmol units - default  = 1
 *  uA[1] = m    units - default  = -nDim(), the number of spatial
 *                                dimensions in the Phase class.
 *  uA[2] = kg   units - default  = 0;
 *  uA[3] = Pa(pressure) units - default = 0;
 *  uA[4] = Temperature units - default = 0;
 *  uA[5] = time units - default = 0
 *
 *  For EOS types other than cIdealSolidSolnPhase1, the default
 *  kmol/m3 holds for standard concentration units. For
 *  cIdealSolidSolnPhase0 type, the standard concentration is
 *  unitless.
 */
void IdealSolidSolnPhase::
getUnitsStandardConc(double* uA, int, int sizeUA) const
{
    int eos = eosType();
    if (eos == cIdealSolidSolnPhase0) {
        for (int i = 0; i < sizeUA; i++) {
            uA[i] = 0.0;
        }
    } else {
        for (int i = 0; i < sizeUA; i++) {
            if (i == 0) {
                uA[0] = 1.0;
            }
            if (i == 1) {
                uA[1] = -int(nDim());
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
}

/*
 * getActivityCoefficients():
 *
 */
void IdealSolidSolnPhase::
getActivityCoefficients(doublereal* ac) const
{
    for (size_t k = 0; k < m_kk; k++) {
        ac[k] = 1.0;
    }
}
//================================================================================================
/*
 *
 * getChemPotentials():
 *
 * This function returns a vector of chemical potentials of the
 * species.
 * \f[
 *    \mu_k = \mu^o_k(T) + V_k * (p - p_o) + R T ln(X_k)
 * \f]
 *  or another way to phrase this is
 * \f[
 *    \mu_k = \mu^o_k(T,p) + R T ln(X_k)
 * \f]
 *  where \f$ \mu^o_k(T,p) = \mu^o_k(T) + V_k * (p - p_o)\f$
 *
 */
void IdealSolidSolnPhase::
getChemPotentials(doublereal* mu) const
{
    doublereal delta_p = m_Pcurrent - m_Pref;
    doublereal xx;
    doublereal RT = temperature() * GasConstant;
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = RT * (g_RT[k] + log(xx))
                + delta_p * m_speciesMolarVolume[k];
    }
}
//================================================================================================
/*
 *
 * getChemPotentials_RT()
 *
 * Get the array of non-dimensional chemical potentials \f$
 * \mu_k / \hat R T \f$, where
 *
 * \f[
 *    \mu_k = \mu^o_k(T) + V_k * (p - p_o) + R T ln(X_k)
 * \f]
 *  or another way to phrase this is
 * \f[
 *    \mu_k = \mu^o_k(T,p) + R T ln(X_k)
 * \f]
 *  where \f$ \mu^o_k(T,p) = \mu^o_k(T) + V_k * (p - p_o)\f$
 *
 */
void IdealSolidSolnPhase::
getChemPotentials_RT(doublereal* mu) const
{
    doublereal RT = temperature() * GasConstant;
    doublereal delta_pdRT = (m_Pcurrent - m_Pref) / RT;
    doublereal xx;
    const vector_fp& g_RT = gibbs_RT_ref();
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        mu[k] = (g_RT[k] + log(xx))
                + delta_pdRT * m_speciesMolarVolume[k];
    }
}

/********************************************************************
 *                    Partial Molar Properties
 ********************************************************************/

/********************************************************************
 *
 * getPartialMolarEnthalpies()
 *
 * For this phase, the partial molar enthalpies are equal to the
 * pure species enthalpies.
 *  \f[
 * \hat h_k(T,P) = \sum_k X_k \hat h^0_k(T) + (p - p_{ref}) (\sum_k X_k \hat V^0_k)
 * \f]
 * The reference-state pure-species enthalpies at the reference
 * pressure p_ref
 * \f$ \hat h^0_k(T) \f$, are computed by the species thermodynamic
 * property manager. They are polynomial functions of temperature.
 * @see SpeciesThermo
 */
void IdealSolidSolnPhase::getPartialMolarEnthalpies(doublereal* hbar) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal rt = GasConstant * temperature();
    scale(_h.begin(), _h.end(), hbar, rt);
}

/********************************************************************
 *
 * getPartialMolarEntropies()
 *
 * Returns an array of partial molar entropies of the species in the
 * solution. Units: J/kmol.
 * For this phase, the partial molar entropies are equal to the
 * pure species entropies plus the ideal solution contribution.
 *  \f[
 * \bar s_k(T,P) =  \hat s^0_k(T) - R log(X_k)
 * \f]
 * The reference-state pure-species entropies,\f$ \hat s^0_k(T) \f$,
 * at the reference pressure, \f$ P_{ref} \f$,  are computed by the
 * species thermodynamic
 * property manager. They are polynomial functions of temperature.
 * @see SpeciesThermo
 */
void IdealSolidSolnPhase::
getPartialMolarEntropies(doublereal* sbar) const
{
    const vector_fp& _s = entropy_R_ref();
    doublereal r = GasConstant;
    doublereal xx;
    for (size_t k = 0; k < m_kk; k++) {
        xx = std::max(SmallNumber, moleFraction(k));
        sbar[k] = r * (_s[k] - log(xx));
    }
}

/********************************************************************
 *
 * getPartialMolarCp()
 *
 * For this phase, the partial molar heat capacities are equal
 * to the standard state heat capacities.

 */
void IdealSolidSolnPhase::
getPartialMolarCp(doublereal* cpbar) const
{
    getCp_R(cpbar);
    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

/******************************************************************
 *
 *  getPartialMolarVolumes()
 *
 * returns an array of partial molar volumes of the species
 * in the solution. Units: m^3 kmol-1.
 *
 * For this solution, thepartial molar volumes are equal to the
 * constant species molar volumes.
 */
void IdealSolidSolnPhase::
getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

/*****************************************************************
 *    Properties of the Standard State of the Species
 *    in the Solution
 *****************************************************************/

/******************************************************************
 *
 * getPureGibbs()
 *
 * Get the Gibbs functions for the pure species
 * at the current <I>T</I> and <I>P</I> of the solution.
 * We assume an incompressible constant partial molar
 * volume here:
 * \f[
 *  \mu^0_k(T,p) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
 * \f]
 * where \f$V_k\f$ is the molar volume of pure species <I>k<\I>.
 * \f$ u^{ref}_k(T)\f$ is the chemical potential of pure
 * species <I>k<\I> at the reference pressure, \f$P_{ref}\f$.
 */
void IdealSolidSolnPhase::
getPureGibbs(doublereal* gpure) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal RT = _RT();
    const doublereal* const gk = DATA_PTR(gibbsrt);
    doublereal delta_p = (m_Pcurrent - m_Pref);
    for (size_t k = 0; k < m_kk; k++) {
        gpure[k] = RT * gk[k] + delta_p * m_speciesMolarVolume[k];
    }
}

/**
 * Get the nondimensional gibbs function for the species
 * standard states at the current T and P of the solution.
 *
 *  \f[
 *  \mu^0_k(T,P) = \mu^{ref}_k(T) + (P - P_{ref}) * V_k
 * \f]
 * where \f$V_k\f$ is the molar volume of pure species <I>k</I>.
 * \f$ \mu^{ref}_k(T)\f$ is the chemical potential of pure
 * species <I>k</I> at the reference pressure, \f$P_{ref}\f$.
 *
 * @param grt Vector of length m_kk, which on return sr[k]
 *           will contain the nondimensional
 *           standard state gibbs function for species k.
 */
void IdealSolidSolnPhase::
getGibbs_RT(doublereal* grt) const
{
    const vector_fp& gibbsrt = gibbs_RT_ref();
    doublereal RT = _RT();
    const doublereal* const gk = DATA_PTR(gibbsrt);
    doublereal delta_prt = (m_Pcurrent - m_Pref)/ RT;
    for (size_t k = 0; k < m_kk; k++) {
        grt[k] = gk[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

/********************************************************************
 *
 * getEnthalpy_RT()
 *
 * Get the array of nondimensional Enthalpy functions for the ss
 * species at the current <I>T</I> and <I>P</I> of the solution.
 * We assume an incompressible constant partial molar
 * volume here:
 * \f[
 *  h^0_k(T,P) = h^{ref}_k(T) + (P - P_{ref}) * V_k
 * \f]
 * where \f$V_k\f$ is the molar volume of pure species <I>k<\I>.
 * \f$ h^{ref}_k(T)\f$ is the enthalpy of the pure
 * species <I>k<\I> at the reference pressure, \f$P_{ref}\f$.
 */
void IdealSolidSolnPhase::
getEnthalpy_RT(doublereal* hrt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal delta_prt = ((m_Pcurrent - m_Pref) /
                            (GasConstant * temperature()));
    for (size_t k = 0; k < m_kk; k++) {
        hrt[k] = _h[k] + delta_prt * m_speciesMolarVolume[k];
    }
}

/**
 * Get the nondimensional Entropies for the species
 * standard states at the current T and P of the solution.
 *
 * Note, this is equal to the reference state entropies
 * due to the zero volume expansivity:
 * i.e., (dS/dp)_T = (dV/dT)_P = 0.0
 *
 * @param sr Vector of length m_kk, which on return sr[k]
 *           will contain the nondimensional
 *           standard state entropy of species k.
 */
void IdealSolidSolnPhase::getEntropy_R(doublereal* sr) const
{
    const vector_fp& _s = entropy_R_ref();
    copy(_s.begin(), _s.end(), sr);
}

/*
 *  Returns the vector of nondimensional
 *  internal Energies of the standard state at the current temperature
 *  of the solution and current pressure for each species.
 * \f[
 *  u^0_k(T,P) = h^{ref}_k(T)  - P_{ref} * V_k
 * \f]
 *
 *  The standard state internal energy is independent of
 *  pressure in this equation of state.
 *  (inherited from ThermoPhase.h)
 */
void IdealSolidSolnPhase::getIntEnergy_RT(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal prefrt = m_Pref / (GasConstant * temperature());
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - prefrt * m_speciesMolarVolume[k];
    }
}

/*
 * Get the nondimensional heat capacity at constant pressure
 * function for the species
 * standard states at the current T and P of the solution.
 *
 *  \f[
 *  Cp^0_k(T,P) = Cp^{ref}_k(T)
 * \f]
 * where \f$V_k\f$ is the molar volume of pure species <I>k<\I>.
 * \f$ Cp^{ref}_k(T)\f$ is the constant pressure heat capacity
 * of species <I>k<\I> at the reference pressure, \f$P_{ref}\f$.
 *
 * @param cpr Vector of length m_kk, which on return cpr[k]
 *           will contain the nondimensional
 *           constant pressure heat capacity for species k.
 */
void IdealSolidSolnPhase::getCp_R(doublereal* cpr) const
{
    const vector_fp& _cpr = cp_R_ref();
    copy(_cpr.begin(), _cpr.end(), cpr);
}

/*
 * Get the molar volumes of each species in their standard
 * states at the current
 * <I>T</I> and <I>P</I> of the solution.
 * units = m^3 / kmol
 */
void IdealSolidSolnPhase::getStandardVolumes(doublereal* vol) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), vol);
}


/*********************************************************************
 *     Thermodynamic Values for the Species Reference States
 *********************************************************************/

/*
 *  Returns the vector of non-dimensional Enthalpy function
 *  of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  Units = unitless
 */
void IdealSolidSolnPhase::getEnthalpy_RT_ref(doublereal* hrt) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        hrt[k] = m_h0_RT[k];
    }
}

/*
 *  Returns the vector of non-dimensional Gibbs function
 *  of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  Units = unitless
 */
void IdealSolidSolnPhase::getGibbs_RT_ref(doublereal* grt) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        grt[k] = m_g0_RT[k];
    }
}

/*
 *  Returns the vector of Gibbs function
 *  of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  Units = J / kmol
 */
void IdealSolidSolnPhase::getGibbs_ref(doublereal* g) const
{
    _updateThermo();
    double tmp = GasConstant * temperature();
    for (size_t k = 0; k != m_kk; k++) {
        g[k] = tmp * m_g0_RT[k];
    }
}

/*
 *  Returns the vector of nondimensional
 *  internal Energies of the standard state at the current temperature
 *  of the solution and current pressure for each species.
 *  (inherited from ThermoPhase.h)
 */
void IdealSolidSolnPhase::getIntEnergy_RT_ref(doublereal* urt) const
{
    const vector_fp& _h = enthalpy_RT_ref();
    doublereal prefrt = m_Pref / (GasConstant * temperature());
    for (size_t k = 0; k < m_kk; k++) {
        urt[k] = _h[k] - prefrt * m_speciesMolarVolume[k];
    }
}

/*
 *  Returns the vector of non-dimensional Entropy function
 *  of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  Units = unitless
 */
void IdealSolidSolnPhase::getEntropy_R_ref(doublereal* er) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        er[k] = m_s0_R[k];
    }
}

/*
 *  Returns the vector of non-dimensional Entropy function
 *  of the reference state at the current temperature
 *  of the solution and the reference pressure for the species.
 *  Units = unitless
 */
void IdealSolidSolnPhase::getCp_R_ref(doublereal* cpr) const
{
    _updateThermo();
    for (size_t k = 0; k != m_kk; k++) {
        cpr[k] = m_cp0_R[k];
    }
}

/*
 *  Returns a reference to the vector of nondimensional
 *  enthalpies of the reference state at the current temperature.
 *  Real reason for its existence is that it also checks
 *  to see if a recalculation of the reference thermodynamics
 *  functions needs to be done.
 */
const vector_fp& IdealSolidSolnPhase::enthalpy_RT_ref() const
{
    _updateThermo();
    return m_h0_RT;
}

/*
 *  Returns a reference to the vector of nondimensional
 *  enthalpies of the reference state at the current temperature.
 *  Real reason for its existence is that it also checks
 *  to see if a recalculation of the reference thermodynamics
 *  functions needs to be done.
 */
const vector_fp& IdealSolidSolnPhase::entropy_R_ref() const
{
    _updateThermo();
    return m_s0_R;
}

/*********************************************************************
 *    Utility Functions
 *********************************************************************/
/*
 * initThermo() function initializes the object for use.
 *
 * Before its invocation, the class isn't ready for calculation.
 */
void IdealSolidSolnPhase::initThermo()
{
}

/*
 * @internal
 *   Import and initialize a ThermoPhase object
 *   using an XML tree.
 *   Here we read extra information about the XML description
 *   of a phase. Regular information about elements and species
 *   and their reference state thermodynamic information
 *   have already been read at this point.
 *   For example, we do not need to call this function for
 *   ideal gas equations of state.
 *   This function is called from importPhase()
 *   after the elements and the
 *   species are initialized with default ideal solution
 *   level data.
 *
 * @param phaseNode This object must be the phase node of a
 *             complete XML tree
 *             description of the phase, including all of the
 *             species data. In other words while "phase" must
 *             point to an XML phase object, it must have
 *             sibling nodes "speciesData" that describe
 *             the species in the phase.
 * @param id   ID of the phase. If nonnull, a check is done
 *             to see if phaseNode is pointing to the phase
 *             with the correct id.
 */
void IdealSolidSolnPhase::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    string subname = "IdealSolidSolnPhase::initThermoXML";
    if (id.size() > 0) {
        string idp = phaseNode.id();
        if (idp != id) {
            throw CanteraError(subname.c_str(),
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Check on the thermo field. Must have:
     * <thermo model="IdealSolidSolution" />
     */
    if (phaseNode.hasChild("thermo")) {
        XML_Node& thNode = phaseNode.child("thermo");
        string mStringa = thNode.attrib("model");
        string mString = lowercase(mStringa);
        if (mString != "idealsolidsolution") {
            throw CanteraError(subname.c_str(),
                               "Unknown thermo model: " + mStringa);
        }
    } else {
        throw CanteraError(subname.c_str(),
                           "Unspecified thermo model");
    }

    /*
     * Form of the standard concentrations. Must have one of:
     *
     *     <standardConc model="unity" />
     *     <standardConc model="molar_volume" />
     *     <standardConc model="solvent_volume" />
     */
    if (phaseNode.hasChild("standardConc")) {
        XML_Node& scNode = phaseNode.child("standardConc");
        string formStringa = scNode.attrib("model");
        string formString = lowercase(formStringa);
        if (formString == "unity") {
            m_formGC = 0;
        } else if (formString == "molar_volume") {
            m_formGC = 1;
        } else if (formString == "solvent_volume") {
            m_formGC = 2;
        } else {
            throw CanteraError(subname.c_str(),
                               "Unknown standardConc model: " + formStringa);
        }
    } else {
        throw CanteraError(subname.c_str(),
                           "Unspecified standardConc model");
    }

    /*
     * Initialize all of the lengths now that we know how many species
     * there are in the phase.
     */
    initLengths();
    /*
     * Now go get the molar volumes
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
                                         &phaseNode.root());
    const vector<string>&sss = speciesNames();

    for (size_t k = 0; k < m_kk; k++) {
        XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
        XML_Node* ss = s->findByName("standardState");
        m_speciesMolarVolume[k] = ctml::getFloat(*ss, "molarVolume", "toSI");
    }

    /*
     * Call the base initThermo, which handles setting the initial
     * state.
     */
    ThermoPhase::initThermoXML(phaseNode, id);
}

/*
 * This internal function adjusts the lengths of arrays
 */
void IdealSolidSolnPhase::
initLengths()
{
    /*
     * Obtain the reference pressure by calling the ThermoPhase
     * function refPressure, which in turn calls the
     * species thermo reference pressure function of the
     * same name.
     */
    m_Pref = refPressure();

    m_h0_RT.resize(m_kk);
    m_g0_RT.resize(m_kk);
    m_expg0_RT.resize(m_kk);
    m_cp0_R.resize(m_kk);
    m_s0_R.resize(m_kk);
    m_pe.resize(m_kk, 0.0);
    m_pp.resize(m_kk);
    m_speciesMolarVolume.resize(m_kk);
}

/*
 * Set mixture to an equilibrium state consistent with specified
 * element potentials and temperature.
 *
 * @param lambda_RT vector of non-dimensional element potentials
 * \f$ \lambda_m/RT \f$.
 *
 */
void IdealSolidSolnPhase::
setToEquilState(const doublereal* lambda_RT)
{
    const vector_fp& grt = gibbs_RT_ref();

    // set the pressure and composition to be consistent with
    // the temperature,
    doublereal pres = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        m_pp[k] = -grt[k];
        for (size_t m = 0; m < nElements(); m++) {
            m_pp[k] += nAtoms(k,m)*lambda_RT[m];
        }
        m_pp[k] = m_Pref * exp(m_pp[k]);
        pres += m_pp[k];
    }
    doublereal* dptr = DATA_PTR(m_pp);
    setState_PX(pres, dptr);
}
//================================================================================================
/*
 *
 * speciesMolarVolume()
 *
 * Report the molar volume of species k
 *
 * units - \f$ m^3 kmol^-1 \f$
 */
double IdealSolidSolnPhase::
speciesMolarVolume(int k) const
{
    return  m_speciesMolarVolume[k];
}

/*
 *
 * getSpeciesMolarVolumes():
 *
 * Fill in a return vector containing the species molar volumes
 * units - \f$ m^3 kmol^-1 \f$
 */
void IdealSolidSolnPhase::
getSpeciesMolarVolumes(doublereal* smv) const
{
    copy(m_speciesMolarVolume.begin(), m_speciesMolarVolume.end(), smv);
}
//================================================================================================
/*
 *
 * _updateThermo()
 *
 * This function gets called for every call to functions in this
 * class. It checks to see whether the temperature has changed and
 * thus the reference thermodynamics functions for all of the species
 * must be recalculated.
 * If the temperature has changed, the species thermo manager is called
 * to recalculate G, Cp, H, and S at the current temperature.
 */
void IdealSolidSolnPhase::
_updateThermo() const
{
    doublereal tnow = temperature();
    if (m_tlast != tnow) {
        /*
         * Update the thermodynamic functions of the reference state.
         */
        m_spthermo->update(tnow, DATA_PTR(m_cp0_R), DATA_PTR(m_h0_RT),
                           DATA_PTR(m_s0_R));
        m_tlast = tnow;
        doublereal rrt = 1.0 / (GasConstant * tnow);
        doublereal deltaE;
        for (size_t k = 0; k < m_kk; k++) {
            deltaE = rrt * m_pe[k];
            m_h0_RT[k] += deltaE;
            m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
        }
        m_tlast = tnow;
    }
}
//================================================================================================
}  // end namespace Cantera
//==================================================================================================
