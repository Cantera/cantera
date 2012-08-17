/**
 *  @file DebyeHuckel.cpp
 *    Declarations for the %DebyeHuckel ThermoPhase object, which models dilute
 *    electrolyte solutions
 *    (see \ref thermoprops and \link Cantera::DebyeHuckel DebyeHuckel \endlink).
 *
 * Class %DebyeHuckel represents a dilute liquid electrolyte phase which
 * obeys the Debye Huckel formulation for nonideality.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/DebyeHuckel.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/thermo/WaterProps.h"
#include "cantera/thermo/PDSS_Water.h"

#include "cantera/base/stringUtils.h"

#include <cstring>
#include <cstdlib>
#include <fstream>

using namespace std;
using namespace ctml;

namespace Cantera
{

/*
 * Default constructor
 */
DebyeHuckel::DebyeHuckel() :
    MolalityVPSSTP(),
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9),   // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0)
{

    m_npActCoeff.resize(3);
    m_npActCoeff[0] = 0.1127;
    m_npActCoeff[1] = -0.01049;
    m_npActCoeff[2] = 1.545E-3;
}

/*
 * Working constructors
 *
 *  The two constructors below are the normal way
 *  the phase initializes itself. They are shells that call
 *  the routine initThermo(), with a reference to the
 *  XML database to get the info for the phase.
 */
DebyeHuckel::DebyeHuckel(std::string inputFile, std::string id) :
    MolalityVPSSTP(),
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9),  // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0)
{
    m_npActCoeff.resize(3);
    m_npActCoeff[0] = 0.1127;
    m_npActCoeff[1] = -0.01049;
    m_npActCoeff[2] = 1.545E-3;
    initThermoFile(inputFile, id);
}

DebyeHuckel::DebyeHuckel(XML_Node& phaseRoot, std::string id) :
    MolalityVPSSTP(),
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(3.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9),   // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0)
{
    m_npActCoeff.resize(3);
    m_npActCoeff[0] = 0.1127;
    m_npActCoeff[1] = -0.01049;
    m_npActCoeff[2] = 1.545E-3;
    importPhase(*findXMLPhase(&phaseRoot, id), this);
}

/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor
 */
DebyeHuckel::DebyeHuckel(const DebyeHuckel& b) :
    MolalityVPSSTP(),
    m_formDH(DHFORM_DILUTE_LIMIT),
    m_formGC(2),
    m_IionicMolality(0.0),
    m_maxIionicStrength(30.0),
    m_useHelgesonFixedForm(false),
    m_IionicMolalityStoich(0.0),
    m_form_A_Debye(A_DEBYE_CONST),
    m_A_Debye(1.172576),   // units = sqrt(kg/gmol)
    m_B_Debye(3.28640E9),   // units = sqrt(kg/gmol) / m
    m_waterSS(0),
    m_densWaterSS(1000.),
    m_waterProps(0)
{
    /*
     * Use the assignment operator to do the brunt
     * of the work for the copy constructor.
     */
    *this = b;
}

/*
 * operator=()
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working assignment operator
 */
DebyeHuckel& DebyeHuckel::
operator=(const DebyeHuckel& b)
{
    if (&b != this) {
        MolalityVPSSTP::operator=(b);
        m_formDH              = b.m_formDH;
        m_formGC              = b.m_formGC;
        m_Aionic              = b.m_Aionic;
        m_npActCoeff          = b.m_npActCoeff;
        m_IionicMolality      = b.m_IionicMolality;
        m_maxIionicStrength   = b.m_maxIionicStrength;
        m_useHelgesonFixedForm= b.m_useHelgesonFixedForm;
        m_IionicMolalityStoich= b.m_IionicMolalityStoich;
        m_form_A_Debye        = b.m_form_A_Debye;
        m_A_Debye             = b.m_A_Debye;
        m_B_Debye             = b.m_B_Debye;
        m_B_Dot               = b.m_B_Dot;
        m_npActCoeff          = b.m_npActCoeff;

        // This is an internal shallow copy of the PDSS_Water pointer
        m_waterSS = dynamic_cast<PDSS_Water*>(providePDSS(0)) ;
        if (!m_waterSS) {
            throw CanteraError("DebyHuckel::operator=()", "Dynamic cast to waterPDSS failed");
        }

        m_densWaterSS         = b.m_densWaterSS;

        if (m_waterProps) {
            delete m_waterProps;
            m_waterProps = 0;
        }
        if (b.m_waterProps) {
            m_waterProps = new WaterProps(m_waterSS);
        }

        m_pp                  = b.m_pp;
        m_tmpV                = b.m_tmpV;
        m_speciesCharge_Stoich= b.m_speciesCharge_Stoich;
        m_Beta_ij             = b.m_Beta_ij;
        m_lnActCoeffMolal     = b.m_lnActCoeffMolal;
        m_d2lnActCoeffMolaldT2= b.m_d2lnActCoeffMolaldT2;
    }
    return *this;
}


/*
 * ~DebyeHuckel():   (virtual)
 *
 *   Destructor for DebyeHuckel. Release objects that
 * it owns.
 */
DebyeHuckel::~DebyeHuckel()
{
    if (m_waterProps) {
        delete m_waterProps;
        m_waterProps = 0;
    }
}

/*
 *  duplMyselfAsThermoPhase():
 *
 *  This routine operates at the ThermoPhase level to
 *  duplicate the current object. It uses the copy constructor
 *  defined above.
 */
ThermoPhase* DebyeHuckel::duplMyselfAsThermoPhase() const
{
    DebyeHuckel* mtp = new DebyeHuckel(*this);
    return (ThermoPhase*) mtp;
}

/*
 * Equation of state type flag. The base class returns
 * zero. Subclasses should define this to return a unique
 * non-zero value. Constants defined for this purpose are
 * listed in mix_defs.h.
 */
int DebyeHuckel::eosType() const
{
    int res;
    switch (m_formGC) {
    case 0:
        res = cDebyeHuckel0;
        break;
    case 1:
        res = cDebyeHuckel1;
        break;
    case 2:
        res = cDebyeHuckel2;
        break;
    default:
        throw CanteraError("eosType", "Unknown type");
        break;
    }
    return res;
}

//
// -------- Molar Thermodynamic Properties of the Solution ---------------
//
/*
 * Molar enthalpy of the solution. Units: J/kmol.
 */
doublereal DebyeHuckel::enthalpy_mole() const
{
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

/*
 * Molar internal energy of the solution. Units: J/kmol.
 *
 * This is calculated from the soln enthalpy and then
 * subtracting pV.
 */
doublereal DebyeHuckel::intEnergy_mole() const
{
    double hh = enthalpy_mole();
    double pres = pressure();
    double molarV = 1.0/molarDensity();
    double uu = hh - pres * molarV;
    return uu;
}

/*
 *  Molar soln entropy at constant pressure. Units: J/kmol/K.
 *
 *  This is calculated from the partial molar entropies.
 */
doublereal DebyeHuckel::entropy_mole() const
{
    getPartialMolarEntropies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

// Molar Gibbs function. Units: J/kmol.
doublereal DebyeHuckel::gibbs_mole() const
{
    getChemPotentials(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

/*
 * Molar heat capacity at constant pressure. Units: J/kmol/K.
 *
 * Returns the solution heat capacition at constant pressure.
 * This is calculated from the partial molar heat capacities.
 */
doublereal DebyeHuckel::cp_mole() const
{
    getPartialMolarCp(DATA_PTR(m_tmpV));
    double val = mean_X(DATA_PTR(m_tmpV));
    return val;
}

/// Molar heat capacity at constant volume. Units: J/kmol/K.
doublereal DebyeHuckel::cv_mole() const
{
    //getPartialMolarCv(m_tmpV.begin());
    //return mean_X(m_tmpV.begin());
    err("not implemented");
    return 0.0;
}

//
// ------- Mechanical Equation of State Properties ------------------------
//

/*
 * Pressure. Units: Pa.
 * For this incompressible system, we return the internally stored
 * independent value of the pressure.
 */
doublereal DebyeHuckel::pressure() const
{
    return m_Pcurrent;
}

void DebyeHuckel::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
}

void DebyeHuckel::setState_TP(doublereal t, doublereal p)
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
    _updateStandardStateThermo();

    /*
     * Calculate all of the other standard volumes
     * -> note these are constant for now
     */
    calcDensity();
}

/*
 * Calculate the density of the mixture using the partial
 * molar volumes and mole fractions as input
 *
 * The formula for this is
 *
 * \f[
 * \rho = \frac{\sum_k{X_k W_k}}{\sum_k{X_k V_k}}
 * \f]
 *
 * where \f$X_k\f$ are the mole fractions, \f$W_k\f$ are
 * the molecular weights, and \f$V_k\f$ are the pure species
 * molar volumes.
 *
 * Note, the basis behind this formula is that in an ideal
 * solution the partial molar volumes are equal to the pure
 * species molar volumes. We have additionally specified
 * in this class that the pure species molar volumes are
 * independent of temperature and pressure.
 *
 */
void DebyeHuckel::calcDensity()
{
    if (m_waterSS) {

        /*
         * Store the internal density of the water SS.
         * Note, we would have to do this for all other
         * species if they had pressure dependent properties.
         */
        m_densWaterSS = m_waterSS->density();
    }
    double* vbar = &m_pp[0];
    getPartialMolarVolumes(vbar);
    double* x = &m_tmpV[0];
    getMoleFractions(x);
    doublereal vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * x[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;
    Phase::setDensity(dd);
}


/*
 * The isothermal compressibility. Units: 1/Pa.
 * The isothermal compressibility is defined as
 * \f[
 * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
 * \f]
 *
 *  It's equal to zero for this model, since the molar volume
 *  doesn't change with pressure or temperature.
 */
doublereal DebyeHuckel::isothermalCompressibility() const
{
    throw CanteraError("DebyeHuckel::isothermalCompressibility",
                       "unimplemented");
    return 0.0;
}

/*
 * The thermal expansion coefficient. Units: 1/K.
 * The thermal expansion coefficient is defined as
 *
 * \f[
 * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
 * \f]
 *
 *  It's equal to zero for this model, since the molar volume
 *  doesn't change with pressure or temperature.
 */
doublereal DebyeHuckel::thermalExpansionCoeff() const
{
    throw CanteraError("DebyeHuckel::thermalExpansionCoeff",
                       "unimplemented");
    return 0.0;
}

/*
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
 *  NOTE: This is an overwritten function from the State.h
 *        class
 */
void DebyeHuckel::setDensity(doublereal rho)
{
    double dens = density();
    if (rho != dens) {
        throw CanteraError("Idea;MolalSoln::setDensity",
                           "Density is not an independent variable");
    }
}

/*
 * Overwritten setMolarDensity() function is necessary because the
 * density is not an independent variable.
 *
 * This function will now throw an error condition.
 *
 *  NOTE: This is a virtual function, overwritten function from the State.h
 *        class
 */
void DebyeHuckel::setMolarDensity(const doublereal conc)
{
    double concI = molarDensity();
    if (conc != concI) {
        throw CanteraError("Idea;MolalSoln::setMolarDensity",
                           "molarDensity/density is not an independent variable");
    }
}

/*
 * Overwritten setTemperature(double) from State.h. This
 * function sets the temperature, and makes sure that
 * the value propagates to underlying objects.
 */
void DebyeHuckel::setTemperature(const doublereal temp)
{
    setState_TP(temp, m_Pcurrent);
}


//
// ------- Activities and Activity Concentrations
//

/*
 * This method returns an array of generalized concentrations
 * \f$ C_k\f$ that are defined such that
 * \f$ a_k = C_k / C^0_k, \f$ where \f$ C^0_k \f$
 * is a standard concentration
 * defined below.  These generalized concentrations are used
 * by kinetics manager classes to compute the forward and
 * reverse rates of elementary reactions.
 *
 * @param c Array of generalized concentrations. The
 *           units depend upon the implementation of the
 *           reaction rate expressions within the phase.
 */
void DebyeHuckel::getActivityConcentrations(doublereal* c) const
{
    double c_solvent = standardConcentration();
    getActivities(c);
    for (size_t k = 0; k < m_kk; k++) {
        c[k] *= c_solvent;
    }
}

/*
 * The standard concentration \f$ C^0_k \f$ used to normalize
 * the generalized concentration. In many cases, this quantity
 * will be the same for all species in a phase - for example,
 * for an ideal gas \f$ C^0_k = P/\hat R T \f$. For this
 * reason, this method returns a single value, instead of an
 * array.  However, for phases in which the standard
 * concentration is species-specific (e.g. surface species of
 * different sizes), this method may be called with an
 * optional parameter indicating the species.
 *
 * For the time being we will use the concentration of pure
 * solvent for the the standard concentration of all species.
 * This has the effect of making reaction rates
 * based on the molality of species proportional to the
 * molality of the species.
 */
doublereal DebyeHuckel::standardConcentration(size_t k) const
{
    double mvSolvent = m_speciesSize[m_indexSolvent];
    return 1.0 / mvSolvent;
}

/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal DebyeHuckel::logStandardConc(size_t k) const
{
    double c_solvent = standardConcentration(k);
    return log(c_solvent);
}

/*
 * Returns the units of the standard and general concentrations
 * Note they have the same units, as their divisor is
 * defined to be equal to the activity of the kth species
 * in the solution, which is unitless.
 *
 * This routine is used in print out applications where the
 * units are needed. Usually, MKS units are assumed throughout
 * the program and in the XML input files.
 *
 * On return uA contains the powers of the units (MKS assumed)
 * of the standard concentrations and generalized concentrations
 * for the kth species.
 *
 *  uA[0] = kmol units - default  = 1
 *  uA[1] = m    units - default  = -nDim(), the number of spatial
 *                                dimensions in the Phase class.
 *  uA[2] = kg   units - default  = 0;
 *  uA[3] = Pa(pressure) units - default = 0;
 *  uA[4] = Temperature units - default = 0;
 *  uA[5] = time units - default = 0
 */
void DebyeHuckel::getUnitsStandardConc(double* uA, int k, int sizeUA) const
{
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


/*
 * Get the array of non-dimensional activities at
 * the current solution temperature, pressure, and
 * solution concentration.
 * (note solvent activity coefficient is on the molar scale).
 *
 */
void DebyeHuckel::getActivities(doublereal* ac) const
{
    _updateStandardStateThermo();
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    s_update_lnMolalityActCoeff();
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            ac[k] = m_molalities[k] * exp(m_lnActCoeffMolal[k]);
        }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    ac[m_indexSolvent] =
        exp(m_lnActCoeffMolal[m_indexSolvent]) * xmolSolvent;
}

/*
 * getMolalityActivityCoefficients()             (virtual, const)
 *
 * Get the array of non-dimensional Molality based
 * activity coefficients at
 * the current solution temperature, pressure, and
 * solution concentration.
 * (note solvent activity coefficient is on the molar scale).
 *
 *  Note, most of the work is done in an internal private routine
 */
void DebyeHuckel::
getMolalityActivityCoefficients(doublereal* acMolality) const
{
    _updateStandardStateThermo();
    A_Debye_TP(-1.0, -1.0);
    s_update_lnMolalityActCoeff();
    copy(m_lnActCoeffMolal.begin(), m_lnActCoeffMolal.end(), acMolality);
    for (size_t k = 0; k < m_kk; k++) {
        acMolality[k] = exp(acMolality[k]);
    }
}

//
// ------ Partial Molar Properties of the Solution -----------------
//
/*
 * Get the species chemical potentials. Units: J/kmol.
 *
 * This function returns a vector of chemical potentials of the
 * species in solution.
 *
 * \f[
 *    \mu_k = \mu^{o}_k(T,P) + R T ln(m_k)
 * \f]
 *
 * \f[
 *    \mu_solvent = \mu^{o}_solvent(T,P) +
 *            R T ((X_solvent - 1.0) / X_solvent)
 * \f]
 */
void DebyeHuckel::getChemPotentials(doublereal* mu) const
{
    double xx;
    const double xxSmall = 1.0E-150;
    /*
     * First get the standard chemical potentials in
     * molar form.
     *  -> this requires updates of standard state as a function
     *     of T and P
     */
    getStandardChemPotentials(mu);
    /*
     * Update the activity coefficients
     * This also updates the internal molality array.
     */
    s_update_lnMolalityActCoeff();
    doublereal RT = GasConstant * temperature();
    double xmolSolvent = moleFraction(m_indexSolvent);
    for (size_t k = 0; k < m_kk; k++) {
        if (m_indexSolvent != k) {
            xx = std::max(m_molalities[k], xxSmall);
            mu[k] += RT * (log(xx) + m_lnActCoeffMolal[k]);
        }
    }
    xx = std::max(xmolSolvent, xxSmall);
    mu[m_indexSolvent] +=
        RT * (log(xx) + m_lnActCoeffMolal[m_indexSolvent]);
}


/*
 * Returns an array of partial molar enthalpies for the species
 * in the mixture.
 * Units (J/kmol)
 *
 * We calculate this quantity partially from the relation and
 * partially by calling the standard state enthalpy function.
 *
 *     hbar_i = - T**2 * d(chemPot_i/T)/dT
 *
 * We calculate
 */
void DebyeHuckel::getPartialMolarEnthalpies(doublereal* hbar) const
{
    /*
     * Get the nondimensional standard state enthalpies
     */
    getEnthalpy_RT(hbar);
    /*
     * Dimensionalize it.
     */
    double T = temperature();
    double RT = GasConstant * T;
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT;
    }
    /*
     * Check to see whether activity coefficients are temperature
     * dependent. If they are, then calculate the their temperature
     * derivatives and add them into the result.
     */
    double dAdT = dA_DebyedT_TP();
    if (dAdT != 0.0) {
        /*
         * Update the activity coefficients, This also update the
         * internally stored molalities.
         */
        s_update_lnMolalityActCoeff();
        s_update_dlnMolalityActCoeff_dT();
        double RTT = GasConstant * T * T;
        for (size_t k = 0; k < m_kk; k++) {
            hbar[k] -= RTT * m_dlnActCoeffMolaldT[k];
        }
    }
}

/*
 *
 * getPartialMolarEntropies()        (virtual, const)
 *
 * Returns an array of partial molar entropies of the species in the
 * solution. Units: J/kmol.
 *
 * Maxwell's equations provide an insight in how to calculate this
 * (p.215 Smith and Van Ness)
 *
 *      d(chemPot_i)/dT = -sbar_i
 *
 * For this phase, the partial molar entropies are equal to the
 * SS species entropies plus the ideal solution contribution.following
 * contribution:
 *  \f[
 * \bar s_k(T,P) =  \hat s^0_k(T) - R log(M0 * molality[k])
 * \f]
 * \f[
 * \bar s_solvent(T,P) =  \hat s^0_solvent(T)
 *             - R ((xmolSolvent - 1.0) / xmolSolvent)
 * \f]
 *
 * The reference-state pure-species entropies,\f$ \hat s^0_k(T) \f$,
 * at the reference pressure, \f$ P_{ref} \f$,  are computed by the
 * species thermodynamic
 * property manager. They are polynomial functions of temperature.
 * @see SpeciesThermo
 */
void DebyeHuckel::
getPartialMolarEntropies(doublereal* sbar) const
{
    /*
     * Get the standard state entropies at the temperature
     * and pressure of the solution.
     */
    getEntropy_R(sbar);
    /*
     * Dimensionalize the entropies
     */
    doublereal R = GasConstant;
    for (size_t k = 0; k < m_kk; k++) {
        sbar[k] *= R;
    }
    /*
     * Update the activity coefficients, This also update the
     * internally stored molalities.
     */
    s_update_lnMolalityActCoeff();
    /*
     * First we will add in the obvious dependence on the T
     * term out front of the log activity term
     */
    doublereal mm;
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            mm = std::max(SmallNumber, m_molalities[k]);
            sbar[k] -= R * (log(mm) + m_lnActCoeffMolal[k]);
        }
    }
    double xmolSolvent = moleFraction(m_indexSolvent);
    mm = std::max(SmallNumber, xmolSolvent);
    sbar[m_indexSolvent] -= R *(log(mm) + m_lnActCoeffMolal[m_indexSolvent]);
    /*
     * Check to see whether activity coefficients are temperature
     * dependent. If they are, then calculate the their temperature
     * derivatives and add them into the result.
     */
    double dAdT = dA_DebyedT_TP();
    if (dAdT != 0.0) {
        s_update_dlnMolalityActCoeff_dT();
        double RT = R * temperature();
        for (size_t k = 0; k < m_kk; k++) {
            sbar[k] -= RT * m_dlnActCoeffMolaldT[k];
        }
    }
}

/*
 * getPartialMolarVolumes()                (virtual, const)
 *
 * returns an array of partial molar volumes of the species
 * in the solution. Units: m^3 kmol-1.
 *
 * For this solution, the partial molar volumes are normally
 *  equal to theconstant species molar volumes, except
 * when the activity coefficients depend on pressure.
 *
 * The general relation is
 *
 *       vbar_i = d(chemPot_i)/dP at const T, n
 *
 *              = V0_i + d(Gex)/dP)_T,M
 *
 *              = V0_i + RT d(lnActCoeffi)dP _T,M
 *
 */
void DebyeHuckel::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
    /*
     * Update the derivatives wrt the activity coefficients.
     */
    s_update_lnMolalityActCoeff();
    s_update_dlnMolalityActCoeff_dP();
    double T = temperature();
    double RT = GasConstant * T;
    for (size_t k = 0; k < m_kk; k++) {
        vbar[k] += RT * m_dlnActCoeffMolaldP[k];
    }
}

/*
 * Partial molar heat capacity of the solution:
 *   The kth partial molar heat capacity is  equal to
 *   the temperature derivative of the partial molar
 *   enthalpy of the kth species in the solution at constant
 *   P and composition (p. 220 Smith and Van Ness).
 *
 *     Cp = -T d2(chemPot_i)/dT2
 */
void DebyeHuckel::getPartialMolarCp(doublereal* cpbar) const
{
    /*
     * Get the nondimensional gibbs standard state of the
     * species at the T and P of the solution.
     */
    getCp_R(cpbar);

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }

    /*
     * Check to see whether activity coefficients are temperature
     * dependent. If they are, then calculate the their temperature
     * derivatives and add them into the result.
     */
    double dAdT = dA_DebyedT_TP();
    if (dAdT != 0.0) {
        /*
         * Update the activity coefficients, This also update the
         * internally stored molalities.
         */
        s_update_lnMolalityActCoeff();
        s_update_dlnMolalityActCoeff_dT();
        s_update_d2lnMolalityActCoeff_dT2();
        double T = temperature();
        double RT = GasConstant * T;
        double RTT = RT * T;
        for (size_t k = 0; k < m_kk; k++) {
            cpbar[k] -= (2.0 * RT * m_dlnActCoeffMolaldT[k] +
                         RTT * m_d2lnActCoeffMolaldT2[k]);
        }
    }
}




/*
 *  -------------- Utilities -------------------------------
 */

/*
 *  Initialization routine for a DebyeHuckel phase.
 *
 * This is a virtual routine. This routine will call initThermo()
 * for the parent class as well.
 */
void DebyeHuckel::initThermo()
{
    MolalityVPSSTP::initThermo();
    initLengths();
}

//! Utility function to assign an integer value from a string for the ElectrolyteSpeciesType field.
/*!
 *  @param estString  input string that will be interpreted
 */
static int interp_est(std::string estString)
{
    const char* cc = estString.c_str();
    string lc = lowercase(estString);
    const char* ccl = lc.c_str();
    if (!strcmp(ccl, "solvent")) {
        return cEST_solvent;
    } else if (!strcmp(ccl, "chargedspecies")) {
        return cEST_chargedSpecies;
    } else if (!strcmp(ccl, "weakacidassociated")) {
        return cEST_weakAcidAssociated;
    } else if (!strcmp(ccl, "strongacidassociated")) {
        return cEST_strongAcidAssociated;
    } else if (!strcmp(ccl, "polarneutral")) {
        return cEST_polarNeutral;
    } else if (!strcmp(ccl, "nonpolarneutral")) {
        return cEST_nonpolarNeutral;
    }
    int retn, rval;
    if ((retn = sscanf(cc, "%d", &rval)) != 1) {
        return -1;
    }
    return rval;
}

/*
 * Process the XML file after species are set up.
 *
 *  This gets called from importPhase(). It processes the XML file
 *  after the species are set up. This is the main routine for
 *  reading in activity coefficient parameters.
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
void DebyeHuckel::
initThermoXML(XML_Node& phaseNode, std::string id)
{
    if (id.size() > 0) {
        std::string idp = phaseNode.id();
        if (idp != id) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("DebyeHuckel::initThermoXML",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Determine the form of the Debye-Huckel model,
     * m_formDH.  We will use this information to size arrays below.
     */
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& scNode = thermoNode.child("activityCoefficients");
        m_formDH = DHFORM_DILUTE_LIMIT;
        std::string formString = scNode.attrib("model");
        if (formString != "") {
            if (formString == "Dilute_limit") {
                m_formDH = DHFORM_DILUTE_LIMIT;
            } else if (formString == "Bdot_with_variable_a") {
                m_formDH = DHFORM_BDOT_AK  ;
            } else if (formString == "Bdot_with_common_a") {
                m_formDH = DHFORM_BDOT_ACOMMON;
            } else if (formString == "Beta_ij") {
                m_formDH = DHFORM_BETAIJ;
            } else if (formString == "Pitzer_with_Beta_ij") {
                m_formDH = DHFORM_PITZER_BETAIJ;
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Unknown standardConc model: " + formString);
            }
        }
    } else {
        /*
         * If there is no XML node named "activityCoefficients", assume
         * that we are doing the extreme dilute limit assumption
         */
        m_formDH = DHFORM_DILUTE_LIMIT;
    }

    std::string stemp;

    /*
     * Possibly change the form of the standard concentrations
     */
    if (thermoNode.hasChild("standardConc")) {
        XML_Node& scNode = thermoNode.child("standardConc");
        m_formGC = 2;
        std::string formString = scNode.attrib("model");
        if (formString != "") {
            if (formString == "unity") {
                m_formGC = 0;
                printf("exit standardConc = unity not done\n");
                exit(EXIT_FAILURE);
            } else if (formString == "molar_volume") {
                m_formGC = 1;
                printf("exit standardConc = molar_volume not done\n");
                exit(EXIT_FAILURE);
            } else if (formString == "solvent_volume") {
                m_formGC = 2;
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Unknown standardConc model: " + formString);
            }
        }
    }

    /*
     * Reconcile the solvent name and index.
     */
    /*
     * Get the Name of the Solvent:
     *      <solvent> solventName </solvent>
     */
    std::string solventName = "";
    if (thermoNode.hasChild("solvent")) {
        XML_Node& scNode = thermoNode.child("solvent");
        vector<std::string> nameSolventa;
        getStringArray(scNode, nameSolventa);
        int nsp = static_cast<int>(nameSolventa.size());
        if (nsp != 1) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "badly formed solvent XML node");
        }
        solventName = nameSolventa[0];
    }
    for (size_t k = 0; k < m_kk; k++) {
        std::string sname = speciesName(k);
        if (solventName == sname) {
            m_indexSolvent = k;
            break;
        }
    }
    if (m_indexSolvent == npos) {
        cout << "DebyeHuckel::initThermoXML: Solvent Name not found"
             << endl;
        throw CanteraError("DebyeHuckel::initThermoXML",
                           "Solvent name not found");
    }
    if (m_indexSolvent != 0) {
        throw CanteraError("DebyeHuckel::initThermoXML",
                           "Solvent " + solventName +
                           " should be first species");
    }

    /*
     * Initialize all of the lengths of arrays in the object
     * now that we know what species are in the phase.
     */
    initThermo();

    /*
     * Now go get the specification of the standard states for
     * species in the solution. This includes the molar volumes
     * data blocks for incompressible species.
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB =
        get_XML_NameID("speciesData", speciesList["datasrc"],
                       &phaseNode.root());
    const vector<string>&sss = speciesNames();

    for (size_t k = 0; k < m_kk; k++) {
        XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
        if (!s) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "Species Data Base " + sss[k] + " not found");
        }
        XML_Node* ss = s->findByName("standardState");
        if (!ss) {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "Species " + sss[k] +
                               " standardState XML block  not found");
        }
        std::string modelStringa = ss->attrib("model");
        if (modelStringa == "") {
            throw CanteraError("DebyeHuckel::initThermoXML",
                               "Species " + sss[k] +
                               " standardState XML block model attribute not found");
        }
        std::string modelString = lowercase(modelStringa);

        if (k == 0) {
            if (modelString == "wateriapws" || modelString == "real_water" ||
                    modelString == "waterpdss") {
                /*
                 * Initialize the water standard state model
                 */
                m_waterSS = dynamic_cast<PDSS_Water*>(providePDSS(0)) ;
                if (!m_waterSS) {
                    throw CanteraError("HMWSoln::installThermoXML",
                                       "Dynamic cast to PDSS_Water failed");
                }
                /*
                 * Fill in the molar volume of water (m3/kmol)
                 * at standard conditions to fill in the m_speciesSize entry
                 * with something reasonable.
                 */
                m_waterSS->setState_TP(300., OneAtm);
                double dens = m_waterSS->density();
                double mw = m_waterSS->molecularWeight();
                m_speciesSize[0] = mw / dens;
#ifdef DEBUG_MODE_NOT
                cout << "Solvent species " << sss[k] << " has volume " <<
                     m_speciesSize[k] << endl;
#endif
            } else if (modelString == "constant_incompressible") {
                m_speciesSize[k] = getFloat(*ss, "molarVolume", "toSi");
#ifdef DEBUG_MODE_NOT
                cout << "species " << sss[k] << " has volume " <<
                     m_speciesSize[k] << endl;
#endif
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Solvent SS Model \"" + modelStringa +
                                   "\" is not known");
            }
        } else {
            if (modelString != "constant_incompressible") {
                throw CanteraError("DebyeHuckel::initThermoXML",
                                   "Solute SS Model \"" + modelStringa +
                                   "\" is not known");
            }
            m_speciesSize[k] = getFloat(*ss, "molarVolume", "toSI");
#ifdef DEBUG_MODE_NOT
            cout << "species " << sss[k] << " has volume " <<
                 m_speciesSize[k] << endl;
#endif
        }

    }

    /*
     * Go get all of the coefficients and factors in the
     * activityCoefficients XML block
     */
    XML_Node* acNodePtr = 0;
    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        acNodePtr = &acNode;
        /*
         * Look for parameters for A_Debye
         */
        if (acNode.hasChild("A_Debye")) {
            XML_Node* ss = acNode.findByName("A_Debye");
            string modelStringa = ss->attrib("model");
            string modelString = lowercase(modelStringa);
            if (modelString != "") {
                if (modelString == "water") {
                    m_form_A_Debye = A_DEBYE_WATER;
                } else {
                    throw CanteraError("DebyeHuckel::initThermoXML",
                                       "A_Debye Model \"" + modelStringa +
                                       "\" is not known");
                }
            } else {
                m_A_Debye = getFloat(acNode, "A_Debye");
#ifdef DEBUG_HKM_NOT
                cout << "A_Debye = " << m_A_Debye << endl;
#endif
            }
        }

        /*
         * Initialize the water property calculator. It will share
         * the internal eos water calculator.
         */
        if (m_form_A_Debye == A_DEBYE_WATER) {
            if (m_waterProps) {
                delete m_waterProps;
            }
            m_waterProps = new WaterProps(m_waterSS);
        }

        /*
         * Look for parameters for B_Debye
         */
        if (acNode.hasChild("B_Debye")) {
            m_B_Debye = getFloat(acNode, "B_Debye");
#ifdef DEBUG_HKM_NOT
            cout << "B_Debye = " << m_B_Debye << endl;
#endif
        }

        /*
         * Look for parameters for B_dot
         */
        if (acNode.hasChild("B_dot")) {
            if (m_formDH == DHFORM_BETAIJ ||
                    m_formDH == DHFORM_DILUTE_LIMIT ||
                    m_formDH == DHFORM_PITZER_BETAIJ) {
                throw CanteraError("DebyeHuckel:init",
                                   "B_dot entry in the wrong DH form");
            }
            double bdot_common = getFloat(acNode, "B_dot");
#ifdef DEBUG_HKM_NOT
            cout << "B_dot = " << bdot_common << endl;
#endif
            /*
             * Set B_dot parameters for charged species
             */
            for (size_t k = 0; k < m_kk; k++) {
                double z_k = charge(k);
                if (fabs(z_k) > 0.0001) {
                    m_B_Dot[k] = bdot_common;
                } else {
                    m_B_Dot[k] = 0.0;
                }
            }
        }

        /*
         * Look for Parameters for the Maximum Ionic Strength
         */
        if (acNode.hasChild("maxIonicStrength")) {
            m_maxIionicStrength = getFloat(acNode, "maxIonicStrength");
#ifdef DEBUG_HKM_NOT
            cout << "m_maxIionicStrength = "
                 <<m_maxIionicStrength << endl;
#endif
        }

        /*
         * Look for Helgeson Parameters
         */
        if (acNode.hasChild("UseHelgesonFixedForm")) {
            m_useHelgesonFixedForm = true;
        } else {
            m_useHelgesonFixedForm = false;
        }

        /*
         * Look for parameters for the Ionic radius
         */
        if (acNode.hasChild("ionicRadius")) {
            XML_Node& irNode = acNode.child("ionicRadius");

            std::string Aunits = "";
            double Afactor = 1.0;
            if (irNode.hasAttrib("units")) {
                std::string Aunits = irNode.attrib("units");
                Afactor = toSI(Aunits);
            }

            if (irNode.hasAttrib("default")) {
                std::string ads = irNode.attrib("default");
                double ad = fpValue(ads);
                for (size_t k = 0; k < m_kk; k++) {
                    m_Aionic[k] = ad * Afactor;
                }
            }

            /*
             * If the Debye-Huckel form is BDOT_AK, we can
             * have separate values for the denominator's ionic
             * size. -> That's how the activity coefficient is
             * parameterized. In this case only do we allow the
             * code to read in these parameters.
             */
            if (m_formDH == DHFORM_BDOT_AK) {
                /*
                 * Define a string-string map, and interpret the
                 * value of the xml element as binary pairs separated
                 * by colons, e.g.:
                 *      Na+:3.0
                 *      Cl-:4.0
                 *      H+:9.0
                 *      OH-:3.5
                 * Read them into the map.
                 */
                map<string, string> m;
                getMap(irNode, m);
                /*
                 * Iterate over the map pairs, interpreting the
                 * first string as a species in the current phase.
                 * If no match is made, silently ignore the
                 * lack of agreement (HKM -> may be changed in the
                 * future).
                 */
                map<std::string,std::string>::const_iterator _b = m.begin();
                for (; _b != m.end(); ++_b) {
                    size_t kk = speciesIndex(_b->first);
                    m_Aionic[kk] = fpValue(_b->second) * Afactor;

                }
            }
        }
        /*
         * Get the matrix of coefficients for the Beta
         * binary interaction parameters. We assume here that
         * this matrix is symmetric, so that we only have to
         * input 1/2 of the values.
         */
        if (acNode.hasChild("DHBetaMatrix")) {
            if (m_formDH == DHFORM_BETAIJ ||
                    m_formDH == DHFORM_PITZER_BETAIJ) {
                XML_Node& irNode = acNode.child("DHBetaMatrix");
                const vector<string>& sn = speciesNames();
                getMatrixValues(irNode, sn, sn, m_Beta_ij, true, true);
            } else {
                throw CanteraError("DebyeHuckel::initThermoXML:",
                                   "DHBetaMatrix found for wrong type");
            }
        }

        /*
         * Fill in parameters for the calculation of the
         * stoichiometric Ionic Strength
         *
         * The default is that stoich charge is the same as the
         * regular charge.
         */
        m_speciesCharge_Stoich.resize(m_kk, 0.0);
        for (size_t k = 0; k < m_kk; k++) {
            m_speciesCharge_Stoich[k] = m_speciesCharge[k];
        }
        /*
         * First look at the species database.
         *  -> Look for the subelement "stoichIsMods"
         *     in each of the species SS databases.
         */
        std::vector<const XML_Node*> xspecies= speciesData();
        std::string kname, jname;
        size_t jj = xspecies.size();
        for (size_t k = 0; k < m_kk; k++) {
            size_t jmap = npos;
            kname = speciesName(k);
            for (size_t j = 0; j < jj; j++) {
                const XML_Node& sp = *xspecies[j];
                jname = sp["name"];
                if (jname == kname) {
                    jmap = j;
                    break;
                }
            }
            if (jmap != npos) {
                const XML_Node& sp = *xspecies[jmap];
                if (sp.hasChild("stoichIsMods")) {
                    double val = getFloat(sp, "stoichIsMods");
                    m_speciesCharge_Stoich[k] = val;
                }
            }
        }

        /*
         * Now look at the activity coefficient database
         */
        if (acNodePtr) {
            if (acNodePtr->hasChild("stoichIsMods")) {
                XML_Node& sIsNode = acNodePtr->child("stoichIsMods");

                map<std::string, std::string> msIs;
                getMap(sIsNode, msIs);
                map<std::string,std::string>::const_iterator _b = msIs.begin();
                for (; _b != msIs.end(); ++_b) {
                    size_t kk = speciesIndex(_b->first);
                    double val = fpValue(_b->second);
                    m_speciesCharge_Stoich[kk] = val;
                }
            }
        }
    }

    /*
     * Fill in the vector specifying the electrolyte species
     * type
     *
     *   First fill in default values. Everything is either
     *   a charge species, a nonpolar neutral, or the solvent.
     */
    for (size_t k = 0; k < m_kk; k++) {
        if (fabs(m_speciesCharge[k]) > 0.0001) {
            m_electrolyteSpeciesType[k] = cEST_chargedSpecies;
            if (fabs(m_speciesCharge_Stoich[k] - m_speciesCharge[k])
                    > 0.0001) {
                m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
            }
        } else if (fabs(m_speciesCharge_Stoich[k]) > 0.0001) {
            m_electrolyteSpeciesType[k] = cEST_weakAcidAssociated;
        } else {
            m_electrolyteSpeciesType[k] = cEST_nonpolarNeutral;
        }
    }
    m_electrolyteSpeciesType[m_indexSolvent] = cEST_solvent;
    /*
     * First look at the species database.
     *  -> Look for the subelement "stoichIsMods"
     *     in each of the species SS databases.
     */
    std::vector<const XML_Node*> xspecies= speciesData();
    const XML_Node* spPtr = 0;
    std::string kname;
    for (size_t k = 0; k < m_kk; k++) {
        kname = speciesName(k);
        spPtr = xspecies[k];
        if (!spPtr) {
            if (spPtr->hasChild("electrolyteSpeciesType")) {
                std::string est = getChildValue(*spPtr, "electrolyteSpeciesType");
                if ((m_electrolyteSpeciesType[k] = interp_est(est)) == -1) {
                    throw CanteraError("DebyeHuckel:initThermoXML",
                                       "Bad electrolyte type: " + est);
                }
            }
        }
    }
    /*
     * Then look at the phase thermo specification
     */
    if (acNodePtr) {
        if (acNodePtr->hasChild("electrolyteSpeciesType")) {
            XML_Node& ESTNode = acNodePtr->child("electrolyteSpeciesType");
            map<std::string, std::string> msEST;
            getMap(ESTNode, msEST);
            map<std::string,std::string>::const_iterator _b = msEST.begin();
            for (; _b != msEST.end(); ++_b) {
                size_t kk = speciesIndex(_b->first);
                std::string est = _b->second;
                if ((m_electrolyteSpeciesType[kk] = interp_est(est))  == -1) {
                    throw CanteraError("DebyeHuckel:initThermoXML",
                                       "Bad electrolyte type: " + est);
                }
            }
        }
    }



    /*
     * Lastly set the state
     */
    if (phaseNode.hasChild("state")) {
        XML_Node& stateNode = phaseNode.child("state");
        setStateFromXML(stateNode);
    }

}

/*
 * @internal
 * Set equation of state parameters. The number and meaning of
 * these depends on the subclass.
 * @param n number of parameters
 * @param c array of \i n coefficients
 *
 */
void DebyeHuckel::setParameters(int n, doublereal* const c)
{
}

void DebyeHuckel::getParameters(int& n, doublereal* const c) const
{
}

/*
 * Set equation of state parameter values from XML
 * entries. This method is called by function importPhase in
 * file importCTML.cpp when processing a phase definition in
 * an input file. It should be overloaded in subclasses to set
 * any parameters that are specific to that particular phase
 * model.
 *
 * @param eosdata An XML_Node object corresponding to
 * the "thermo" entry for this phase in the input file.
 *
 * HKM -> Right now, the parameters are set elsewhere (initThermoXML)
 *        It just didn't seem to fit.
 */
void DebyeHuckel::setParametersFromXML(const XML_Node& eosdata)
{
}

/*
 * Report the molar volume of species k
 *
 * units - \f$ m^3 kmol^-1 \f$
 */
// double DebyeHuckel::speciesMolarVolume(int k) const {
// return m_speciesSize[k];
//}


/*
 * A_Debye_TP()                              (virtual)
 *
 *   Returns the A_Debye parameter as a function of temperature
 *  and pressure.
 *
 *  The default is to assume that it is constant, given
 *  in the initialization process and stored in the
 *  member double, m_A_Debye
 */
double DebyeHuckel::A_Debye_TP(double tempArg, double presArg) const
{
    double T = temperature();
    double A;
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }

    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        A = m_A_Debye;
        break;
    case A_DEBYE_WATER:
        A = m_waterProps->ADebye(T, P, 0);
        m_A_Debye = A;
        break;
    default:
        printf("shouldn't be here\n");
        exit(EXIT_FAILURE);
    }
    return A;
}

/*
 * dA_DebyedT_TP()                              (virtual)
 *
 *  Returns the derivative of the A_Debye parameter with
 *  respect to temperature as a function of temperature
 *  and pressure.
 *
 * units = A_Debye has units of sqrt(gmol kg-1).
 *         Temp has units of Kelvin.
 */
double DebyeHuckel::dA_DebyedT_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double dAdT;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdT = 0.0;
        break;
    case A_DEBYE_WATER:
        dAdT = m_waterProps->ADebye(T, P, 1);
        break;
    default:
        printf("shouldn't be here\n");
        exit(EXIT_FAILURE);
    }
    return dAdT;
}

/*
 * d2A_DebyedT2_TP()                              (virtual)
 *
 *  Returns the 2nd derivative of the A_Debye parameter with
 *  respect to temperature as a function of temperature
 *  and pressure.
 *
 * units = A_Debye has units of sqrt(gmol kg-1).
 *         Temp has units of Kelvin.
 */
double DebyeHuckel::d2A_DebyedT2_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double d2AdT2;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        d2AdT2 = 0.0;
        break;
    case A_DEBYE_WATER:
        d2AdT2 = m_waterProps->ADebye(T, P, 2);
        break;
    default:
        printf("shouldn't be here\n");
        exit(EXIT_FAILURE);
    }
    return d2AdT2;
}

/*
 * dA_DebyedP_TP()                              (virtual)
 *
 *  Returns the derivative of the A_Debye parameter with
 *  respect to pressure, as a function of temperature
 *  and pressure.
 *
 * units = A_Debye has units of sqrt(gmol kg-1).
 *         Pressure has units of pascals.
 */
double DebyeHuckel::dA_DebyedP_TP(double tempArg, double presArg) const
{
    double T = temperature();
    if (tempArg != -1.0) {
        T = tempArg;
    }
    double P = pressure();
    if (presArg != -1.0) {
        P = presArg;
    }
    double dAdP;
    switch (m_form_A_Debye) {
    case A_DEBYE_CONST:
        dAdP = 0.0;
        break;
    case A_DEBYE_WATER:
        dAdP = m_waterProps->ADebye(T, P, 3);
        break;
    default:
        printf("shouldn't be here\n");
        exit(EXIT_FAILURE);
    }
    return dAdP;
}

/*
 * ----------- Critical State Properties --------------------------
 */

/*
 * ---------- Other Property Functions
 */
double DebyeHuckel::AionicRadius(int k) const
{
    return m_Aionic[k];
}

/*
 * ------------ Private and Restricted Functions ------------------
 */

/*
 * Bail out of functions with an error exit if they are not
 * implemented.
 */
doublereal DebyeHuckel::err(std::string msg) const
{
    throw CanteraError("DebyeHuckel",
                       "Unfinished func called: " + msg);
    return 0.0;
}


/*
 * initLengths():
 *
 * This internal function adjusts the lengths of arrays based on
 * the number of species
 */
void DebyeHuckel::initLengths()
{
    m_kk = nSpecies();

    /*
     * Obtain the limits of the temperature from the species
     * thermo handler's limits.
     */
    m_electrolyteSpeciesType.resize(m_kk, cEST_polarNeutral);
    m_speciesSize.resize(m_kk);
    m_Aionic.resize(m_kk, 0.0);
    m_lnActCoeffMolal.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldT.resize(m_kk, 0.0);
    m_d2lnActCoeffMolaldT2.resize(m_kk, 0.0);
    m_dlnActCoeffMolaldP.resize(m_kk, 0.0);
    m_B_Dot.resize(m_kk, 0.0);
    m_pp.resize(m_kk, 0.0);
    m_tmpV.resize(m_kk, 0.0);
    if (m_formDH == DHFORM_BETAIJ ||
            m_formDH == DHFORM_PITZER_BETAIJ) {
        m_Beta_ij.resize(m_kk, m_kk, 0.0);
    }
}

/*
 * nonpolarActCoeff()                    (private)
 *
 *   Static function that implements the non-polar species
 *   salt-out modifications.
 *   Returns the calculated activity coefficients.
 */
double DebyeHuckel::_nonpolarActCoeff(double IionicMolality) const
{
    double I2 = IionicMolality * IionicMolality;
    double l10actCoeff =
        m_npActCoeff[0] * IionicMolality +
        m_npActCoeff[1] * I2 +
        m_npActCoeff[2] * I2 * IionicMolality;
    return pow(10.0 , l10actCoeff);
}


/**
 *  _osmoticCoeffHelgesonFixedForm()
 *
 *      Formula for the osmotic coefficient that occurs in
 *      the GWB. It is originally from Helgeson for a variable
 *      NaCl brine. It's to be used with extreme caution.
 */
double DebyeHuckel::
_osmoticCoeffHelgesonFixedForm() const
{
    const double a0 = 1.454;
    const double b0 = 0.02236;
    const double c0 = 9.380E-3;
    const double d0 = -5.362E-4;
    double Is = m_IionicMolalityStoich;
    if (Is <= 0.0) {
        return 0.0;
    }
    double Is2 = Is * Is;
    double bhat = 1.0 + a0 * sqrt(Is);
    double func = bhat - 2.0 * log(bhat) - 1.0/bhat;
    double v1 = m_A_Debye / (a0 * a0 * a0 * Is) * func;
    double oc = 1.0 - v1 + b0 * Is / 2.0 + 2.0 * c0 * Is2 / 3.0
                + 3.0 * d0 * Is2 * Is / 4.0;
    return oc;
}


/*
 *  _activityWaterHelgesonFixedForm()
 *
 *      Formula for the log of the activity of the water
 *      solvent that occurs in
 *      the GWB. It is originally from Helgeson for a variable
 *      NaCl brine. It's to be used with extreme caution.
 */
double DebyeHuckel::
_lnactivityWaterHelgesonFixedForm() const
{
    /*
     * Update the internally stored vector of molalities
     */
    calcMolalities();
    double oc = _osmoticCoeffHelgesonFixedForm();
    double sum = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        if (k != m_indexSolvent) {
            sum += std::max(m_molalities[k], 0.0);
        }
    }
    if (sum > 2.0 * m_maxIionicStrength) {
        sum = 2.0 *  m_maxIionicStrength;
    };
    double lac = - m_Mnaught * sum * oc;
    return lac;
}

/*
 * s_update_lnMolalityActCoeff():
 *
 *   Using internally stored values, this function calculates
 *   the activity coefficients for all species.
 *
 *   The ln(activity_solvent) is first calculated for the
 *   solvent. Then the molar based activity coefficient
 *   is calculated and returned.
 *
 *   ( Note this is the main routine for implementing the
 *     activity coefficient formulation.)
 */
void DebyeHuckel::s_update_lnMolalityActCoeff() const
{
    double z_k, zs_k1, zs_k2;
    /*
     * Update the internally stored vector of molalities
     */
    calcMolalities();
    /*
     * Calculate the apparent (real) ionic strength.
     *
     * Note this is not the stoichiometric ionic strengh,
     * where reactions of ions forming neutral salts
     * are ignorred in calculating the ionic strength.
     */
    m_IionicMolality = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        z_k = m_speciesCharge[k];
        m_IionicMolality += m_molalities[k] * z_k * z_k;
    }
    m_IionicMolality /= 2.0;

    if (m_IionicMolality > m_maxIionicStrength) {
        m_IionicMolality = m_maxIionicStrength;
    }

    /*
     * Calculate the stoichiometric ionic charge
     */
    m_IionicMolalityStoich = 0.0;
    for (size_t k = 0; k < m_kk; k++) {
        z_k = m_speciesCharge[k];
        zs_k1 =  m_speciesCharge_Stoich[k];
        if (z_k == zs_k1) {
            m_IionicMolalityStoich += m_molalities[k] * z_k * z_k;
        } else {
            zs_k2 = z_k - zs_k1;
            m_IionicMolalityStoich
            += m_molalities[k] * (zs_k1 * zs_k1 + zs_k2 * zs_k2);
        }
    }
    m_IionicMolalityStoich /= 2.0;

    if (m_IionicMolalityStoich > m_maxIionicStrength) {
        m_IionicMolalityStoich = m_maxIionicStrength;
    }

    /*
     * Possibly update the stored value of the
     * Debye-Huckel parameter A_Debye
     * This parameter appears on the top of the activity
     * coefficient expression.
     * It depends on T (and P), as it depends explicitly
     * on the temperature. Also, the dielectric constant
     * is usually a fairly strong function of T, also.
     */
    m_A_Debye = A_Debye_TP();

    /*
     * Calculate a safe value for the mole fraction
     * of the solvent
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);

    int est;
    double ac_nonPolar = 1.0;
    double numTmp = m_A_Debye * sqrt(m_IionicMolality);
    double denomTmp = m_B_Debye * sqrt(m_IionicMolality);
    double coeff;
    double lnActivitySolvent = 0.0;
    double tmp;
    double tmpLn;
    double y, yp1, sigma;
    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_lnActCoeffMolal[k] = - z_k * z_k * numTmp;
        }
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 / 3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality);
        break;

    case DHFORM_BDOT_AK:
        ac_nonPolar = _nonpolarActCoeff(m_IionicMolality);
        for (size_t k = 0; k < m_kk; k++) {
            est = m_electrolyteSpeciesType[k];
            if (est == cEST_nonpolarNeutral) {
                m_lnActCoeffMolal[k] = log(ac_nonPolar);
            } else {
                z_k = m_speciesCharge[k];
                m_lnActCoeffMolal[k] =
                    - z_k * z_k * numTmp / (1.0 + denomTmp * m_Aionic[k])
                    + log(10.0) * m_B_Dot[k] * m_IionicMolality;
            }
        }

        lnActivitySolvent = (xmolSolvent - 1.0)/xmolSolvent;
        coeff = 2.0 / 3.0 * m_A_Debye * m_Mnaught
                * sqrt(m_IionicMolality);
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                if (k != m_indexSolvent || m_Aionic[k] != 0.0) {
                    y = denomTmp * m_Aionic[k];
                    yp1 = y + 1.0;
                    sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                    z_k = m_speciesCharge[k];
                    tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
                }
            }
        }
        lnActivitySolvent += coeff * tmp;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            if ((k != m_indexSolvent) && (z_k != 0.0)) {
                tmp +=  m_B_Dot[k] * m_molalities[k];
            }
        }
        lnActivitySolvent -=
            m_Mnaught * log(10.0) * m_IionicMolality * tmp / 2.0;

        /*
         * Special section to implement the Helgeson fixed form
         * for the water brine activity coefficient.
         */
        if (m_useHelgesonFixedForm) {
            lnActivitySolvent = _lnactivityWaterHelgesonFixedForm();
        }
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_lnActCoeffMolal[k] =
                - z_k * z_k * numTmp / (1.0 + denomTmp)
                + log(10.0) * m_B_Dot[k] * m_IionicMolality;
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 /3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality) * sigma;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            if ((k != m_indexSolvent) && (z_k != 0.0)) {
                tmp +=  m_B_Dot[k] * m_molalities[k];
            }
        }
        lnActivitySolvent -=
            m_Mnaught * log(10.0) * m_IionicMolality * tmp / 2.0;

        break;

    case DHFORM_BETAIJ:
        denomTmp = m_B_Debye * m_Aionic[0];
        denomTmp *= sqrt(m_IionicMolality);
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent;

        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_lnActCoeffMolal[k] =
                    - z_k * z_k * numTmp / (1.0 + denomTmp);
                for (size_t j = 0; j < m_kk; j++) {
                    double beta =   m_Beta_ij.value(k, j);
#ifdef DEBUG_HKM_NOT
                    if (beta != 0.0) {
                        printf("b: k = %d, j = %d, betakj = %g\n",
                               k, j, beta);
                    }
#endif
                    m_lnActCoeffMolal[k] += 2.0 * m_molalities[j] * beta;
                }
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 -2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 /3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality) * sigma;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            for (size_t j = 0; j < m_kk; j++) {
                tmp +=
                    m_Beta_ij.value(k, j) * m_molalities[k] * m_molalities[j];
            }
        }
        lnActivitySolvent -= m_Mnaught * tmp;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp = m_B_Debye * sqrt(m_IionicMolality);
        denomTmp *= m_Aionic[0];
        numTmp = m_A_Debye * sqrt(m_IionicMolality);
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_lnActCoeffMolal[k] =
                    - z_k * z_k * numTmp / 3.0 / (1.0 + denomTmp);
                m_lnActCoeffMolal[k] +=
                    - 2.0 * z_k * z_k * m_A_Debye * tmpLn /
                    (3.0 * m_B_Debye * m_Aionic[0]);
                for (size_t j = 0; j < m_kk; j++) {
                    m_lnActCoeffMolal[k] += 2.0 * m_molalities[j] *
                                            m_Beta_ij.value(k, j);
                }
            }
        }
        sigma = 1.0 / (1.0 + denomTmp);
        lnActivitySolvent =
            (xmolSolvent - 1.0)/xmolSolvent +
            2.0 /3.0 * m_A_Debye * m_Mnaught *
            m_IionicMolality * sqrt(m_IionicMolality) * sigma;
        tmp = 0.0;
        for (size_t k = 0; k < m_kk; k++) {
            for (size_t j = 0; j < m_kk; j++) {
                tmp +=
                    m_Beta_ij.value(k, j) * m_molalities[k] * m_molalities[j];
            }
        }
        lnActivitySolvent -= m_Mnaught * tmp;
        break;

    default:
        printf("ERROR\n");
        exit(EXIT_FAILURE);
    }
    /*
     * Above, we calculated the ln(activitySolvent). Translate that
     * into the molar-based activity coefficient by dividing by
     * the solvent mole fraction. Solvents are not on the molality
     * scale.
     */
    xmolSolvent = moleFraction(m_indexSolvent);
    m_lnActCoeffMolal[m_indexSolvent] =
        lnActivitySolvent - log(xmolSolvent);
}

/*
 * s_update_dMolalityActCoeff_dT()         (private, const )
 *
 *   Using internally stored values, this function calculates
 *   the temperature derivative of the logarithm of the
 *   activity coefficient for all species in the mechanism.
 *
 *   We assume that the activity coefficients are current.
 *
 *   solvent activity coefficient is on the molality
 *   scale. Its derivative is too.
 */
void DebyeHuckel::s_update_dlnMolalityActCoeff_dT() const
{
    double z_k, coeff, tmp, y, yp1, sigma, tmpLn;
    // First we store dAdT explicitly here
    double dAdT =  dA_DebyedT_TP();
    if (dAdT == 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            m_dlnActCoeffMolaldT[k] = 0.0;
        }
        return;
    }
    /*
     * Calculate a safe value for the mole fraction
     * of the solvent
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);


    double sqrtI  = sqrt(m_IionicMolality);
    double numdAdTTmp = dAdT * sqrtI;
    double denomTmp = m_B_Debye * sqrtI;
    double d_lnActivitySolvent_dT = 0;

    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 1; k < m_kk; k++) {
            m_dlnActCoeffMolaldT[k] =
                m_lnActCoeffMolal[k] * dAdT / m_A_Debye;
        }
        d_lnActivitySolvent_dT = 2.0 / 3.0 * dAdT * m_Mnaught *
                                 m_IionicMolality * sqrt(m_IionicMolality);
        m_dlnActCoeffMolaldT[m_indexSolvent] =  d_lnActivitySolvent_dT;
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] =
                - z_k * z_k * numdAdTTmp / (1.0 + denomTmp * m_Aionic[k]);
        }

        m_dlnActCoeffMolaldT[m_indexSolvent] = 0.0;

        coeff = 2.0 / 3.0 * dAdT * m_Mnaught * sqrtI;
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                y = denomTmp * m_Aionic[k];
                yp1 = y + 1.0;
                sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                z_k = m_speciesCharge[k];
                tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
            }
        }
        m_dlnActCoeffMolaldT[m_indexSolvent] += coeff * tmp;
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldT[k] =
                - z_k * z_k * numdAdTTmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldT[m_indexSolvent] =
            2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldT[k] =
                    - z_k * z_k * numdAdTTmp / (1.0 + denomTmp);
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldT[m_indexSolvent] =
            2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldT[k] =
                    - z_k * z_k * numdAdTTmp / (1.0 + denomTmp)
                    - 2.0 * z_k * z_k * dAdT * tmpLn
                    / (m_B_Debye * m_Aionic[0]);
                m_dlnActCoeffMolaldT[k] /= 3.0;
            }
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_dlnActCoeffMolaldT[m_indexSolvent] =
            2.0 /3.0 * dAdT * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        printf("ERROR\n");
        exit(EXIT_FAILURE);
        break;
    }


}

/*
 * s_update_d2lnMolalityActCoeff_dT2()         (private, const )
 *
 *   Using internally stored values, this function calculates
 *   the temperature 2nd derivative of the logarithm of the
 *   activity coefficient
 *   for all species in the mechanism.
 *
 *   We assume that the activity coefficients are current.
 *
 *   solvent activity coefficient is on the molality
 *   scale. Its derivatives are too.
 */
void DebyeHuckel::s_update_d2lnMolalityActCoeff_dT2() const
{
    double z_k, coeff, tmp, y, yp1, sigma, tmpLn;
    double dAdT =  dA_DebyedT_TP();
    double d2AdT2 = d2A_DebyedT2_TP();
    if (d2AdT2 == 0.0 && dAdT == 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            m_d2lnActCoeffMolaldT2[k] = 0.0;
        }
        return;
    }

    /*
     * Calculate a safe value for the mole fraction
     * of the solvent
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);


    double sqrtI  = sqrt(m_IionicMolality);
    double numd2AdT2Tmp = d2AdT2 * sqrtI;
    double denomTmp = m_B_Debye * sqrtI;

    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 0; k < m_kk; k++) {
            m_d2lnActCoeffMolaldT2[k] =
                m_lnActCoeffMolal[k] * d2AdT2 / m_A_Debye;
        }
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_d2lnActCoeffMolaldT2[k] =
                - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp * m_Aionic[k]);
        }

        m_d2lnActCoeffMolaldT2[m_indexSolvent] = 0.0;

        coeff = 2.0 / 3.0 * d2AdT2 * m_Mnaught * sqrtI;
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                y = denomTmp * m_Aionic[k];
                yp1 = y + 1.0;
                sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                z_k = m_speciesCharge[k];
                tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
            }
        }
        m_d2lnActCoeffMolaldT2[m_indexSolvent] += coeff * tmp;
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_d2lnActCoeffMolaldT2[k] =
                - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_d2lnActCoeffMolaldT2[m_indexSolvent] =
            2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_d2lnActCoeffMolaldT2[k] =
                    - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp);
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 -2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_d2lnActCoeffMolaldT2[m_indexSolvent] =
            2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_d2lnActCoeffMolaldT2[k] =
                    - z_k * z_k * numd2AdT2Tmp / (1.0 + denomTmp)
                    - 2.0 * z_k * z_k * d2AdT2 * tmpLn
                    / (m_B_Debye * m_Aionic[0]);
                m_d2lnActCoeffMolaldT2[k] /= 3.0;
            }
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_d2lnActCoeffMolaldT2[m_indexSolvent] =
            2.0 /3.0 * d2AdT2 * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        printf("ERROR\n");
        exit(EXIT_FAILURE);
        break;
    }
}

/*
 * s_update_dlnMolalityActCoeff_dP()         (private, const )
 *
 *   Using internally stored values, this function calculates
 *   the pressure derivative of the logarithm of the
 *   activity coefficient for all species in the mechanism.
 *
 *   We assume that the activity coefficients, molalities,
 *   and A_Debye are current.
 *
 *   solvent activity coefficient is on the molality
 *   scale. Its derivatives are too.
 */
void DebyeHuckel::s_update_dlnMolalityActCoeff_dP() const
{
    double z_k, coeff, tmp, y, yp1, sigma, tmpLn;
    int est;
    double dAdP =  dA_DebyedP_TP();
    if (dAdP == 0.0) {
        for (size_t k = 0; k < m_kk; k++) {
            m_dlnActCoeffMolaldP[k] = 0.0;
        }
        return;
    }
    /*
     * Calculate a safe value for the mole fraction
     * of the solvent
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    xmolSolvent = std::max(8.689E-3, xmolSolvent);


    double sqrtI  = sqrt(m_IionicMolality);
    double numdAdPTmp = dAdP * sqrtI;
    double denomTmp = m_B_Debye * sqrtI;

    switch (m_formDH) {
    case DHFORM_DILUTE_LIMIT:
        for (size_t k = 0; k < m_kk; k++) {
            m_dlnActCoeffMolaldP[k] =
                m_lnActCoeffMolal[k] * dAdP / m_A_Debye;
        }
        break;

    case DHFORM_BDOT_AK:
        for (size_t k = 0; k < m_kk; k++) {
            est = m_electrolyteSpeciesType[k];
            if (est == cEST_nonpolarNeutral) {
                m_lnActCoeffMolal[k] = 0.0;
            } else {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldP[k] =
                    - z_k * z_k * numdAdPTmp / (1.0 + denomTmp * m_Aionic[k]);
            }
        }

        m_dlnActCoeffMolaldP[m_indexSolvent] = 0.0;

        coeff = 2.0 / 3.0 * dAdP * m_Mnaught * sqrtI;
        tmp = 0.0;
        if (denomTmp > 0.0) {
            for (size_t k = 0; k < m_kk; k++) {
                y = denomTmp * m_Aionic[k];
                yp1 = y + 1.0;
                sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
                z_k = m_speciesCharge[k];
                tmp += m_molalities[k] * z_k * z_k * sigma / 2.0;
            }
        }
        m_dlnActCoeffMolaldP[m_indexSolvent] += coeff * tmp;
        break;

    case DHFORM_BDOT_ACOMMON:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            z_k = m_speciesCharge[k];
            m_dlnActCoeffMolaldP[k] =
                - z_k * z_k * numdAdPTmp / (1.0 + denomTmp);
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldP[m_indexSolvent] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_BETAIJ:
        denomTmp *= m_Aionic[0];
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldP[k] =
                    - z_k * z_k * numdAdPTmp / (1.0 + denomTmp);
            }
        }
        if (denomTmp > 0.0) {
            y = denomTmp;
            yp1 = y + 1.0;
            sigma = 3.0 / (y * y * y) * (yp1 - 1.0/yp1 - 2.0*log(yp1));
        } else {
            sigma = 0.0;
        }
        m_dlnActCoeffMolaldP[m_indexSolvent] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    case DHFORM_PITZER_BETAIJ:
        denomTmp *= m_Aionic[0];
        tmpLn = log(1.0 + denomTmp);
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                z_k = m_speciesCharge[k];
                m_dlnActCoeffMolaldP[k] =
                    - z_k * z_k * numdAdPTmp / (1.0 + denomTmp)
                    - 2.0 * z_k * z_k * dAdP * tmpLn
                    / (m_B_Debye * m_Aionic[0]);
                m_dlnActCoeffMolaldP[k] /= 3.0;
            }
        }

        sigma = 1.0 / (1.0 + denomTmp);
        m_dlnActCoeffMolaldP[m_indexSolvent] =
            2.0 /3.0 * dAdP * m_Mnaught *
            m_IionicMolality * sqrtI * sigma;
        break;

    default:
        printf("ERROR\n");
        exit(EXIT_FAILURE);
        break;
    }
}

/*
 * Updates the standard state thermodynamic functions at the current T and P of the solution.
 *
 * @internal
 *
 * This function gets called for every call to functions in this
 * class. It checks to see whether the temperature or pressure has changed and
 * thus the ss thermodynamics functions for all of the species
 * must be recalculated.
 */
//  void DebyeHuckel::_updateStandardStateThermo() const {
// doublereal tnow = temperature();
// doublereal pnow = m_Pcurrent;
// if (m_waterSS) {
//   m_waterSS->setTempPressure(tnow, pnow);
// }
// m_VPSS_ptr->setState_TP(tnow, pnow);
// VPStandardStateTP::updateStandardStateThermo();

//}

}
