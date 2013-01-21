/**
 *  @file IdealMolalSoln.cpp
 *   ThermoPhase object for the ideal molal equation of
 * state (see \ref thermoprops
 * and class \link Cantera::IdealMolalSoln IdealMolalSoln\endlink).
 *
 * Definition file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon
 * activities on the molality scale. The Ideal molal
 * solution assumes that all molality-based activity
 * coefficients are equal to one. This turns out, actually, to be
 * highly nonlinear when the solvent densities get low.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/thermo/ThermoFactory.h"

#include <cmath>
#include <fstream>

using namespace ctml;

namespace Cantera
{

/*
 * Default constructor
 */
IdealMolalSoln::IdealMolalSoln() :
    MolalityVPSSTP(),
    m_formGC(2),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.20),
    IMS_gamma_o_min_(0.00001),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(.05),
    IMS_slopefCut_(0.6),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_slopegCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0)
{
}

/*
 * Copy Constructor:
 *
 *  Note this stuff will not work until the underlying phase
 *  has a working copy constructor
 */
IdealMolalSoln::IdealMolalSoln(const IdealMolalSoln& b) :
    MolalityVPSSTP(b)
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
IdealMolalSoln& IdealMolalSoln::
operator=(const IdealMolalSoln& b)
{
    if (&b != this) {
        MolalityVPSSTP::operator=(b);
        m_speciesMolarVolume  = b.m_speciesMolarVolume;
        m_formGC              = b.m_formGC;
        IMS_typeCutoff_       = b.IMS_typeCutoff_;
        IMS_X_o_cutoff_       = b.IMS_X_o_cutoff_;
        IMS_gamma_o_min_      = b.IMS_gamma_o_min_;
        IMS_gamma_k_min_      = b.IMS_gamma_k_min_;
        IMS_cCut_             = b.IMS_cCut_;
        IMS_slopefCut_        = b.IMS_slopefCut_;
        IMS_dfCut_            = b.IMS_dfCut_;
        IMS_efCut_            = b.IMS_efCut_;
        IMS_afCut_            = b.IMS_afCut_;
        IMS_bfCut_            = b.IMS_bfCut_;
        IMS_slopegCut_        = b.IMS_slopegCut_;
        IMS_dgCut_            = b.IMS_dgCut_;
        IMS_egCut_            = b.IMS_egCut_;
        IMS_agCut_            = b.IMS_agCut_;
        IMS_bgCut_            = b.IMS_bgCut_;
        m_pp                  = b.m_pp;
        m_tmpV                = b.m_tmpV;
        IMS_lnActCoeffMolal_  = b.IMS_lnActCoeffMolal_;
    }
    return *this;
}

IdealMolalSoln::IdealMolalSoln(const std::string& inputFile,
                               const std::string& id) :
    MolalityVPSSTP(),
    m_formGC(2),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(0.00001),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(.05),
    IMS_slopefCut_(0.6),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_slopegCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0)
{
    initThermoFile(inputFile, id);
}

IdealMolalSoln::IdealMolalSoln(XML_Node& root, const std::string& id) :
    MolalityVPSSTP(),
    m_formGC(2),
    IMS_typeCutoff_(0),
    IMS_X_o_cutoff_(0.2),
    IMS_gamma_o_min_(0.00001),
    IMS_gamma_k_min_(10.0),
    IMS_cCut_(.05),
    IMS_slopefCut_(0.6),
    IMS_dfCut_(0.0),
    IMS_efCut_(0.0),
    IMS_afCut_(0.0),
    IMS_bfCut_(0.0),
    IMS_slopegCut_(0.0),
    IMS_dgCut_(0.0),
    IMS_egCut_(0.0),
    IMS_agCut_(0.0),
    IMS_bgCut_(0.0)
{
    importPhase(*findXMLPhase(&root, id), this);
}

/*
 *
 * ~IdealMolalSoln():   (virtual)
 *
 * Destructor: does nothing:
 *
 */
IdealMolalSoln::~IdealMolalSoln()
{
}

/**
 *
 */
thermo_t* IdealMolalSoln::duplMyselfAsThermoPhase() const
{
    return new IdealMolalSoln(*this);
}

//
// -------- Molar Thermodynamic Properties of the Solution ---------------
//
/*
 * Molar enthalpy of the solution: Units: J/kmol.
 *
 * Returns the amount of enthalpy per mole of solution.
 * For an ideal molal solution,
 * \f[
 * \bar{h}(T, P, X_k) = \sum_k X_k \bar{h}_k(T)
 * \f]
 * The formula is written in terms of the partial molar enthalpies.
 * \f$ \bar{h}_k(T, p, m_k) \f$.
 * See the partial molar enthalpy function, getPartialMolarEnthalpies(),
 * for details.
 *
 * Units: J/kmol
 */
doublereal IdealMolalSoln::enthalpy_mole() const
{
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    getMoleFractions(DATA_PTR(m_pp));
    double val = mean_X(DATA_PTR(m_tmpV));
    return val;
}

/*
 * Molar internal energy of the solution: Units: J/kmol.
 *
 * Returns the amount of internal energy per mole of solution.
 * For an ideal molal solution,
 * \f[
 * \bar{u}(T, P, X_k) = \sum_k X_k \bar{u}_k(T)
 * \f]
 * The formula is written in terms of the partial molar internal energy.
 * \f$ \bar{u}_k(T, p, m_k) \f$.
 */
doublereal IdealMolalSoln::intEnergy_mole() const
{
    getPartialMolarEnthalpies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

/*
 * Molar entropy of the solution: Units J/kmol/K.
 *
 * Returns the amount of entropy per mole of solution.
 * For an ideal molal solution,
 * \f[
 * \bar{s}(T, P, X_k) = \sum_k X_k \bar{s}_k(T)
 * \f]
 * The formula is written in terms of the partial molar entropies.
 * \f$ \bar{s}_k(T, p, m_k) \f$.
 * See the partial molar entropies function, getPartialMolarEntropies(),
 * for details.
 *
 * Units: J/kmol/K.
 */
doublereal IdealMolalSoln::entropy_mole() const
{
    getPartialMolarEntropies(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

/*
 * Molar Gibbs function for the solution: Units J/kmol.
 *
 * Returns the gibbs free energy of the solution per mole
 * of the solution.
 *
 * \f[
 * \bar{g}(T, P, X_k) = \sum_k X_k \mu_k(T)
 * \f]
 *
 * Units: J/kmol
 */
doublereal IdealMolalSoln::gibbs_mole() const
{
    getChemPotentials(DATA_PTR(m_tmpV));
    return mean_X(DATA_PTR(m_tmpV));
}

/*
 * Molar heat capacity at constant pressure: Units: J/kmol/K.
 *  * \f[
 * \bar{c}_p(T, P, X_k) = \sum_k X_k \bar{c}_{p,k}(T)
 * \f]
 *
 * Units: J/kmol/K
 */
doublereal IdealMolalSoln::cp_mole() const
{
    getPartialMolarCp(DATA_PTR(m_tmpV));
    double val = mean_X(DATA_PTR(m_tmpV));
    return val;
}

/*
 * Molar heat capacity at constant volume: Units: J/kmol/K.
 * NOT IMPLEMENTED.
 * Units: J/kmol/K
 */
doublereal IdealMolalSoln::cv_mole() const
{
    return err("not implemented");
}

//
// ------- Mechanical Equation of State Properties ------------------------
//



/*
 * Set the pressure at constant temperature. Units: Pa.
 * This method sets a constant within the object.
 * The mass density is not a function of pressure.
 */
void IdealMolalSoln::setPressure(doublereal p)
{
    setState_TP(temperature(), p);
}

void IdealMolalSoln::calcDensity()
{
    double* vbar = &m_pp[0];
    getPartialMolarVolumes(vbar);
    double* x = &m_tmpV[0];
    getMoleFractions(x);
    doublereal vtotal = 0.0;
    for (size_t i = 0; i < m_kk; i++) {
        vtotal += vbar[i] * x[i];
    }
    doublereal dd = meanMolecularWeight() / vtotal;
    phase_t::setDensity(dd);
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
doublereal IdealMolalSoln::isothermalCompressibility() const
{
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
doublereal IdealMolalSoln::thermalExpansionCoeff() const
{
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
void IdealMolalSoln::setDensity(const doublereal rho)
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
void IdealMolalSoln::setMolarDensity(const doublereal conc)
{
    double concI = phase_t::molarDensity();
    if (conc != concI) {
        throw CanteraError("IdealMolalSoln::setMolarDensity",
                           "molarDensity/denisty is not an independent variable");
    }
}

void IdealMolalSoln::setState_TP(doublereal temp, doublereal pres)
{
    phase_t::setTemperature(temp);
    m_Pcurrent = pres;
    updateStandardStateThermo();
    //m_densWaterSS = m_waterSS->density();
    calcDensity();
}

//
// ------- Activities and Activity Concentrations
//

/*
 * This method returns an array of activity concentrations \f$ C^a_k\f$.
 * \f$ C^a_k\f$ are defined such that
 * \f$ a_k = C^a_k / C^s_k, \f$ where \f$ C^s_k \f$
 * is a standard concentration
 * defined below.  These activity concentrations are used
 * by kinetics manager classes to compute the forward and
 * reverse rates of elementary reactions.
 *
 * @param c Array of activity concentrations. The
 *           units depend upon the implementation of the
 *           reaction rate expressions within the phase.
 */
void IdealMolalSoln::getActivityConcentrations(doublereal* c) const
{
    if (m_formGC != 1) {
        double c_solvent = standardConcentration();
        getActivities(c);
        for (size_t k = 0; k < m_kk; k++) {
            c[k] *= c_solvent;
        }
    } else {
        getActivities(c);
        for (size_t k = 0; k < m_kk; k++) {
            double c0 = standardConcentration(k);
            c[k] *= c0;
        }
    }
}

/*
 * The standard concentration \f$ C^s_k \f$ used to normalize
 * the activity concentration. In many cases, this quantity
 * will be the same for all species in a phase - for example,
 * for an ideal gas \f$ C^s_k = P/\hat R T \f$. For this
 * reason, this method returns a single value, instead of an
 * array.  However, for phases in which the standard
 * concentration is species-specific (e.g. surface species of
 * different sizes), this method may be called with an
 * optional parameter indicating the species.
 *
 */
doublereal IdealMolalSoln::standardConcentration(size_t k) const
{
    double c0 = 1.0, mvSolvent;
    switch (m_formGC) {
    case 0:
        break;
    case 1:
        c0 = 1.0 /m_speciesMolarVolume[m_indexSolvent];
        break;
    case 2:
        mvSolvent = m_speciesMolarVolume[m_indexSolvent];
        c0 = 1.0 / mvSolvent;
        break;
    }
    return c0;
}

/*
 * Returns the natural logarithm of the standard
 * concentration of the kth species
 */
doublereal IdealMolalSoln::logStandardConc(size_t k) const
{
    double c0 = standardConcentration(k);
    return log(c0);
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
void IdealMolalSoln::getUnitsStandardConc(double* uA, int k, int sizeUA) const
{
    int eos = eosType();
    if (eos == 0) {
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
 * Get the array of non-dimensional molality-based
 * activities at the current solution temperature,
 * pressure, and solution concentration.
 *
 *  The max against xmolSolventMIN is to limit the activity
 *  coefficient to be finite as the solvent mf goes to zero.
 */
void IdealMolalSoln::getActivities(doublereal* ac) const
{
    _updateStandardStateThermo();
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    if (IMS_typeCutoff_ == 0) {
        calcMolalities();
        for (size_t k = 0; k < m_kk; k++) {
            ac[k] = m_molalities[k];
        }
        double xmolSolvent = moleFraction(m_indexSolvent);
        xmolSolvent = std::max(m_xmolSolventMIN, xmolSolvent);
        ac[m_indexSolvent] =
            exp((xmolSolvent - 1.0)/xmolSolvent);
    } else {

        s_updateIMS_lnMolalityActCoeff();
        /*
         * Now calculate the array of activities.
         */
        for (size_t k = 1; k < m_kk; k++) {
            ac[k] = m_molalities[k] * exp(IMS_lnActCoeffMolal_[k]);
        }
        double xmolSolvent = moleFraction(m_indexSolvent);
        ac[m_indexSolvent] =
            exp(IMS_lnActCoeffMolal_[m_indexSolvent]) * xmolSolvent;

    }
}

/*
 * Get the array of non-dimensional Molality based
 * activity coefficients at
 * the current solution temperature, pressure, and
 * solution concentration.
 * See Denbigh
 * (note solvent activity coefficient is on the molar scale).
 *
 *  The max against xmolSolventMIN is to limit the activity
 *  coefficient to be finite as the solvent mf goes to zero.
 */
void IdealMolalSoln::
getMolalityActivityCoefficients(doublereal* acMolality) const
{
    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 0; k < m_kk; k++) {
            acMolality[k] = 1.0;
        }
        double xmolSolvent = moleFraction(m_indexSolvent);
        xmolSolvent = std::max(m_xmolSolventMIN, xmolSolvent);
        acMolality[m_indexSolvent] =
            exp((xmolSolvent - 1.0)/xmolSolvent) / xmolSolvent;
    } else {
        s_updateIMS_lnMolalityActCoeff();
        std::copy(IMS_lnActCoeffMolal_.begin(), IMS_lnActCoeffMolal_.end(), acMolality);
        for (size_t k = 0; k < m_kk; k++) {
            acMolality[k] = exp(acMolality[k]);
        }
    }
}

//
// ------ Partial Molar Properties of the Solution -----------------
//

/*
 * Get the species chemical potentials: Units: J/kmol.
 *
 * This function returns a vector of chemical potentials of the
 * species in solution.
 *
 * \f[
 *    \mu_k = \mu^{o}_k(T,P) + R T \ln(\frac{m_k}{m^\Delta})
 * \f]
 * \f[
 *    \mu_w = \mu^{o}_w(T,P) +
 *            R T ((X_w - 1.0) / X_w)
 * \f]
 *
 * \f$ w \f$ refers to the solvent species.
 * \f$ X_w \f$ is the mole fraction of the solvent.
 * \f$ m_k \f$ is the molality of the kth solute.
 * \f$ m^\Delta is 1 gmol solute per kg solvent. \f$
 *
 * Units: J/kmol.
 */
void IdealMolalSoln::getChemPotentials(doublereal* mu) const
{
    double xx;

    // Assertion is made for speed
    AssertThrow(m_indexSolvent == 0, "solvent not the first species");

    /*
     * First get the standard chemical potentials
     *  -> this requires updates of standard state as a function
     *     of T and P
     * These are defined at unit molality.
     */
    getStandardChemPotentials(mu);
    /*
     * Update the molality array, m_molalities()
     *   This requires an update due to mole fractions
     */
    calcMolalities();
    /*
     * get the solvent mole fraction
     */
    double xmolSolvent = moleFraction(m_indexSolvent);
    doublereal RT = GasConstant * temperature();

    if (IMS_typeCutoff_ == 0 || xmolSolvent > 3.* IMS_X_o_cutoff_/2.0) {

        for (size_t k = 1; k < m_kk; k++) {
            xx = std::max(m_molalities[k], SmallNumber);
            mu[k] += RT * log(xx);
        }
        /*
         * Do the solvent
         *  -> see my notes
         */

        xx = std::max(xmolSolvent, SmallNumber);
        mu[m_indexSolvent] +=
            (RT * (xmolSolvent - 1.0) / xx);
    } else {
        /*
         * Update the activity coefficients
         * This also updates the internal molality array.
         */
        s_updateIMS_lnMolalityActCoeff();


        for (size_t k = 1; k < m_kk; k++) {
            xx = std::max(m_molalities[k], SmallNumber);
            mu[k] += RT * (log(xx) + IMS_lnActCoeffMolal_[k]);
        }
        xx = std::max(xmolSolvent, SmallNumber);
        mu[m_indexSolvent] +=
            RT * (log(xx) + IMS_lnActCoeffMolal_[m_indexSolvent]);
    }

}

/*
 * Returns an array of partial molar enthalpies for the species
 * in the mixture: Units (J/kmol).
 *
 */
void IdealMolalSoln::getPartialMolarEnthalpies(doublereal* hbar) const
{
    getEnthalpy_RT(hbar);
    doublereal RT = _RT();
    for (size_t k = 0; k < m_kk; k++) {
        hbar[k] *= RT;
    }
}

/*
 * Returns an array of partial molar entropies of the species in the
 * solution: Units: J/kmol.
 *
 * Maxwell's equations provide an insight in how to calculate this
 * (p.215 Smith and Van Ness)
 * \f[
 *      \frac{d(\mu_k)}{dT} = -\bar{s}_i
 * \f]
 * For this phase, the partial molar entropies are equal to the
 * standard state species entropies plus the ideal molal solution contribution.
 *
 * \f[
 *   \bar{s}_k(T,P) =  s^0_k(T) - R log( m_k )
 * \f]
 * \f[
 *   \bar{s}_w(T,P) =  s^0_w(T) - R ((X_w - 1.0) / X_w)
 * \f]
 *
 * The subscript, w, refers to the solvent species. \f$ X_w \f$ is
 * the mole fraction of solvent.
 * The reference-state pure-species entropies,\f$ s^0_k(T) \f$,
 * at the reference pressure, \f$ P_{ref} \f$, are computed by the
 * species thermodynamic
 * property manager. They are polynomial functions of temperature.
 * @see SpeciesThermo
 */
void IdealMolalSoln::getPartialMolarEntropies(doublereal* sbar) const
{
    getEntropy_R(sbar);
    doublereal R = GasConstant;
    doublereal mm;
    calcMolalities();
    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                mm = std::max(SmallNumber, m_molalities[k]);
                sbar[k] -= R * log(mm);
            }
        }
        double xmolSolvent = moleFraction(m_indexSolvent);
        sbar[m_indexSolvent] -= (R * (xmolSolvent - 1.0) / xmolSolvent);
    } else {
        /*
         * Update the activity coefficients, This also update the
         * internally stored molalities.
         */
        s_updateIMS_lnMolalityActCoeff();
        /*
         * First we will add in the obvious dependence on the T
         * term out front of the log activity term
         */
        doublereal mm;
        for (size_t k = 0; k < m_kk; k++) {
            if (k != m_indexSolvent) {
                mm = std::max(SmallNumber, m_molalities[k]);
                sbar[k] -= R * (log(mm) + IMS_lnActCoeffMolal_[k]);
            }
        }
        double xmolSolvent = moleFraction(m_indexSolvent);
        mm = std::max(SmallNumber, xmolSolvent);
        sbar[m_indexSolvent] -= R *(log(mm) + IMS_lnActCoeffMolal_[m_indexSolvent]);

    }
}

/*
 * Returns an array of partial molar volumes of the species
 * in the solution: Units: m^3 kmol-1.
 *
 * For this solution, the partial molar volumes are equal to the
 * constant species molar volumes.
 *
 * Units: m^3 kmol-1.
 */
void IdealMolalSoln::getPartialMolarVolumes(doublereal* vbar) const
{
    getStandardVolumes(vbar);
}

/*
 * Partial molar heat capacity of the solution: Units: J/kmol/K.
 *
 *   The kth partial molar heat capacity is equal to
 *   the temperature derivative of the partial molar
 *   enthalpy of the kth species in the solution at constant
 *   P and composition (p. 220 Smith and Van Ness).
 *  \f[
 *  \bar{Cp}_k(T,P) =  {Cp}^0_k(T)
 * \f]
 *
 *   For this solution, this is equal to the reference state
 *   heat capacities.
 *
 *  Units: J/kmol/K
 */
void IdealMolalSoln::getPartialMolarCp(doublereal* cpbar) const
{
    /*
     * Get the nondimensional gibbs standard state of the
     * species at the T and P of the solution.
     */
    getCp_R(cpbar);

    for (size_t k = 0; k < m_kk; k++) {
        cpbar[k] *= GasConstant;
    }
}

/*
 * -------- Properties of the Standard State of the Species
 *           in the Solution ------------------
 */



/*
 * ------ Thermodynamic Values for the Species Reference States ---
 */

// -> This is handled by VPStandardStatesTP

/*
 *  -------------- Utilities -------------------------------
 */

/*
 *  Initialization routine for an IdealMolalSoln phase.
 *
 * This is a virtual routine. This routine will call initThermo()
 * for the parent class as well.
 */
void IdealMolalSoln::initThermo()
{
    initLengths();
    MolalityVPSSTP::initThermo();
}

/*
 *   Import and initialize an IdealMolalSoln phase
 *   specification in an XML tree into the current object.
 *
 *   This routine is called from importPhase() to finish
 *   up the initialization of the thermo object. It reads in the
 *   species molar volumes.
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
void IdealMolalSoln::initThermoXML(XML_Node& phaseNode, const std::string& id)
{
    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("IdealMolalSoln::initThermoXML",
                           "no thermo XML node");
    }

    /*
     * Initialize the whole thermo object, using a virtual function.
     */
    initThermo();

    if (id.size() > 0) {
        std::string idp = phaseNode.id();
        if (idp != id) {
            throw CanteraError("IdealMolalSoln::initThermo",
                               "phasenode and Id are incompatible");
        }
    }

    /*
     * Find the Thermo XML node
     */
    if (!phaseNode.hasChild("thermo")) {
        throw CanteraError("IdealMolalSoln::initThermo",
                           "no thermo XML node");
    }
    XML_Node& thermoNode = phaseNode.child("thermo");

    /*
     * Possible change the form of the standard concentrations
     */
    if (thermoNode.hasChild("standardConc")) {
        XML_Node& scNode = thermoNode.child("standardConc");
        m_formGC = 2;
        std::string formString = scNode.attrib("model");
        if (formString != "") {
            if (formString == "unity") {
                m_formGC = 0;
            } else if (formString == "molar_volume") {
                m_formGC = 1;
            } else if (formString == "solvent_volume") {
                m_formGC = 2;
            } else {
                throw CanteraError("IdealMolalSoln::initThermo",
                                   "Unknown standardConc model: " + formString);
            }
        }
    }

    /*
     * Get the Name of the Solvent:
     *      <solvent> solventName </solvent>
     */
    std::string solventName = "";
    if (thermoNode.hasChild("solvent")) {
        XML_Node& scNode = thermoNode.child("solvent");
        std::vector<std::string> nameSolventa;
        getStringArray(scNode, nameSolventa);
        int nsp = static_cast<int>(nameSolventa.size());
        if (nsp != 1) {
            throw CanteraError("IdealMolalSoln::initThermoXML",
                               "badly formed solvent XML node");
        }
        solventName = nameSolventa[0];
    }

    if (thermoNode.hasChild("activityCoefficients")) {
        XML_Node& acNode = thermoNode.child("activityCoefficients");
        std::string modelString = acNode.attrib("model");
        IMS_typeCutoff_ = 0;
        if (modelString != "IdealMolalSoln") {
            throw CanteraError("IdealMolalSoln::initThermoXML",
                               "unknown ActivityCoefficient model: " + modelString);
        }
        if (acNode.hasChild("idealMolalSolnCutoff")) {
            XML_Node& ccNode = acNode.child("idealMolalSolnCutoff");
            modelString = ccNode.attrib("model");
            if (modelString != "") {
                if (modelString == "polyExp") {
                    IMS_typeCutoff_ = 2;
                } else if (modelString == "poly") {
                    IMS_typeCutoff_ = 1;
                } else {
                    throw CanteraError("IdealMolalSoln::initThermoXML",
                                       "Unknown  idealMolalSolnCutoff form: " + modelString);
                }

                if (ccNode.hasChild("gamma_o_limit")) {
                    IMS_gamma_o_min_ = getFloat(ccNode, "gamma_o_limit");
                }
                if (ccNode.hasChild("gamma_k_limit")) {
                    IMS_gamma_k_min_ = getFloat(ccNode, "gamma_k_limit");
                }
                if (ccNode.hasChild("X_o_cutoff")) {
                    IMS_X_o_cutoff_ = getFloat(ccNode, "X_o_cutoff");
                }
                if (ccNode.hasChild("c_0_param")) {
                    IMS_cCut_ = getFloat(ccNode, "c_0_param");
                }
                if (ccNode.hasChild("slope_f_limit")) {
                    IMS_slopefCut_ = getFloat(ccNode, "slope_f_limit");
                }
                if (ccNode.hasChild("slope_g_limit")) {
                    IMS_slopegCut_ = getFloat(ccNode, "slope_g_limit");
                }

            }
        }
    }


    /*
     * Reconcile the solvent name and index.
     */
    for (size_t k = 0; k < m_kk; k++) {
        std::string sname = speciesName(k);
        if (solventName == sname) {
            m_indexSolvent = k;
            break;
        }
    }
    if (m_indexSolvent == npos) {
        std::cout << "IdealMolalSoln::initThermo: Solvent Name not found"
                  << std::endl;
        throw CanteraError("IdealMolalSoln::initThermo",
                           "Solvent name not found");
    }
    if (m_indexSolvent != 0) {
        throw CanteraError("IdealMolalSoln::initThermo",
                           "Solvent " + solventName +
                           " should be first species");
    }


    /*
     * Now go get the molar volumes
     */
    XML_Node& speciesList = phaseNode.child("speciesArray");
    XML_Node* speciesDB =
        get_XML_NameID("speciesData", speciesList["datasrc"],
                       &phaseNode.root());
    const std::vector<std::string> &sss = speciesNames();

    for (size_t k = 0; k < m_kk; k++) {
        XML_Node* s =  speciesDB->findByAttr("name", sss[k]);
        XML_Node* ss = s->findByName("standardState");
        m_speciesMolarVolume[k] = getFloat(*ss, "molarVolume", "toSI");
    }

    IMS_typeCutoff_ = 2;
    if (IMS_typeCutoff_ == 2) {
        calcIMSCutoffParams_();
    }

    MolalityVPSSTP::initThermoXML(phaseNode, id);


    setMoleFSolventMin(1.0E-5);
    /*
     * Set the state
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
void IdealMolalSoln::setParameters(int n, doublereal* const c)
{
}

void IdealMolalSoln::getParameters(int& n, doublereal* const c) const
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
 * HKM -> Right now, the parameters are set elsewhere (initThermo)
 *        It just didn't seem to fit.
 */
void IdealMolalSoln::setParametersFromXML(const XML_Node& eosdata)
{
}

/*
 * ----------- Critical State Properties --------------------------
 */

/*
 * ------------ Private and Restricted Functions ------------------
 */

/*
 * Bail out of functions with an error exit if they are not
 * implemented.
 */
doublereal IdealMolalSoln::err(const std::string& msg) const
{
    throw CanteraError("IdealMolalSoln",
                       "Unfinished func called: " + msg);
    return 0.0;
}



// This function will be called to update the internally stored
// natural logarithm of the molality activity coefficients
/*
 * Normally they are all one. However, sometimes they are not,
 * due to stability schemes
 *
 *    gamma_k_molar =  gamma_k_molal / Xmol_solvent
 *
 *    gamma_o_molar = gamma_o_molal
 */
void  IdealMolalSoln::s_updateIMS_lnMolalityActCoeff() const
{
    double tmp;
    /*
     * Calculate the molalities. Currently, the molalities
     * may not be current with respect to the contents of the
     * State objects' data.
     */
    calcMolalities();

    double xmolSolvent = moleFraction(m_indexSolvent);
    double xx = std::max(m_xmolSolventMIN, xmolSolvent);

    if (IMS_typeCutoff_ == 0) {
        for (size_t k = 1; k < m_kk; k++) {
            IMS_lnActCoeffMolal_[k]= 0.0;
        }
        IMS_lnActCoeffMolal_[m_indexSolvent] = - log(xx) + (xx - 1.0)/xx;
        return;
    } else if (IMS_typeCutoff_ == 1) {
        if (xmolSolvent > 3.0 * IMS_X_o_cutoff_/2.0) {
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= 0.0;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = - log(xx) + (xx - 1.0)/xx;
            return;
        } else if (xmolSolvent < IMS_X_o_cutoff_/2.0) {
            tmp = log(xx * IMS_gamma_k_min_);
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = log(IMS_gamma_o_min_);
            return;
        } else {

            /*
             * If we are in the middle region, calculate the connecting polynomials
             */
            double xminus  = xmolSolvent - IMS_X_o_cutoff_/2.0;
            double xminus2 = xminus * xminus;
            double xminus3 = xminus2 * xminus;
            double x_o_cut2 = IMS_X_o_cutoff_ * IMS_X_o_cutoff_;
            double x_o_cut3 =  x_o_cut2 * IMS_X_o_cutoff_;

            double h2 = 3.5 * xminus2 /  IMS_X_o_cutoff_ - 2.0 * xminus3 / x_o_cut2;
            double h2_prime = 7.0 * xminus /  IMS_X_o_cutoff_ - 6.0 * xminus2 /  x_o_cut2;

            double h1 = (1.0 - 3.0 * xminus2 /  x_o_cut2 + 2.0 *  xminus3/ x_o_cut3);
            double h1_prime = (- 6.0 * xminus /  x_o_cut2 + 6.0 *  xminus2/ x_o_cut3);

            double h1_g = h1 / IMS_gamma_o_min_;
            double h1_g_prime  = h1_prime / IMS_gamma_o_min_;

            double alpha = 1.0 / (exp(1.0) * IMS_gamma_k_min_);
            double h1_f = h1 * alpha;
            double h1_f_prime  = h1_prime * alpha;

            double f = h2 + h1_f;
            double f_prime = h2_prime + h1_f_prime;

            double g = h2 + h1_g;
            double g_prime = h2_prime + h1_g_prime;

            tmp = (xmolSolvent/ g * g_prime + (1.0-xmolSolvent) / f * f_prime);
            double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
            double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

            tmp = log(xmolSolvent) + lngammak;
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = lngammao;
        }
    }

    // Exponentials - trial 2
    else if (IMS_typeCutoff_ == 2) {
        if (xmolSolvent > IMS_X_o_cutoff_) {
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= 0.0;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = - log(xx) + (xx - 1.0)/xx;
            return;
        } else {

            double xoverc = xmolSolvent/IMS_cCut_;
            double eterm = std::exp(-xoverc);

            double fptmp = IMS_bfCut_  - IMS_afCut_ / IMS_cCut_ - IMS_bfCut_*xoverc
                           + 2.0*IMS_dfCut_*xmolSolvent - IMS_dfCut_*xmolSolvent*xoverc;
            double f_prime = 1.0 + eterm*fptmp;
            double f = xmolSolvent + IMS_efCut_ + eterm * (IMS_afCut_ + xmolSolvent * (IMS_bfCut_ + IMS_dfCut_*xmolSolvent));

            double gptmp = IMS_bgCut_  - IMS_agCut_ / IMS_cCut_ - IMS_bgCut_*xoverc
                           + 2.0*IMS_dgCut_*xmolSolvent - IMS_dgCut_*xmolSolvent*xoverc;
            double g_prime = 1.0 + eterm*gptmp;
            double g = xmolSolvent + IMS_egCut_ + eterm * (IMS_agCut_ + xmolSolvent * (IMS_bgCut_ + IMS_dgCut_*xmolSolvent));

            tmp = (xmolSolvent / g * g_prime + (1.0 - xmolSolvent) / f * f_prime);
            double lngammak = -1.0 - log(f) + tmp * xmolSolvent;
            double lngammao =-log(g) - tmp * (1.0-xmolSolvent);

            tmp = log(xx) + lngammak;
            for (size_t k = 1; k < m_kk; k++) {
                IMS_lnActCoeffMolal_[k]= tmp;
            }
            IMS_lnActCoeffMolal_[m_indexSolvent] = lngammao;
        }
    }
    return;
}

/*
 * This internal function adjusts the lengths of arrays.
 *
 * This function is not virtual nor is it inherited
 */
void IdealMolalSoln::initLengths()
{
    m_kk = nSpecies();
    /*
     * Obtain the limits of the temperature from the species
     * thermo handler's limits.
     */
    m_pp.resize(m_kk);
    m_speciesMolarVolume.resize(m_kk);
    m_tmpV.resize(m_kk);
    IMS_lnActCoeffMolal_.resize(m_kk);
}


void  IdealMolalSoln::calcIMSCutoffParams_()
{
    IMS_afCut_ = 1.0 / (std::exp(1.0) *  IMS_gamma_k_min_);
    IMS_efCut_ = 0.0;
    bool converged = false;
    double oldV = 0.0;
    int its;
    for (its = 0; its < 100 && !converged; its++) {
        oldV = IMS_efCut_;
        IMS_afCut_ = 1.0 / (std::exp(1.0) * IMS_gamma_k_min_)  - IMS_efCut_;
        IMS_bfCut_ = IMS_afCut_ / IMS_cCut_ + IMS_slopefCut_ - 1.0;
        IMS_dfCut_ = ((- IMS_afCut_/IMS_cCut_ + IMS_bfCut_ - IMS_bfCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_afCut_ + IMS_X_o_cutoff_*(IMS_bfCut_ + IMS_dfCut_ * IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_efCut_ = - eterm * (tmp);
        if (fabs(IMS_efCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError(" IdealMolalSoln::calcCutoffParams_()",
                           " failed to converge on the f polynomial");
    }
    converged = false;
    double f_0 = IMS_afCut_ + IMS_efCut_;
    double f_prime_0 = 1.0 - IMS_afCut_ / IMS_cCut_ + IMS_bfCut_;
    IMS_egCut_ = 0.0;
    for (its = 0; its < 100 && !converged; its++) {
        oldV = IMS_egCut_;
        double lng_0 = -log(IMS_gamma_o_min_) -  f_prime_0 / f_0;
        IMS_agCut_ = exp(lng_0) - IMS_egCut_;
        IMS_bgCut_ = IMS_agCut_ / IMS_cCut_ + IMS_slopegCut_ - 1.0;
        IMS_dgCut_ = ((- IMS_agCut_/IMS_cCut_ + IMS_bgCut_ - IMS_bgCut_*IMS_X_o_cutoff_/IMS_cCut_)
                      /
                      (IMS_X_o_cutoff_*IMS_X_o_cutoff_/IMS_cCut_ - 2.0 * IMS_X_o_cutoff_));
        double tmp = IMS_agCut_ + IMS_X_o_cutoff_*(IMS_bgCut_ + IMS_dgCut_ *IMS_X_o_cutoff_);
        double eterm = std::exp(-IMS_X_o_cutoff_/IMS_cCut_);
        IMS_egCut_ = - eterm * (tmp);
        if (fabs(IMS_egCut_ - oldV) < 1.0E-14) {
            converged = true;
        }
    }
    if (!converged) {
        throw CanteraError(" IdealMolalSoln::calcCutoffParams_()",
                           " failed to converge on the g polynomial");
    }
}

}

