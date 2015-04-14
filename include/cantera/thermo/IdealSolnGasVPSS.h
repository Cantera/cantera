/**
 *  @file IdealSolnGasVPSS.h
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::IdealSolnGasVPSS IdealSolnGasVPSS\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_IDEALSOLNGASVPSS_H
#define CT_IDEALSOLNGASVPSS_H

#include "VPStandardStateTP.h"

namespace Cantera
{
/*!
  * @name CONSTANTS
  * Models for the Standard State of an IdealSolnPhase
  */
//@{
const int cIdealSolnGasPhaseG = 6009;
const int cIdealSolnGasPhase0 = 6010;
const int cIdealSolnGasPhase1 = 6011;
const int cIdealSolnGasPhase2 = 6012;

/**
 * @ingroup thermoprops
 *
 * An ideal solution or an ideal gas approximation of a phase. Uses variable
 * pressure standard state methods for calculating thermodynamic properties.
 */
class IdealSolnGasVPSS : public VPStandardStateTP
{
public:
    /*!
     * @name Constructors and Duplicators for IdealSolnGasVPSS
     */
    //! @{

    /// Constructor.
    IdealSolnGasVPSS();

    /// Create an object from an XML input file
    IdealSolnGasVPSS(const std::string& infile, std::string id="");

    /// Copy Constructor.
    IdealSolnGasVPSS(const IdealSolnGasVPSS&);

    /// Assignment operator
    IdealSolnGasVPSS& operator=(const IdealSolnGasVPSS&);

    //! Duplication routine
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    //@}
    //! @name  Utilities (IdealSolnGasVPSS)
    //@{
    /**
     * Equation of state type flag. The base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Constants defined for this purpose are
     * listed in mix_defs.h.
     */
    virtual int eosType() const;

    //! @}
    //! @name Molar Thermodynamic Properties
    //! @{

    /// Molar enthalpy. Units: J/kmol.
    doublereal enthalpy_mole() const;

    /// Molar entropy. Units: J/kmol/K.
    doublereal entropy_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
    doublereal cv_mole() const;

    //! @}
    //! @name Mechanical Properties
    //! @{

    //! Set the pressure in the fluid
    /*!
     * @param p pressure in pascals.
     */
    void setPressure(doublereal p);

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

protected:
    /**
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
     * solution the partial molar volumes are equal to the
     * species standard state molar volumes.
     * The species molar volumes may be functions
     * of temperature and pressure.
     */
    virtual void calcDensity();
    //! @}

public:
    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in the
     * thermodynamic functions.  These activity (or generalized)
     * concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions. Note that they may
     * or may not have units of concentration --- they might be
     * partial pressures, mole fractions, or surface coverages,
     * for example.
     *
     * @param c Output array of generalized concentrations. The
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Returns the standard concentration \f$ C^0_k \f$, which is used to normalize
    //! the generalized concentration.
    /*!
     * This is defined as the concentration  by which the generalized
     * concentration is normalized to produce the activity.
     * In many cases, this quantity will be the same for all species in a phase.
     * Since the activity for an ideal gas mixture is
     * simply the mole fraction, for an ideal gas \f$ C^0_k = P/\hat R T \f$.
     *
     * @param k Optional parameter indicating the species. The default
     *          is to assume this refers to species 0.
     * @return
     *   Returns the standard Concentration in units of m3 kmol-1.
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Returns the units of the standard and generalized concentrations.
    /*!
     * Note they have the same units, as their
     * ratio is defined to be equal to the activity of the kth
     * species in the solution, which is unitless.
     *
     * This routine is used in print out applications where the
     * units are needed. Usually, MKS units are assumed throughout
     * the program and in the XML input files.
     *
     * The base ThermoPhase class assigns the default quantities
     * of (kmol/m3) for all species.
     * Inherited classes are responsible for overriding the default
     * values if necessary.
     *
     * @param uA Output vector containing the units
     *  uA[0] = kmol units - default  = 1
     *  uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                dimensions in the Phase class.
     *  uA[2] = kg   units - default  = 0;
     *  uA[3] = Pa(pressure) units - default = 0;
     *  uA[4] = Temperature units - default = 0;
     *  uA[5] = time units - default = 0
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     * @deprecated To be removed after Cantera 2.2.
     */
    virtual void getUnitsStandardConc(double* uA, int k = 0,
                                      int sizeUA = 6) const;

    //! Get the array of non-dimensional activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     *  For ideal gases, the activity coefficients are all equal to one.
     *
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;


    /// @name  Partial Molar Properties of the Solution
    //@{

    //! Get the array of non-dimensional species chemical potentials
    //! These are partial molar Gibbs free energies.
    /*!
     * \f$ \mu_k / \hat R T \f$.
     * Units: unitless
     *
     * We close the loop on this function, here, calling
     * getChemPotentials() and then dividing by RT. No need for child
     * classes to handle.
     *
     * @param mu    Output vector of  non-dimensional species chemical potentials
     *              Length: m_kk.
     */
    void getChemPotentials_RT(doublereal* mu) const;

    //! Get the species chemical potentials. Units: J/kmol.
    /*!
     * This function returns a vector of chemical potentials of the
     * species in solution at the current temperature, pressure
     * and mole fraction of the solution.
     *
     * @param mu  Output vector of species chemical
     *            potentials. Length: m_kk. Units: J/kmol
     */
    virtual void getChemPotentials(doublereal* mu) const;

    //!  Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param hbar    Output vector of species partial molar enthalpies.
     *                Length: m_kk. units are J/kmol.
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Get the species partial molar entropies. Units: J/kmol/K.
    /*!
     * @param sbar    Output vector of species partial molar entropies.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Get the species partial molar enthalpies. Units: J/kmol.
    /*!
     * @param ubar    Output vector of species partial molar internal energies.
     *                Length = m_kk. units are J/kmol.
     */
    virtual void getPartialMolarIntEnergies(doublereal* ubar) const;

    //! Get the partial molar heat capacities Units: J/kmol/K
    /*!
     * @param cpbar   Output vector of species partial molar heat capacities
     *                at constant pressure.
     *                Length = m_kk. units are J/kmol/K.
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Get the species partial molar volumes. Units: m^3/kmol.
    /*!
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;
    //@}

public:
    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() when processing a phase
     * definition in an input file. It should be overloaded in subclasses to
     * set any parameters that are specific to that particular phase model.
     *
     * @param thermoNode An XML_Node object corresponding to
     *                   the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& thermoNode);

    //! @internal Initialize the object
    /*!
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase().
     */
    virtual void initThermo();

    //!This method is used by the ChemEquil equilibrium solver.
    /*!
     * It sets the state such that the chemical potentials satisfy
     * \f[ \frac{\mu_k}{\hat R T} = \sum_m A_{k,m}
     * \left(\frac{\lambda_m} {\hat R T}\right) \f] where
     * \f$ \lambda_m \f$ is the element potential of element m. The
     * temperature is unchanged.  Any phase (ideal or not) that
     * implements this method can be equilibrated by ChemEquil.
     *
     * @param lambda_RT Input vector of dimensionless element potentials
     *                  The length is equal to nElements().
     */
    void setToEquilState(const doublereal* lambda_RT);

    //!   Initialize a ThermoPhase object, potentially reading activity
    //!   coefficient information from an XML database.
    /*!
     * This routine initializes the lengths in the current object and
     * then calls the parent routine.
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called just prior to returning
     * from function importPhase().
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
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

private:
    //!  @internal Initialize the internal lengths in this object.
    /*!
     * Note this is not a virtual function and only handles this object
     */
    void initLengths();

    //@}

protected:
    //! boolean indicating what ideal solution this is
    /*!
     *  - 1 = ideal gas
     *  - 0 = ideal soln
     */
    int m_idealGas;

    //! form of the generalized concentrations
    /*!
     *    - 0 unity
     *    - 1 1/V_k
     *    - 2 1/V_0
     */
    int m_formGC;

    //! Temporary storage - length = m_kk.
    vector_fp m_pp;
};
}

#endif
