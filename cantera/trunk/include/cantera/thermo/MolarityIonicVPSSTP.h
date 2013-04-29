/**
 *  @file MolarityIonicVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ gibbs excess free energy based formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MolarityIonicVPSSTP MolarityIonicVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon activities
 * based on the molarity scale.  In this class, we expect that there are
 * ions, but they are treated on the molarity scale.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_MOLARITYIONICVPSSTP_H
#define CT_MOLARITYIONICVPSSTP_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 */

/*!
 *  MolarityIonicVPSSTP is a derived class of GibbsExcessVPSSTP that handles
 *  variable pressure standard state methods for calculating
 *  thermodynamic properties that are further based on
 *  expressing the Excess Gibbs free energy as a function of
 *  the mole fractions (or pseudo mole fractions) of the constituents.
 *  This category is the workhorse for describing ionic systems which are not on the molality scale.
 *
 *  This class adds additional functions onto the ThermoPhase interface
 *  that handles the calculation of the excess Gibbs free energy. The ThermoPhase
 *  class includes a member function, ThermoPhase::activityConvention()
 *  that indicates which convention the activities are based on. The
 *  default is to assume activities are based on the molar convention.
 *  That default is used here.
 *
 *  All of the Excess Gibbs free energy formulations in this area employ
 *  symmetrical formulations.
 *
 *  This layer will massage the mole fraction vector to implement
 *  cation and anion based mole numbers in an optional manner, such that
 *  it is expected that there exists a charge balance at all times.
 *  One of the ions must be a "special ion" in the sense that its' thermodynamic
 *  functions are set to zero, and the thermo functions of all other
 *  ions are based on a valuation relative to that special ion.
 *
 */
class MolarityIonicVPSSTP : public GibbsExcessVPSSTP
{

public:
    /// Constructor
    /*!
     * This doesn't do much more than initialize constants with
     * default values for water at 25C. Water molecular weight
     * comes from the default elements.xml file. It actually
     * differs slightly from the IAPWS95 value of 18.015268. However,
     * density conservation and therefore element conservation
     * is the more important principle to follow.
     */
    MolarityIonicVPSSTP();

    //! Construct and initialize a MolarityIonicVPSSTP ThermoPhase object
    //! directly from an XML input file
    /*!
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    MolarityIonicVPSSTP(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize a MolarityIonicVPSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    MolarityIonicVPSSTP(XML_Node& phaseRef, const std::string& id = "");

    //! Copy constructor
    /*!
     * @param b class to be copied
     */
    MolarityIonicVPSSTP(const  MolarityIonicVPSSTP& b);

    /// Assignment operator
    /*!
     * @param b class to be copied.
     */
    MolarityIonicVPSSTP& operator=(const MolarityIonicVPSSTP& b);

    //! Duplication routine for objects which inherit from  ThermoPhase.
    /*!
     *  This virtual routine can be used to duplicate thermophase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    //! @name  Utilities
    //! @{

    //! Equation of state type flag.
    /*!
     * The ThermoPhase base class returns
     * zero. Subclasses should define this to return a unique
     * non-zero value. Known constants defined for this purpose are
     * listed in mix_defs.h. The MolalityVPSSTP class also returns
     * zero, as it is a non-complete class.
     */
    virtual int eosType() const;

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     * The activity \f$a_k\f$ of a species in solution is
     * related to the chemical potential by \f[ \mu_k = \mu_k^0(T)
     * + \hat R T \log a_k. \f] The quantity \f$\mu_k^0(T,P)\f$ is
     * the chemical potential at unit activity, which depends only
     * on temperature and pressure.
     * @{
     */

    //! Get the array of non-dimensional molar-based ln activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param lnac Output vector of ln activity coefficients. Length: m_kk.
     */
    virtual void getLnActivityCoefficients(doublereal* lnac) const;

    //@}
    /// @name  Partial Molar Properties of the Solution
    //@{

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

    /**
     * Get the species electrochemical potentials.
     * These are partial molar quantities.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the
     * to each chemical potential.
     *
     * Units: J/kmol
     *
     * @param mu     output vector containing the species electrochemical potentials.
     *               Length: m_kk.
     */
    void getElectrochemPotentials(doublereal* mu) const;

    //! Returns an array of partial molar enthalpies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the
     * standard state enthalpies modified by the derivative of the
     * molality-based activity coefficient wrt temperature
     *
     *  \f[
     *   \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *  \f]
     *
     * @param hbar  Vector of returned partial molar enthalpies
     *              (length m_kk, units = J/kmol)
     */
    virtual void getPartialMolarEnthalpies(doublereal* hbar) const;

    //! Returns an array of partial molar entropies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the
     * standard state enthalpies modified by the derivative of the
     * activity coefficient wrt temperature
     *
     *  \f[
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *  \f]
     *
     * @param sbar  Vector of returned partial molar entropies
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;

    //! Returns an array of partial molar entropies for the species
    //! in the mixture.
    /*!
     * Units (J/kmol)
     *
     * For this phase, the partial molar enthalpies are equal to the
     * standard state enthalpies modified by the derivative of the
     * activity coefficient wrt temperature
     *
     *  \f[
     *   ???????????????
     *   \bar s_k(T,P) = s^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     *                              - R \ln( \gamma_k X_k)
     *                              - R T \frac{d \ln(\gamma_k) }{dT}
     *   ???????????????
     *  \f]
     *
     * @param cpbar  Vector of returned partial molar heat capacities
     *              (length m_kk, units = J/kmol/K)
     */
    virtual void getPartialMolarCp(doublereal* cpbar) const;

    //! Return an array of partial molar volumes for the
    //! species in the mixture. Units: m^3/kmol.
    /*!
     *  Frequently, for this class of thermodynamics representations,
     *  the excess Volume due to mixing is zero. Here, we set it as
     *  a default. It may be overridden in derived classes.
     *
     *  @param vbar   Output vector of species partial molar volumes.
     *                Length = m_kk. units are m^3/kmol.
     */
    virtual void getPartialMolarVolumes(doublereal* vbar) const;

    //@}

    //! Calculate pseudo binary mole fractions
    virtual void calcPseudoBinaryMoleFractions() const;

    /// @name Initialization
    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see files importCTML.cpp and
    /// ThermoFactory.cpp.
    /// @{

    /*!
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
     * @see importCTML.cpp
     */
    virtual void initThermo();

    /**
     *   Import and initialize a ThermoPhase object
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
    void initThermoXML(XML_Node& phaseNode, const std::string& id);
    //! @}

    //! returns a summary of the state of the phase as a string
    /*!
     * @param show_thermo If true, extra information is printed out
     *                    about the thermodynamic state of the system.
     */
    virtual std::string report(bool show_thermo = true) const;

private:
    //! Initialize lengths of local variables after all species have been identified.
    void initLengths();

    //! Process an XML node called "binaryNeutralSpeciesParameters"
    /*!
     * This node contains all of the parameters necessary to describe
     * the Redlich-Kister model for a particular binary interaction.
     * This function reads the XML file and writes the coefficients
     * it finds to an internal data structures.
     *
     * @param xmlBinarySpecies  Reference to the XML_Node named "binaryNeutralSpeciesParameters"
     *                          containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node& xmlBinarySpecies);

    //! Update the activity coefficients
    /*!
     * This function will be called to update the internally stored
     * natural logarithm of the activity coefficients
     */
    void s_update_lnActCoeff() const;

    //! Update the derivative of the log of the activity coefficients wrt T
    /*!
     * This function will be called to update the internally stored
     * derivative of the natural logarithm of the activity coefficients
     * wrt temperature.
     */
    void s_update_dlnActCoeff_dT() const;

    //! Internal routine that calculates the derivative of the activity coefficients wrt
    //! the mole fractions.
    /*!
     *  This routine calculates the the derivative of the activity coefficients wrt to mole fraction
     *  with all other mole fractions held constant. This is strictly not permitted. However, if the
     *  resulting matrix is multiplied by a permissible deltaX vector then everything is ok.
     *
     *  This is the natural way to handle concentration derivatives in this routine.
     */
    void s_update_dlnActCoeff_dX_() const;

private:
    //! Error function
    /*!
     *  Print an error string and exit
     *
     * @param msg  Message to be printed
     */
    doublereal err(const std::string& msg) const;

protected:
    // Pseudobinary type
    /*!
     *  - `PBTYPE_PASSTHROUGH` - All species are passthrough species
     *  - `PBTYPE_SINGLEANION` - there is only one anion in the mixture
     *  - `PBTYPE_SINGLECATION` - there is only one cation in the mixture
     *  - `PBTYPE_MULTICATIONANION` - Complex mixture
     */
    int PBType_;

    //! Number of pseudo binary species
    size_t numPBSpecies_;

    //! index of special species
    size_t indexSpecialSpecies_;

    mutable std::vector<doublereal> PBMoleFractions_;

    //! Vector of cation indices in the mixture
    std::vector<size_t> cationList_;

    //! Number of cations in the mixture
    size_t numCationSpecies_;

    std::vector<size_t> anionList_;
    size_t numAnionSpecies_;

    std::vector<size_t> passThroughList_;
    size_t numPassThroughSpecies_;
    size_t neutralPBindexStart;

    mutable std::vector<doublereal> moleFractionsTmp_;
};

#define  PBTYPE_PASSTHROUGH        0
#define  PBTYPE_SINGLEANION        1
#define  PBTYPE_SINGLECATION       2
#define  PBTYPE_MULTICATIONANION   3

}

#endif
