/**
 *  @file IonsFromNeutralVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   consist of ions whose thermodynamics is calculated from neutral molecule thermodynamics.
 *  (see \ref thermoprops
 * and class \link Cantera::IonsFromNeutralVPSSTP IonsFromNeutralVPSSTP\endlink).
 *
 * Header file for a derived class of ThermoPhase that handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties that are further based upon activities
 * based on the molality scale.  These include most of the methods for
 * calculating liquid electrolyte thermodynamics.
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_IONSFROMNEUTRALVPSSTP_H
#define CT_IONSFROMNEUTRALVPSSTP_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera
{

//! enums for molten salt ion solution types
/*!
 *  Types identify how complicated the solution is. If there
 *  is just mixing on one of the sublattices but not the other,
 *  then the math is considerably simpler.
 */
enum IonSolnType_enumType {
    cIonSolnType_PASSTHROUGH = 2000 ,
    cIonSolnType_SINGLEANION ,
    cIonSolnType_SINGLECATION ,
    cIonSolnType_MULTICATIONANION
};

/*!
 * The IonsFromNeutralVPSSTP is a derived class of ThermoPhase
 * that handles the specification of the chemical potentials for
 * ionic species, given a specification of the chemical potentials
 * for the same phase expressed in terms of combinations of the
 * ionic species that represent neutral molecules. It's expected
 * that the neutral molecules will be represented in terms of
 * an excess gibbs free energy approximation that is a derivative
 * of the GbbsExcessVPSSTP object. All of the e Excess Gibbs free
 *  energy formulations in this area employ
 *  symmetrical formulations.
 *
 *  This class is used for molten salts.
 *
 *  This object actually employs 4 different mole fraction types.
 *
 *  1. There is a mole fraction associated the the cations and
 *     anions and neutrals from this ThermoPhase object. This
 *     is the normal mole fraction vector for this object.
 *     Note, however, it isn't the appropriate mole fraction
 *     vector to use even for obtaining the correct ideal
 *     free energies of mixing.
 *  2. There is a mole fraction vector associated with the
 *     neutral molecule ThermoPhase object.
 *  3. There is a mole fraction vector associated with the
 *     cation lattice.
 *  4. There is a mole fraction vector associated with the
 *     anion lattice
 *
 *  This object can translate between any of the four mole
 *  fraction representations.
 */
class IonsFromNeutralVPSSTP : public GibbsExcessVPSSTP
{
public:

    //! @name Constructors
    //! @{

    /*!
     * Default constructor
     */
    IonsFromNeutralVPSSTP();

    //! Construct and initialize an IonsFromNeutralVPSSTP object
    //! directly from an ASCII input file
    /*!
     *  This constructor is a shell around the routine initThermo(), with a
     *  reference to the XML database to get the info for the phase.
     *
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     * @param neutralPhase   The object takes a neutralPhase ThermoPhase
     *                       object as input. It can either take a pointer
     *                       to an existing object in the parameter list,
     *                       in which case it does not own the object, or
     *                       it can construct a neutral Phase as a slave
     *                       object, in which case, it does own the slave
     *                       object, for purposes of who gets to destroy
     *                       the object.
     *                       If this parameter is zero, then a slave
     *                       neutral phase object is created and used.
     */
    IonsFromNeutralVPSSTP(const std::string& inputFile,
                          const std::string& id = "",
                          ThermoPhase* neutralPhase = 0);

    //! Construct and initialize an IonsFromNeutralVPSSTP object
    //! directly from an XML database
    /*!
     *  @param phaseRoot XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     * @param neutralPhase   The object takes a neutralPhase ThermoPhase
     *                       object as input. It can either take a pointer
     *                       to an existing object in the parameter list,
     *                       in which case it does not own the object, or
     *                       it can construct a neutral Phase as a slave
     *                       object, in which case, it does own the slave
     *                       object, for purposes of who gets to destroy
     *                       the object.
     *                       If this parameter is zero, then a slave
     *                       neutral phase object is created and used.
     */
    IonsFromNeutralVPSSTP(XML_Node& phaseRoot, const std::string& id = "",
                          ThermoPhase* neutralPhase = 0);


    //! Copy constructor
    /*!
     * @param b class to be copied
     */
    IonsFromNeutralVPSSTP(const IonsFromNeutralVPSSTP& b);

    /// Assignment operator
    /*!
     * @param b class to be copied.
     */
    IonsFromNeutralVPSSTP& operator=(const IonsFromNeutralVPSSTP& b);

    /// Destructor.
    virtual ~IonsFromNeutralVPSSTP();

    //! Duplication routine for objects which inherit from ThermoPhase.
    /*!
     *  This virtual routine can be used to duplicate ThermoPhase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    // @}

    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an
    /// input file.

    //! Initialization of an IonsFromNeutralVPSSTP phase using an xml file
    /*!
     * This routine is a precursor to initThermo(XML_Node*)
     * routine, which does most of the work.
     *
     * @param inputFile XML file containing the description of the phase
     *
     * @param id  Optional parameter identifying the name of the
     *            phase. If none is given, the first XML
     *            phase element will be used.
     */
    void constructPhaseFile(std::string inputFile, std::string id);

    //!   Import and initialize an IonsFromNeutralVPSSTP phase
    //!   specification in an XML tree into the current object.
    /*!
     *   Here we read an XML description of the phase.
     *   We import descriptions of the elements that make up the
     *   species in a phase.
     *   We import information about the species, including their
     *   reference state thermodynamic polynomials. We then freeze
     *   the state of the species.
     *
     *   Then, we read the species molar volumes from the xml
     *   tree to finish the initialization.
     *
     * @param phaseNode This object must be the phase node of a complete XML tree
     *             description of the phase, including all of the
     *             species data. In other words while "phase" must
     *             point to an XML phase object, it must have
     *             sibling nodes "speciesData" that describe
     *             the species in the phase.
     *
     * @param id   ID of the phase. If nonnull, a check is done
     *             to see if phaseNode is pointing to the phase
     *             with the correct id.
     */
    void constructPhaseXML(XML_Node& phaseNode, std::string id);

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

    //! @}
    //! @name Molar Thermodynamic Properties
    //! @{

    //! Return the Molar enthalpy. Units: J/kmol.
    /*!
     * This is calculated from the partial molar enthalpies of the species.
     */
    virtual doublereal enthalpy_mole() const;

    /**
     * Molar internal energy. J/kmol.
     *
     * This is calculated from the soln enthalpy and then
     * subtracting pV.
     */
    virtual doublereal intEnergy_mole() const;

    //! Molar entropy. Units: J/kmol/K.
    virtual doublereal entropy_mole() const;

    //! Molar Gibbs free Energy for an ideal gas. Units =  J/kmol.
    virtual doublereal gibbs_mole() const;

    //! Molar heat capacity at constant pressure. Units: J/kmol/K.
    virtual doublereal cp_mole() const;

    //! Molar heat capacity at constant volume. Units: J/kmol/K.
    virtual doublereal cv_mole() const;

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

    //! Get the array of non-dimensional molar-based activity coefficients at
    //! the current solution temperature, pressure, and solution concentration.
    /*!
     * @param ac Output vector of activity coefficients. Length: m_kk.
     */
    virtual void getActivityCoefficients(doublereal* ac) const;

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
     * \bar h_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
     * \f]
     *
     *  @param hbar  Output vector of species partial molar enthalpies.
     *               Length: m_kk. Units: J/kmol
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
     *
     *  @param sbar  Output vector of species partial molar entropies.
     *               Length: m_kk. Units: J/kmol/K
     */
    virtual void getPartialMolarEntropies(doublereal* sbar) const;


    //! Get the change in activity coefficients w.r.t. change in state (temp, mole fraction, etc.) along
    //! a line in parameter space or along a line in physical space
    /*!
     *
     * @param dTds           Input of temperature change along the path
     * @param dXds           Input vector of changes in mole fraction along the path. length = m_kk
     *                       Along the path length it must be the case that the mole fractions sum to one.
     * @param dlnActCoeffds  Output vector of the directional derivatives of the
     *                       log Activity Coefficients along the path. length = m_kk
     */
    virtual void getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds,
                                  doublereal* dlnActCoeffds) const;

    //! Get the array of log concentration-like derivatives of the
    //! log activity coefficients - diagonal component
    /*!
     * For ideal mixtures (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to the logarithm of the mole
     * fraction. This quantity is to be used in conjunction with derivatives
     * of that concentration-like variable when the derivative of the chemical
     * potential is taken.
     *
     * units = dimensionless
     *
     * @param dlnActCoeffdlnX_diag    Output vector of log(mole fraction)
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const;

    //! Get the array of log concentration-like derivatives of the
    //! log activity coefficients - diagonal components
    /*!
     * For ideal mixtures (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the logarithm of the
     * activity coefficient with respect to the logarithm of the species mole
     * numbe. This routine just does the diagonal entries.
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of diagonal components of the log(mole fraction)
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnN_diag(doublereal* dlnActCoeffdlnN_diag) const;

    //! Get the array of derivatives of the ln activity coefficients with respect to the ln species mole numbers
    /*!
     * Implementations should take the derivative of the logarithm of the activity coefficient with respect to a
     * log of a species mole number (with all other species mole numbers held constant)
     *
     *  units = 1 / kmol
     *
     *  dlnActCoeffdlnN[ ld * k  + m]  will contain the derivative of log act_coeff for the <I>m</I><SUP>th</SUP>
     *                                 species with respect to the number of moles of the <I>k</I><SUP>th</SUP> species.
     *
     * \f[
     *        \frac{d \ln(\gamma_m) }{d \ln( n_k ) }\Bigg|_{n_i}
     * \f]
     *
     * @param ld               Number of rows in the matrix
     * @param dlnActCoeffdlnN    Output vector of derivatives of the
     *                         log Activity Coefficients. length = m_kk * m_kk
     */
    virtual void getdlnActCoeffdlnN(const size_t ld, doublereal* const dlnActCoeffdlnN) ;
    //! @}

    //! Get the Salt Dissociation Coefficients
    //! Returns the vector of dissociation coefficients and vector of charges
    /*!
     *  @param fm_neutralMolec_ions Returns the formula matrix for the composition of neutral molecules
     *                              in terms of the ions.
     *  @param charges              Returns a vector containing the charges of all species in this phase
     *  @param neutMolIndex         Returns the vector fm_invert_ionForNeutral
     *                               This is the mapping between ion species and neutral molecule for quick invert.
     */
    void getDissociationCoeffs(vector_fp& fm_neutralMolec_ions, vector_fp& charges, std::vector<size_t>& neutMolIndex) const;

    //! Return the current value of the neutral mole fraction vector
    /*!
     *   @param neutralMoleculeMoleFractions  Vector of neutral molecule mole fractions.
     */
    void getNeutralMolecMoleFractions(vector_fp& neutralMoleculeMoleFractions) const {
        neutralMoleculeMoleFractions = NeutralMolecMoleFractions_;
    }

    //! Calculate neutral molecule mole fractions
    /*!
     *  This routine calculates the neutral molecule mole
     *  fraction given the vector of ion mole fractions,
     *  i.e., the mole fractions from this ThermoPhase.
     *  Note, this routine basically assumes that there
     *  is charge neutrality. If there isn't, then it wouldn't
     *  make much sense.
     *
     *  for the case of  cIonSolnType_SINGLEANION, some slough
     *  in the charge neutrality is allowed. The cation number
     *  is followed, while the difference in charge neutrality
     *  is dumped into the anion mole number to fix the imbalance.
     *
     *  @param  dx  input vector of ion mole fraction gradients
     *  @param  dy  output Vector of neutral molecule mole fraction gradients
     */
    void getNeutralMoleculeMoleGrads(const doublereal* const dx, doublereal* const dy) const;

    //! Get the list of cations in this object
    /*!
     *  @param cation  List of cations
     */
    void getCationList(std::vector<size_t>& cation) const {
        cation=cationList_;
    }

    //! Get the list of anions in this object
    /*!
     *  @param anion  List of anions
     */
    void getAnionList(std::vector<size_t>& anion) const {
        anion=anionList_;
    }

    /**
     * @name Setting the State
     * These methods set all or part of the thermodynamic state.
     * @{
     */

    virtual void setTemperature(const doublereal t);
    virtual void setPressure(doublereal p);

    //! Set the temperature (K) and pressure (Pa)
    /*!
     * Setting the pressure may involve the solution of a nonlinear equation.
     *
     * @param t    Temperature (K)
     * @param p    Pressure (Pa)
     */
    virtual void setState_TP(doublereal t, doublereal p);


    //! Calculate ion mole fractions from neutral molecule mole fractions.
    /*!
     *  @param mf Dump the mole fractions into this vector.
     */
    virtual void calcIonMoleFractions(doublereal* const mf) const;

    //! Calculate neutral molecule mole fractions
    /*!
     *  This routine calculates the neutral molecule mole
     *  fraction given the vector of ion mole fractions,
     *  i.e., the mole fractions from this ThermoPhase.
     *  Note, this routine basically assumes that there
     *  is charge neutrality. If there isn't, then it wouldn't
     *  make much sense.
     *
     *  for the case of  cIonSolnType_SINGLEANION, some slough
     *  in the charge neutrality is allowed. The cation number
     *  is followed, while the difference in charge neutrality
     *  is dumped into the anion mole number to fix the imbalance.
     */
    virtual void calcNeutralMoleculeMoleFractions() const;

    /**
      * Set the mass fractions to the specified values, and then
      * normalize them so that they sum to 1.0.
      * @param y Array of unnormalized mass fraction values (input).
      * Must have a length greater than or equal to the number of
      * species.
      *
      * @param y  Input vector of mass fractions.
      *           Length is m_kk.
      */
    virtual void setMassFractions(const doublereal* const y);

    /**
     * Set the mass fractions to the specified values without
     * normalizing. This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of
     * equations.
     *
     * @param y  Input vector of mass fractions.
     *           Length is m_kk.
     */
    virtual void setMassFractions_NoNorm(const doublereal* const y);

    /**
     * Set the mole fractions to the specified values, and then
     * normalize them so that they sum to 1.0.
     * @param x Array of unnormalized mole fraction values (input).
     * Must have a length greater than or equal to the number of
     * species.
     *
     * @param x  Input vector of mole fractions.
     *           Length is m_kk.
     */
    virtual void setMoleFractions(const doublereal* const x);

    /**
     * Set the mole fractions to the specified values without
     * normalizing. This is useful when the normalization
     * condition is being handled by some other means, for example
     * by a constraint equation as part of a larger set of
     * equations.
     *
     * @param x  Input vector of mole fractions.
     *           Length is m_kk.
     */
    virtual void setMoleFractions_NoNorm(const doublereal* const x);

    /**
     * Set the concentrations to the specified values within the phase.
     *
     * @param c The input vector to this routine is in dimensional
     *        units. For volumetric phases c[k] is the
     *        concentration of the kth species in kmol/m3.
     *        For surface phases, c[k] is the concentration
     *        in kmol/m2. The length of the vector is the number
     *        of species in the phase.
     */
    virtual void setConcentrations(const doublereal* const c);

    //@}

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


private:
    //! Initialize lengths of local variables after all species have
    //! been identified.
    void initLengths();

    //! Update the activity coefficients
    /*!
     * This function will be called to update the internally stored
     * natural logarithm of the activity coefficients
     */
    void s_update_lnActCoeff() const;

    //! Update the temperature derivative of the ln activity coefficients
    /*!
     * This function will be called to update the internally stored
     * temperature derivative of the natural logarithm of the activity coefficients
     */
    void s_update_dlnActCoeffdT() const;

    //! Update the change in the ln activity coefficients
    /*!
     * This function will be called to update the internally stored
     * change of the natural logarithm of the activity coefficients
     * w.r.t a change in state (temp, mole fraction, etc)
     */
    void s_update_dlnActCoeff() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(mole fraction)
    /*!
     * This function will be called to update the internally stored
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the mole fractions.
     */
    void s_update_dlnActCoeff_dlnX_diag() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(number of moles) - diagonal components
    /*!
     * This function will be called to update the internally stored
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the number of moles of given species.
     */
    void s_update_dlnActCoeff_dlnN_diag() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(number of moles) - diagonal components
    /*!
     * This function will be called to update the internally stored
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the number of moles of given species.
     */
    void s_update_dlnActCoeff_dlnN() const;


private:
    //! Error function
    /*!
     *  Print an error string and exit
     *
     * @param msg  Message to be printed
     */
    doublereal err(const std::string& msg) const;

protected:

    //! Ion solution type
    /*!
     *  There is either mixing on the anion, cation, or both lattices.
     *  There is also a passthrough option
     *
     *  Defaults to cIonSolnType_SINGLEANION, so that LiKCl can be hardwired
     */
    IonSolnType_enumType ionSolnType_;

    //! Number of neutral molecule species
    /*!
     *  This is equal to the number of species in the
     *  neutralMoleculePhase_ ThermoPhase.
     */
    size_t numNeutralMoleculeSpecies_;

    //! Index of special species
    size_t indexSpecialSpecies_;

    //! Index of special species
    size_t indexSecondSpecialSpecies_;

    //! Formula Matrix for composition of neutral molecules
    //! in terms of the molecules in this ThermoPhase
    /*!
     *       fm_neutralMolec_ions[ i + jNeut * m_kk ]
     *
     *  This is the number of ions of type i in the neutral
     *  molecule jNeut.
     */
    std::vector<double> fm_neutralMolec_ions_;

    //! Mapping between ion species and neutral molecule for quick invert.
    /*!
     * fm_invert_ionForNeutral returns vector of int. Each element represents
     * an ionic species and stores the value of the corresponding neutral
     * molecule
     *
     *  For the case of fm_invert_simple_ = true, we assume that there
     *  is a quick way to invert the formula matrix so that we can
     *  quickly calculate the neutral molecule mole fraction
     *  given the ion mole fraction vector.
     *
     *  We assume that for a selected set of ion species, that that
     *  ion is only in the neutral molecule, jNeut.
     *
     *  therefore,
     *
     *      NeutralMolecMoleFractions_[jNeut] += moleFractions_[i_ion] / fmij;
     *
     *  where fmij is the number of ions in neutral molecule jNeut.
     *
     *  Thus, we formulate the neutral molecule mole fraction NeutralMolecMoleFractions_[]
     *  vector from this association. We further assume that there are
     *  no other associations.  If  fm_invert_simple_ is not true,
     *  then we need to do a formal inversion which takes a great
     *  deal of time and is not currently implemented.
     */
    std::vector<size_t> fm_invert_ionForNeutral;

    //! Mole fractions using the Neutral Molecule Mole fraction basis
    mutable std::vector<doublereal> NeutralMolecMoleFractions_;

    //! List of the species in this ThermoPhase which are cation species
    std::vector<size_t> cationList_;

    //! Number of cation species
    int numCationSpecies_;

    //! List of the species in this ThermoPhase which are anion species
    std::vector<size_t> anionList_;

    //! Number of anion species
    int numAnionSpecies_;

    //! List of the species in this ThermoPhase which are passed
    //! through to the neutralMoleculePhase ThermoPhase.
    /*!
     *  These have neutral charges.
     */
    std::vector<size_t> passThroughList_;

    //! Number of the species in this ThermoPhase which are passed
    //! through to the neutralMoleculePhase ThermoPhase
    size_t numPassThroughSpecies_;

public:
    //! This is a pointer to the neutral Molecule Phase
    /*!
     *  If the variable, IOwnNThermoPhase_ is true, then we own
     *  the pointer. If not, then this is considered a shallow pointer.
     */
    ThermoPhase* neutralMoleculePhase_;

private:
    GibbsExcessVPSSTP* geThermo;
    // Temporary vectors that I don't want to allocate every time the function is called
    mutable vector_fp y_;
    mutable vector_fp dlnActCoeff_NeutralMolecule_;
    mutable vector_fp dX_NeutralMolecule_;

    //! If true then we own the underlying neutral Molecule Phase
    /*!
     *  If this is false, then the neutral molecule phase is considered
     *  as a shallow pointer.
     */
    bool IOwnNThermoPhase_;

    //! Temporary mole fraction vector
    mutable std::vector<doublereal> moleFractionsTmp_;

    //! Storage vector for the neutral molecule chemical potentials
    /*!
     *  This vector is used as a temporary storage area when calculating the ion chemical
     *  potentials.
     *
     *  - Units = Joules/kmol
     *  - Length =  numNeutralMoleculeSpecies_
     */
    mutable std::vector<doublereal> muNeutralMolecule_;

    //! Storage vector for the neutral molecule ln activity coefficients
    /*!
     *  This vector is used as a temporary storage area when calculating the ion chemical
     *  potentials and activity coefficients
     *
     *  - Units = none
     *  - Length =  numNeutralMoleculeSpecies_
     */
    mutable std::vector<doublereal> lnActCoeff_NeutralMolecule_;

    //! Storage vector for the neutral molecule d ln activity coefficients dT
    /*!
     *  This vector is used as a temporary storage area when calculating the ion derivatives
     *
     *  - Units =  1/Kelvin
     *  - Length =  numNeutralMoleculeSpecies_
     */
    mutable std::vector<doublereal> dlnActCoeffdT_NeutralMolecule_;

    //! Storage vector for the neutral molecule d ln activity coefficients dX - diagonal component
    /*!
     *  This vector is used as a temporary storage area when calculating the ion derivatives
     *
     *  - Units =  none
     *  - Length =  numNeutralMoleculeSpecies_
     */
    mutable std::vector<doublereal> dlnActCoeffdlnX_diag_NeutralMolecule_;

    //! Storage vector for the neutral molecule d ln activity coefficients dlnN - diagonal component
    /*!
     *  This vector is used as a temporary storage area when calculating the ion derivatives
     *
     *  - Units =  none
     *  - Length =  numNeutralMoleculeSpecies_
     */
    mutable std::vector<doublereal> dlnActCoeffdlnN_diag_NeutralMolecule_;

    //! Storage vector for the neutral molecule d ln activity coefficients dlnN
    /*!
     *  This vector is used as a temporary storage area when calculating the ion derivatives
     *
     *  - Units =  none
     *  - Length =  numNeutralMoleculeSpecies_
     */
    mutable Array2D dlnActCoeffdlnN_NeutralMolecule_;
};

}

#endif
