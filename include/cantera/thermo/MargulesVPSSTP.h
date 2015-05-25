/**
 *  @file  MargulesVPSSTP.h
 *   Header for intermediate ThermoPhase object for phases which
 *   employ Gibbs excess free energy based formulations
 *  (see \ref thermoprops
 * and class \link Cantera::MargulesVPSSTP MargulesVPSSTP\endlink).
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
#ifndef CT_MARGULESVPSSTP_H
#define CT_MARGULESVPSSTP_H

#include "GibbsExcessVPSSTP.h"

namespace Cantera
{

/**
 * @ingroup thermoprops
 */

//!  MargulesVPSSTP is a derived class of GibbsExcessVPSSTP that employs
//!  the Margules approximation for the excess Gibbs free energy
/*!
 *
 * MargulesVPSSTP derives from class GibbsExcessVPSSTP which is derived
 * from VPStandardStateTP,
 * and overloads the virtual methods defined there with ones that
 * use expressions appropriate for the Margules Excess Gibbs free energy
 * approximation.
 *
 * The independent unknowns are pressure, temperature, and mass fraction.
 *
 *  Several concepts are introduced. The first concept is there are temporary
 *  variables for holding the species standard state values
 *  of Cp, H, S, G, and V at the
 *  last temperature and pressure called. These functions are not recalculated
 *  if a new call is made using the previous temperature and pressure. Currently,
 *  these variables and the calculation method are handled by the VPSSMgr class,
 *  for which VPStandardStateTP owns a pointer to.
 *
 *  To support the above functionality, pressure and temperature variables,
 *  m_plast_ss and m_tlast_ss, are kept which store the last pressure and temperature
 *  used in the evaluation of standard state properties.
 *
 *  This class is usually used for nearly incompressible phases. For those phases, it
 *  makes sense to change the equation of state independent variable from
 *  density to pressure. The variable m_Pcurrent contains the current value of the
 *  pressure within the phase.
 *
 *
 * <HR>
 * <H2> Specification of Species Standard State Properties </H2>
 * <HR>
 *
 *  All species are defined to have standard states that depend upon both
 *  the temperature and the pressure. The Margules approximation assumes
 *  symmetric standard states, where all of the standard state assume
 *  that the species are in pure component states at the temperature
 *  and pressure of the solution.  I don't think it prevents, however,
 *  some species from being dilute in the solution.
 *
 *
 * <HR>
 * <H2> Specification of Solution Thermodynamic Properties </H2>
 * <HR>
 *
 * The molar excess Gibbs free energy is given by the following formula which is a sum over interactions <I>i</I>.
 * Each of the interactions are binary interactions involving two of the species in the phase, denoted, <I>Ai</I>
 * and <I>Bi</I>.
 * This is the generalization of the Margules formulation for a phase
 * that has more than 2 species.
 *
 *      \f[
 *          G^E = \sum_i \left(  H_{Ei} - T S_{Ei} \right)
 *      \f]
 *      \f[
 *         H^E_i = n X_{Ai} X_{Bi} \left( h_{o,i} +  h_{1,i} X_{Bi} \right)
 *      \f]
 *      \f[
 *         S^E_i = n X_{Ai} X_{Bi} \left( s_{o,i} +  s_{1,i} X_{Bi} \right)
 *      \f]
 *
 * where n is the total moles in the solution.
 *
 * The activity of a species defined in the phase is given by an excess
 * Gibbs free energy formulation.
 *
 *       \f[
 *            a_k = \gamma_k  X_k
 *       \f]
 *
 * where
 *
 *      \f[
 *           R T \ln( \gamma_k )= \frac{d(n G^E)}{d(n_k)}\Bigg|_{n_i}
 *      \f]
 *
 * Taking the derivatives results in the following expression
 *
 *      \f[
 *           R T \ln( \gamma_k )= \sum_i \left( \left( \delta_{Ai,k} X_{Bi} + \delta_{Bi,k} X_{Ai}  - X_{Ai} X_{Bi} \right)
 *            \left( g^E_{o,i} +  g^E_{1,i} X_{Bi} \right) +
 *            \left( \delta_{Bi,k} - X_{Bi} \right)      X_{Ai} X_{Bi}  g^E_{1,i} \right)
 *      \f]
 * where
 *        \f$  g^E_{o,i} =  h_{o,i} - T s_{o,i} \f$ and \f$ g^E_{1,i} =  h_{1,i} - T s_{1,i} \f$
 * and where \f$ X_k \f$ is the mole fraction of species <I>k</I>.
 *
 * This object inherits from the class VPStandardStateTP. Therefore, the specification and
 * calculation of all standard state and reference state values are handled at that level. Various functional
 * forms for the standard state are permissible.
 * The chemical potential for species <I>k</I> is equal to
 *
 *       \f[
 *            \mu_k(T,P) = \mu^o_k(T, P) + R T \ln(\gamma_k X_k)
 *       \f]
 *
 * The partial molar entropy for species <I>k</I> is given by the following relation,
 *
 *       \f[
 *             \tilde{s}_k(T,P) =  s^o_k(T,P)  - R \ln( \gamma_k X_k )
 *                    - R T \frac{d \ln(\gamma_k) }{dT}
 *       \f]
 *
 * The partial molar enthalpy for species <I>k</I> is given by
 *
 *       \f[
 *            \tilde{h}_k(T,P) = h^o_k(T,P) - R T^2 \frac{d \ln(\gamma_k)}{dT}
 *       \f]
 *
 * The partial molar volume for  species <I>k</I> is
 *
 *       \f[
 *              \tilde V_k(T,P)  = V^o_k(T,P)  + R T \frac{d \ln(\gamma_k) }{dP}
 *       \f]
 *
 * The partial molar Heat Capacity for species <I>k</I> is
 *
 *       \f[
 *            \tilde{C}_{p,k}(T,P) = C^o_{p,k}(T,P)   - 2 R T \frac{d \ln( \gamma_k )}{dT}
 *                    - R T^2 \frac{d^2 \ln(\gamma_k) }{{dT}^2}
 *       \f]
 *
 * <HR>
 * <H2> %Application within Kinetics Managers </H2>
 * <HR>
 *
 *   \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
 *   C^s_k, \f$ where \f$ C^s_k \f$ is a standard concentration
 *   defined below and \f$ a_k \f$ are activities used in the
 *   thermodynamic functions.  These activity (or generalized)
 *   concentrations are used
 *   by kinetics manager classes to compute the forward and
 *   reverse rates of elementary reactions.
 *   The activity concentration,\f$  C^a_k \f$,is given by the following expression.
 *
 *       \f[
 *            C^a_k = C^s_k  X_k  = \frac{P}{R T} X_k
 *       \f]
 *
 * The standard concentration for species <I>k</I> is independent of <I>k</I> and equal to
 *
 *        \f[
 *            C^s_k =  C^s = \frac{P}{R T}
 *        \f]
 *
 * For example, a bulk-phase binary gas reaction between species j and k, producing
 * a new gas species l would have the
 * following equation for its rate of progress variable, \f$ R^1 \f$, which has
 * units of kmol m-3 s-1.
 *
 *   \f[
 *    R^1 = k^1 C_j^a C_k^a =  k^1 (C^s a_j) (C^s a_k)
 *   \f]
 *  where
 *   \f[
 *      C_j^a = C^s a_j \mbox{\quad and \quad} C_k^a = C^s a_k
 *   \f]
 *
 *  \f$ C_j^a \f$ is the activity concentration of species j, and
 *  \f$ C_k^a \f$ is the activity concentration of species k. \f$ C^s \f$
 *  is the standard concentration. \f$ a_j \f$ is
 *  the activity of species j which is equal to the mole fraction of j.
 *
 *  The reverse rate constant can then be obtained from the law of microscopic reversibility
 *  and the equilibrium expression for the system.
 *
 *   \f[
 *         \frac{a_j a_k}{ a_l} = K_a^{o,1} = \exp(\frac{\mu^o_l - \mu^o_j - \mu^o_k}{R T} )
 *   \f]
 *
 *  \f$  K_a^{o,1} \f$ is the dimensionless form of the equilibrium constant, associated with
 *  the pressure dependent standard states \f$ \mu^o_l(T,P) \f$ and their associated activities,
 *  \f$ a_l \f$, repeated here:
 *
 *       \f[
 *            \mu_l(T,P) = \mu^o_l(T, P) + R T \log(a_l)
 *       \f]
 *
 *  We can switch over to expressing the equilibrium constant in terms of the reference
 *  state chemical potentials
 *
 *   \f[
 *       K_a^{o,1} = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{P}
 *   \f]
 *
 *   The concentration equilibrium constant, \f$ K_c \f$, may be obtained by changing over
 *   to activity concentrations. When this is done:
 *
 *   \f[
 *         \frac{C^a_j C^a_k}{ C^a_l} = C^o K_a^{o,1} = K_c^1 =
 *             \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} ) * \frac{P_{ref}}{RT}
 *   \f]
 *
 *  %Kinetics managers will calculate the concentration equilibrium constant, \f$ K_c \f$,
 *  using the second and third part of the above expression as a definition for the concentration
 *  equilibrium constant.
 *
 *  For completeness, the pressure equilibrium constant may be obtained as well
 *
 *   \f[
 *         \frac{P_j P_k}{ P_l P_{ref}} = K_p^1 = \exp(\frac{\mu^{ref}_l - \mu^{ref}_j - \mu^{ref}_k}{R T} )
 *   \f]
 *
 *   \f$ K_p \f$ is the simplest form of the equilibrium constant for ideal gases. However, it isn't
 *   necessarily the simplest form of the equilibrium constant for other types of phases; \f$ K_c \f$ is
 *   used instead because it is completely general.
 *
 *   The reverse rate of progress may be written down as
 *   \f[
 *    R^{-1} = k^{-1} C_l^a =  k^{-1} (C^o a_l)
 *   \f]
 *
 *  where we can use the concept of microscopic reversibility to
 *  write the reverse rate constant in terms of the
 *  forward reate constant and the concentration equilibrium
 *   constant, \f$ K_c \f$.
 *
 *    \f[
 *       k^{-1} =  k^1 K^1_c
 *    \f]
 *
 *  \f$k^{-1} \f$ has units of s-1.
 *
 * @ingroup thermoprops
 */
class MargulesVPSSTP : public GibbsExcessVPSSTP
{

public:
    //! Constructor
    /*!
     * This doesn't do much more than initialize constants with
     * default values for water at 25C. Water molecular weight
     * comes from the default elements.xml file. It actually
     * differs slightly from the IAPWS95 value of 18.015268. However,
     * density conservation and therefore element conservation
     * is the more important principle to follow.
     */
    MargulesVPSSTP();

    //! Construct and initialize a MargulesVPSSTP ThermoPhase object
    //! directly from an XML input file
    /*!
     * Working constructors
     *
     *  The two constructors below are the normal way
     *  the phase initializes itself. They are shells that call
     *  the routine initThermo(), with a reference to the
     *  XML database to get the info for the phase.
     *
     * @param inputFile Name of the input file containing the phase XML data
     *                  to set up the object
     * @param id        ID of the phase in the input file. Defaults to the
     *                  empty string.
     */
    MargulesVPSSTP(const std::string& inputFile, const std::string& id = "");

    //! Construct and initialize a MargulesVPSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML phase node containing the description of the phase
     *  @param id     id attribute containing the name of the phase.
     *                (default is the empty string)
     */
    MargulesVPSSTP(XML_Node& phaseRef, const std::string& id = "");

    //! Special constructor for a hard-coded problem
    /*!
     *  @param testProb Hard-coded value. Only the value of 1 is used. It's
     *                  for a LiKCl system test to predict the eutectic and
     *                  liquidus correctly.
     *  @deprecated To be removed after Cantera 2.2.
     */
    MargulesVPSSTP(int testProb);

    //! Copy constructor
    /*!
     *  Note this stuff will not work until the underlying phase
     *  has a working copy constructor
     *
     * @param b class to be copied
     */
    MargulesVPSSTP(const MargulesVPSSTP& b);

    //! Assignment operator
    /*!
     * @param b class to be copied.
     */
    MargulesVPSSTP& operator=(const MargulesVPSSTP& b);

    //! Duplication routine for objects which inherit from  ThermoPhase.
    /*!
     *  This virtual routine can be used to duplicate ThermoPhase objects
     *  inherited from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    //! @name  Molar Thermodynamic Properties
    //! @{

    /// Molar enthalpy. Units: J/kmol.
    virtual doublereal enthalpy_mole() const;

    /// Molar entropy. Units: J/kmol.
    virtual doublereal entropy_mole() const;

    /// Molar heat capacity at constant pressure. Units: J/kmol/K.
    virtual doublereal cp_mole() const;

    /// Molar heat capacity at constant volume. Units: J/kmol/K.
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

    //! Get the species electrochemical potentials.
    /*!
     * These are partial molar quantities.
     * This method adds a term \f$ Fz_k \phi_k \f$ to the
     * to each chemical potential.
     *
     * Units: J/kmol
     *
     * @param mu     output vector containing the species electrochemical potentials.
     *               Length: m_kk., units = J/kmol
     */
    void getElectrochemPotentials(doublereal* mu) const;

    //! Get the array of temperature second derivatives of the log activity coefficients
    /*!
     * This function is a virtual class, but it first appears in GibbsExcessVPSSTP
     * class and derived classes from GibbsExcessVPSSTP.
     *
     *  units = 1/Kelvin
     *
     * @param d2lnActCoeffdT2  Output vector of temperature 2nd derivatives of the
     *                         log Activity Coefficients. length = m_kk
     *
     */
    virtual void getd2lnActCoeffdT2(doublereal* d2lnActCoeffdT2) const;

    //! Get the array of temperature derivatives of the log activity coefficients
    /*!
     * This function is a virtual class, but it first appears in GibbsExcessVPSSTP
     * class and derived classes from GibbsExcessVPSSTP.
     *
     *  units = 1/Kelvin
     *
     * @param dlnActCoeffdT    Output vector of temperature derivatives of the
     *                         log Activity Coefficients. length = m_kk
     *
     */
    virtual void getdlnActCoeffdT(doublereal* dlnActCoeffdT) const;

    /// @}
    /// @name Initialization
    /// The following methods are used in the process of constructing
    /// the phase and setting its parameters from a specification in an
    /// input file. They are not normally used in application programs.
    /// To see how they are used, see importPhase()
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
     * from function importPhase().
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
    //! @name  Derivatives of Thermodynamic Variables needed for Applications
    //! @{

    //! Get the change in activity coefficients w.r.t. change in state (temp, mole fraction, etc.) along
    //! a line in parameter space or along a line in physical space
    /*!
     *
     * @param dTds           Input of temperature change along the path
     * @param dXds           Input vector of changes in mole fraction along the path. length = m_kk
     *                       Along the path length it must be the case that the mole fractions sum to one.
     * @param dlnActCoeffds  Output vector of the directional derivatives of the
     *                       log Activity Coefficients along the path. length = m_kk
     *  units are 1/units(s). if s is a physical coordinate then the units are 1/m.
     */
    virtual void getdlnActCoeffds(const doublereal dTds, const doublereal* const dXds, doublereal* dlnActCoeffds) const;

    //! Get the array of log concentration-like derivatives of the
    //! log activity coefficients - diagonal component
    /*!
     * This function is a virtual method.  For ideal mixtures
     * (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the
     * logarithm of the activity coefficient with respect to the
     * logarithm of the mole fraction.
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnX_diag    Output vector of the diagonal component of the log(mole fraction)
     *                 derivatives of the log Activity Coefficients.
     *                 length = m_kk
     */
    virtual void getdlnActCoeffdlnX_diag(doublereal* dlnActCoeffdlnX_diag) const;

    //! Get the array of  derivatives of the log activity coefficients wrt mole numbers - diagonal only
    /*!
     * This function is a virtual method.  For ideal mixtures
     * (unity activity coefficients), this can return zero.
     * Implementations should take the derivative of the
     * logarithm of the activity coefficient with respect to the
     * logarithm of the concentration-like variable (i.e. mole fraction,
     * molality, etc.) that represents the standard state.
     *
     *  units = dimensionless
     *
     * @param dlnActCoeffdlnN_diag    Output vector of the diagonal entries for the log(mole fraction)
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

    //@}

private:
    //! Process an XML node called "binaryNeutralSpeciesParameters"
    /*!
     * This node contains all of the parameters necessary to describe
     * the Margules model for a particular binary interaction.
     * This function reads the XML file and writes the coefficients
     * it finds to an internal data structures.
     *
     * @param xmlBinarySpecies  Reference to the XML_Node named "binaryNeutralSpeciesParameters"
     *                          containing the binary interaction
     */
    void readXMLBinarySpecies(XML_Node& xmlBinarySpecies);

    //! Resize internal arrays within the object that depend upon the number
    //! of binary Margules interaction terms
    /*!
     *  @param num Number of binary Margules interaction terms
     */
    void resizeNumInteractions(const size_t num);

    //! Initialize lengths of local variables after all species have
    //! been identified.
    void initLengths();

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

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(mole fraction)
    /*!
     * This function will be called to update the internally stored
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the mole fractions.
     */
    void s_update_dlnActCoeff_dlnX_diag() const;

    //! Update the derivative of the log of the activity coefficients
    //!  wrt log(moles) - diagonal only
    /*!
     * This function will be called to update the internally stored diagonal entries for the
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the moles.
     */
    void s_update_dlnActCoeff_dlnN_diag() const;

    //! Update the derivative of the log of the activity coefficients  wrt log(moles_m)
    /*!
     * This function will be called to update the internally stored
     * derivative of the natural logarithm of the activity coefficients
     * wrt logarithm of the mole number of species
     */
    void s_update_dlnActCoeff_dlnN() const;

protected:
    //! number of binary interaction expressions
    size_t numBinaryInteractions_;

    //! Enthalpy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_HE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_HE_c_ij;

    //! Enthalpy term for the quaternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_HE_d_ij;

    //! Entropy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_SE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_SE_c_ij;

    //! Entropy term for the quaternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_SE_d_ij;

    //! Enthalpy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_VHE_b_ij;

    //! Enthalpy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_VHE_c_ij;

    //! Enthalpy term for the quaternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_VHE_d_ij;

    //! Entropy term for the binary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_VSE_b_ij;

    //! Entropy term for the ternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_VSE_c_ij;

    //! Entropy term for the quaternary mole fraction interaction of the
    //! excess Gibbs free energy expression
    mutable vector_fp m_VSE_d_ij;

    //! vector of species indices representing species A in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and B.
     *  This vector identifies species A.
     */
    std::vector<size_t> m_pSpecies_A_ij;

    //! vector of species indices representing species B in the interaction
    /*!
     *  Each Margules excess Gibbs free energy term involves two species, A and B.
     *  This vector identifies species B.
     */
    std::vector<size_t> m_pSpecies_B_ij;

    //! form of the Margules interaction expression
    /*!
     *  Currently there is only one form.
     */
    int formMargules_;

    //! form of the temperature dependence of the Margules interaction expression
    /*!
     *  Currently there is only one form -> constant wrt temperature.
     */
    int formTempModel_;
};

}

#endif
