/**
 * @file StoichSubstanceSSTP.h
 * Header file for the StoichSubstanceSSTP class, which represents a fixed-composition
 * incompressible substance (see \ref thermoprops and
 * class \link Cantera::StoichSubstanceSSTP StoichSubstanceSSTP\endlink)
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_STOICHSUBSTANCESSTP_H
#define CT_STOICHSUBSTANCESSTP_H

#include "SingleSpeciesTP.h"

namespace Cantera
{

//!  Class StoichSubstanceSSTP represents a stoichiometric (fixed
//!   composition)  incompressible substance.
/*!
 * This class internally changes the independent degree of freedom from
 * density to pressure. This is necessary because the phase is
 * incompressible. It uses a constant volume approximation.
 *
 * <b> Specification of Species Standard State Properties </b>
 *
 *  This class inherits from SingleSpeciesTP.
 *  It is assumed that the reference state thermodynamics may be
 *  obtained by a pointer to a populated species thermodynamic property
 *  manager class (see ThermoPhase::m_spthermo). How to relate pressure
 *  changes to the reference state thermodynamics is resolved at this level.
 *
 *  For an incompressible,
 * stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties
 * are specified by giving the standard-state enthalpy, the
 * term \f$ P_0 \hat v\f$ is subtracted from the specified molar
 * enthalpy to compute the molar internal energy. The entropy is
 * assumed to be independent of the pressure.
 *
 * The enthalpy function is given by the following relation.
 *
 *       \f[
 *              h^o_k(T,P) =
 *                  h^{ref}_k(T) + \tilde v \left( P - P_{ref} \right)
 *       \f]
 *
 * For an incompressible,
 * stoichiometric substance, the molar internal energy is
 * independent of pressure. Since the thermodynamic properties
 * are specified by giving the standard-state enthalpy, the
 * term \f$ P_{ref} \tilde v\f$ is subtracted from the specified reference molar
 * enthalpy to compute the molar internal energy.
 *
 *       \f[
 *            u^o_k(T,P) = h^{ref}_k(T) - P_{ref} \tilde v
 *       \f]
 *
 * The standard state heat capacity and entropy are independent
 * of pressure. The standard state Gibbs free energy is obtained
 * from the enthalpy and entropy functions.
 *
 * <b> Specification of Solution Thermodynamic Properties </b>
 *
 *  All solution properties are obtained from the standard state
 *  species functions, since there is only one species in the phase.
 *
 * <b> Application within Kinetics Managers </b>
 *
 * The standard concentration is equal to 1.0. This means that the
 * kinetics operator works on an (activities basis). Since this
 * is a stoichiometric substance, this means that the concentration
 * of this phase drops out of kinetics expressions.
 *
 * An example of a reaction using this is a sticking coefficient
 * reaction of a substance in an ideal gas phase on a surface with a bulk phase
 * species in this phase. In this case, the rate of progress for this
 * reaction, \f$ R_s \f$, may be expressed via the following equation:
 *   \f[
 *    R_s = k_s C_{gas}
 *   \f]
 * where the units for \f$ R_s \f$ are kmol m-2 s-1. \f$ C_{gas} \f$ has units
 * of kmol m-3. Therefore, the kinetic rate constant,  \f$ k_s \f$, has
 * units of m s-1. Nowhere does the concentration of the bulk phase
 * appear in the rate constant expression, since it's a stoichiometric
 * phase and the activity is always equal to 1.0.
 *
 * <b> Instantiation of the Class </b>
 *
 * The constructor for this phase is NOT located in the default ThermoFactory
 * for %Cantera. However, a new StoichSubstanceSSTP may be created by
 * the following code snippets:
 *
 * @code
 *    sprintf(file_ID,"%s#NaCl(S)", iFile);
 *    XML_Node *xm = get_XML_NameID("phase", file_ID, 0);
 *    StoichSubstanceSSTP *solid = new StoichSubstanceSSTP(*xm);
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    sprintf(file_ID,"%s#NaCl(S)", iFile);
 *    XML_Node *xm = get_XML_NameID("phase", file_ID, 0);
 *    StoichSubstanceSSTP solid;
 *    importPhase(*xm, &solid);
 * @endcode
 *
 *   <b> XML Example </b>
 *
 * The phase model name for this is called StoichSubstance. It must be supplied
 * as the model attribute of the thermo XML element entry.
 * Within the phase XML block,
 * the density of the phase must be specified. An example of an XML file
 * this phase is given below.
 *
 * @code
 * <!-- phase NaCl(S)    -->
 * <phase dim="3" id="NaCl(S)">
 *    <elementArray datasrc="elements.xml">
 *       Na Cl
 *    </elementArray>
 *    <speciesArray datasrc="#species_NaCl(S)"> NaCl(S) </speciesArray>
 *    <thermo model="StoichSubstanceSSTP">
 *       <density units="g/cm3">2.165</density>
 *    </thermo>
 *    <transport model="None"/>
 *    <kinetics model="none"/>
 * </phase>
 *
 * <!-- species definitions     -->
 * <speciesData id="species_NaCl(S)">
 *   <!-- species NaCl(S)   -->
 *   <species name="NaCl(S)">
 *      <atomArray> Na:1 Cl:1 </atomArray>
 *      <thermo>
 *         <Shomate Pref="1 bar" Tmax="1075.0" Tmin="250.0">
 *            <floatArray size="7">
 *                50.72389, 6.672267, -2.517167,
 *                10.15934, -0.200675, -427.2115,
 *                130.3973
 *            </floatArray>
 *         </Shomate>
 *      </thermo>
 *      <density units="g/cm3">2.165</density>
 *    </species>
 * </speciesData>  @endcode
 *
 *  The model attribute, "StoichSubstanceSSTP", on the thermo element
 *  identifies the phase as being a StoichSubstanceSSTP object.
 *
 * @ingroup thermoprops
 */
class StoichSubstanceSSTP : public SingleSpeciesTP
{
public:
    //! Default constructor for the StoichSubstanceSSTP class
    StoichSubstanceSSTP() {}

    //! Construct and initialize a StoichSubstanceSSTP ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    StoichSubstanceSSTP(const std::string& infile, std::string id = "");

    //! Construct and initialize a StoichSubstanceSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML node pointing to a StoichSubstanceSSTP description
     *  @param id       Id of the phase.
     */
    StoichSubstanceSSTP(XML_Node& phaseRef, const std::string& id = "");

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    StoichSubstanceSSTP(const StoichSubstanceSSTP&  right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    StoichSubstanceSSTP& operator=(const StoichSubstanceSSTP& right);

    //! Duplication function
    /*!
     * This virtual function is used to create a duplicate of the
     * current phase. It's used to duplicate the phase when given
     * a ThermoPhase pointer to the phase.
     *
     * @return It returns a ThermoPhase pointer.
     */
    ThermoPhase* duplMyselfAsThermoPhase() const;

    /**
     * Equation of state flag.
     *
     * Returns the value cStoichSubstance, defined in mix_defs.h.
     */
    virtual int eosType() const;

    //! @name Mechanical Equation of State
    //! @{

    //! Report the Pressure. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent
     * of pressure. This method simply returns the stored
     * pressure value.
     */
    virtual doublereal pressure() const;

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     * For an incompressible substance, the density is
     * independent of pressure. Therefore, this method only
     * stores the specified pressure value. It does not
     * modify the density.
     *
     * @param p Pressure (units - Pa)
     */
    virtual void setPressure(doublereal p);

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     */
    virtual doublereal isothermalCompressibility() const;

    //! Return the volumetric thermal expansion coefficient. Units: 1/K.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal thermalExpansionCoeff() const ;

    /**
     * @}
     * @name Activities, Standard States, and Activity Concentrations
     *
     *  This section is largely handled by parent classes, since there
     *  is only one species. Therefore, the activity is equal to one.
     * @{
     */

    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k /
     * C^0_k, \f$ where \f$ C^0_k \f$ is a standard concentration
     * defined below and \f$ a_k \f$ are activities used in the
     * thermodynamic functions.  These activity (or generalized)
     * concentrations are used
     * by kinetics manager classes to compute the forward and
     * reverse rates of elementary reactions.
     *
     *  For a stoichiometric substance, there is
     *  only one species, and the generalized concentration is 1.0.
     *
     * @param c Output array of generalized concentrations. The
     *           units depend upon the implementation of the
     *           reaction rate expressions within the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize
     * the activity (i.e., generalized) concentration.
     * This phase assumes that the kinetics operator works on an
     * dimensionless basis. Thus, the standard concentration is
     * equal to 1.0.
     *
     * @param k Optional parameter indicating the species. The default
     *         is to assume this refers to species 0.
     * @return
     *   Returns The standard Concentration as 1.0
     */
    virtual doublereal standardConcentration(size_t k=0) const;

    //! Natural logarithm of the standard concentration of the kth species.
    /*!
     * @param k    index of the species (defaults to zero)
     */
    virtual doublereal logStandardConc(size_t k=0) const;

    //! Get the array of chemical potentials at unit activity for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * For a stoichiometric substance, there is no activity term in
     * the chemical potential expression, and therefore the
     * standard chemical potential and the chemical potential
     * are both equal to the molar Gibbs function.
     *
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current
     * temperature and pressure of the solution
     *
     * @param mu0     Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

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
     *
     *     uA[0] = kmol units - default  = 1
     *     uA[1] = m    units - default  = -nDim(), the number of spatial
     *                                   dimensions in the Phase class.
     *     uA[2] = kg   units - default  = 0;
     *     uA[3] = Pa(pressure) units - default = 0;
     *     uA[4] = Temperature units - default = 0;
     *     uA[5] = time units - default = 0
     *
     * @param k species index. Defaults to 0.
     * @param sizeUA output int containing the size of the vector.
     *        Currently, this is equal to 6.
     * @deprecated To be removed after Cantera 2.2.
     */
    virtual void getUnitsStandardConc(doublereal* uA, int k = 0,
                                      int sizeUA = 6) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

    //! Get the nondimensional Enthalpy functions for the species
    //! at their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt      Output vector of  nondimensional standard state enthalpies.
     *                 Length: m_kk.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the array of nondimensional Entropy functions for the
    //! standard state species at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param sr   Output vector of  nondimensional standard state entropies.
     *             Length: m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Get the nondimensional Gibbs functions for the species
    //! in their standard states at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param grt  Output vector of nondimensional standard state Gibbs free energies
     *             Length: m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the nondimensional Heat Capacities at constant
    //! pressure for the species standard states
    //! at the current <I>T</I> and <I>P</I> of the solution
    /*!
     * @param cpr   Output vector of nondimensional standard state heat capacities
     *              Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const;

    //!  Returns the vector of nondimensional Internal Energies  of the standard
    //!  state species at the current <I>T</I> and <I>P</I> of the solution
    /*!
     *  For an incompressible,
     * stoichiometric substance, the molar internal energy is
     * independent of pressure. Since the thermodynamic properties
     * are specified by giving the standard-state enthalpy, the
     * term \f$ P_{ref} \hat v\f$ is subtracted from the specified reference molar
     * enthalpy to compute the standard state molar internal energy.
     *
     * @param urt  output vector of nondimensional standard state
     *             internal energies of the species. Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    //! Returns the vector of nondimensional
    //!  internal Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state
     *               internal energies of the species.
     *               Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(doublereal* urt) const;
    // @}

    /*
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

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Set the equation of state parameters
    /*!
     * @internal
     *  The number and meaning of these depends on the subclass.
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     *        c[0] = density of phase [ kg/m3 ]
     */
    virtual void setParameters(int n, doublereal* const c);

    //! Get the equation of state parameters in a vector
    /*!
     * @internal
     *
     * @param n number of parameters
     * @param c array of \a n coefficients
     *
     *  For this phase:
     *       -  n = 1
     *       -  c[0] = density of phase [ kg/m3 ]
     */
    virtual void getParameters(int& n, doublereal* const c) const;

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() when processing a phase
     * definition in an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialized with elements and/or species.
     *
     *  For this phase, the density of the phase is specified in this block.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     *
     * eosdata points to the thermo block, and looks like this:
     *
     *   @code
     *   <phase id="stoichsolid" >
     *     <thermo model="StoichSubstance">
     *         <density units="g/cm3">3.52</density>
     *     </thermo>
     *   </phase>
     *   @endcode
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);
};

//!  Class %electrodeElectron represents an electron in a
//!  metal using the Standard hydrogen reference electrode
/*!
 *  The class is based on the electron have a chemical potential
 *  equal to one-half of the entropy of the H2 gas at 1 bar.
 *
 *  @deprecated Deprecated in favor of class MetalSHEelectrons. To be removed
 *      after Cantera 2.2.
 */
class electrodeElectron : public StoichSubstanceSSTP
{
public:
    //! Default constructor for the electrodeElectron class
    electrodeElectron();

    //! Construct and initialize a electrodeElectron ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    electrodeElectron(const std::string& infile, std::string id = "");

    //! Construct and initialize a electrodeElectron ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML node pointing to a electrodeElectron description
     *  @param id       Id of the phase.
     */
    electrodeElectron(XML_Node& phaseRef, const std::string& id = "");

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    electrodeElectron(const electrodeElectron&  right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    electrodeElectron& operator=(const electrodeElectron& right);

    //! Destructor
    virtual ~electrodeElectron() {}

    void setParametersFromXML(const XML_Node& eosdata);

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    void setParameters(int n, doublereal* const c);
};

}

#endif
