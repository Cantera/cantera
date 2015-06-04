/**
 * @file FixedChemPotSSTP.h
 * Header file for the FixedChemPotSSTP class, which represents a fixed-composition
 * incompressible substance with a constant chemical potential (see \ref thermoprops and
 * class \link Cantera::FixedChemPotSSTP FixedChemPotSSTP\endlink)
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_FIXEDCHEMPOTSSTP_H
#define CT_FIXEDCHEMPOTSSTP_H

#include "SingleSpeciesTP.h"

namespace Cantera
{

//!  Class FixedChemPotSSTP represents a stoichiometric (fixed
//!   composition)  incompressible substance.
/*!
 * This class internally changes the independent degree of freedom from
 * density to pressure. This is necessary because the phase is
 * incompressible. It uses a zero volume approximation.
 *
 * <b> Specification of Species Standard State Properties </b>
 *
 *  This class inherits from SingleSpeciesTP.
 *  It uses a single value for the chemical potential which is assumed to be constant
 *  with respect to temperature and pressure.
 *
 *  The reference state thermodynamics is inherited from SingleSpeciesTP. However,
 *  it's only used to set the initial chemical potential to the value
 *  of the chemical potential at the starting conditions. Thereafter,
 *  it is ignored.
 *
 *  For a zero volume material, the internal energy and the enthalpy are
 *  equal to the chemical potential. The entropy, the heat capacity, and the molar volume
 *  are equal to zero.
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
 * phase, and the activity is always equal to 1.0.
 *
 * <b> Instantiation of the Class </b>
 *
 * This phase may be instantiated by calling the default ThermoFactory routine
 * for %Cantera. This new FixedChemPotSSTP object must then have a standalone XML file
 * description an example of which is given below.
 *
 * It may also be created by the following code snippets. The code
 * includes the special member function setChemicalPotential( chempot), which
 * sets the chemical potential to a specific value in J / kmol.
 *
 * @code
 *    sprintf(file_ID,"%s#Li(Fixed)", iFile);
 *    XML_Node *xm = get_XML_NameID("phase", file_ID, 0);
 *    FixedChemPotSSTP *LiFixed = new FixedChemPotSSTP(*xm);
      // Set the chemical potential to -2.3E7 J/kmol
 *    LiFixed->setChemicalPotential(-2.3E7.)
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    sprintf(file_ID,"%s#NaCl(S)", iFile);
 *    XML_Node *xm = get_XML_NameID("phase", file_ID, 0);
 *    FixedChemPotSSTP solid;
 *    importPhase(*xm, &solid);
 * @endcode
 *
 * The phase may also be created by a special constructor so that element
 * potentials may be set. The constructor takes the name of the element and
 * the value of the element chemical potential. An example is given below.
 *
 * @code
 *     FixedChemPotSSTP *LiFixed = new FixedChemPotSSTP("Li", -2.3E7);
 * @endcode
 *
 *   <b> XML Example </b>
 *
 * The phase model name for this is called FixedChemPot. It must be supplied
 * as the model attribute of the thermo XML element entry.
 *
 * @code
 * <?xml version="1.0"?>
 * <ctml>
 *   <validate reactions="yes" species="yes"/>
 *
 *   <!-- phase NaCl(S)    -->
 *   <phase dim="3" id="LiFixed">
 *     <elementArray datasrc="elements.xml">
 *       Li
 *     </elementArray>
 *     <speciesArray datasrc="#species_Li(Fixed)">
 *       LiFixed
 *     </speciesArray>
 *     <thermo model="FixedChemPot">
 *       <chemicalPotential units="J/kmol"> -2.3E7  </chemicalPotential>
 *     </thermo>
 *     <transport model="None"/>
 *     <kinetics model="none"/>
 *   </phase>
 *
 *   <!-- species definitions     -->
 *   <speciesData id="species_Li(Fixed)">
 *     <species name="LiFixed">
 *       <atomArray> Li:1 </atomArray>
 *       <thermo>
 *         <Shomate Pref="1 bar" Tmax="1075.0" Tmin="250.0">
 *           <floatArray size="7">
 *             50.72389, 6.672267, -2.517167,
 *             10.15934, -0.200675, -427.2115,
 *             130.3973
 *           </floatArray>
 *         </Shomate>
 *       </thermo>
 *     </species>
 *   </speciesData>
 * </ctml>
 * @endcode
 *
 * The model attribute, "FixedChemPot", on the thermo element
 * identifies the phase as being a FixedChemPotSSTP object.
 *
 * @ingroup thermoprops
 */
class FixedChemPotSSTP : public SingleSpeciesTP
{
public:
    //! Default constructor for the FixedChemPotSSTP class
    FixedChemPotSSTP();

    //! Construct and initialize a FixedChemPotSSTP ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    FixedChemPotSSTP(const std::string& infile, std::string id = "");

    //! Construct and initialize a FixedChemPotSSTP ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML node pointing to a FixedChemPotSSTP description
     *  @param id       Id of the phase.
     */
    FixedChemPotSSTP(XML_Node& phaseRef, const std::string& id = "");

    //! Copy constructor
    /*!
     * @param right Object to be copied
     */
    FixedChemPotSSTP(const FixedChemPotSSTP&  right);

    //! Special constructor for the FixecChemPotSSTP class setting an element chemical
    //! potential directly
    /*!
     *  This will create a FixedChemPotSSTP consisting of a single species with the
     *  stoichiometry of one of the specified atom. It will have a chemical potential
     *  that is given by the second argument.
     *
     *  @param Ename String name of the element
     *  @param chemPot  Value of the chemical potential of that element (J/kmol)
     */
    FixedChemPotSSTP(const std::string& Ename, doublereal chemPot);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    FixedChemPotSSTP& operator=(const FixedChemPotSSTP& right);

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

    //! @}
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
     * @param uA Output vector containing the units:
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
    /// @name Partial Molar Properties of the Solution
    /// These properties are handled by the parent class, SingleSpeciesTP
    //@{

    //! Get the species partial molar volumes. Units: m^3/kmol.
    /*!
     * This is the phase molar volume.  \f$ V(T,P) = V_o(T,P) \f$.
     *
     * set to zero.
     *
     *  @param vbar On return, contains the molar volume of the single species
     *              and the phase. Units are m^3 / kmol. Length = 1
     */
    void getPartialMolarVolumes(doublereal* vbar) const;

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

    //! Get the molar volumes of each species in their standard
    //! states at the current  <I>T</I> and <I>P</I> of the solution.
    /*
     *   units = m^3 / kmol
     *
     * We set this to zero
     *
     * @param vbar On output this contains the standard volume of the species
     *             and phase (m^3/kmol). Vector of length 1
     */
    virtual void getStandardVolumes(doublereal* vbar) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    //! Returns the vector of nondimensional
    //!  internal Energies of the reference state at the current temperature
    //!  of the solution and the reference pressure for each species.
    /*!
     * @param urt    Output vector of nondimensional reference state internal
     *               energies of the species. Length: m_kk
     */
    virtual void getIntEnergy_RT_ref(doublereal* urt) const;

    /*!
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param hrt     Output vector containing the nondimensional reference state enthalpies
     *                Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /*!
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param grt     Output vector containing the nondimensional reference state
     *                Gibbs Free energies.  Length: m_kk.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    /*!
     *  Returns the vector of the
     *  Gibbs function of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *  units = J/kmol
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param g       Output vector containing the  reference state
     *                Gibbs Free energies.  Length: m_kk. Units: J/kmol.
     */
    virtual void  getGibbs_ref(doublereal* g) const;

    /*!
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for each species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param er      Output vector containing the nondimensional reference state
     *                entropies.  Length: m_kk.
       */
    virtual void getEntropy_R_ref(doublereal* er) const;

    /*!
     *  Returns the vector of nondimensional
     *  constant pressure heat capacities of the reference state
     *  at the current temperature of the solution
     *  and reference pressure for each species.
     *
     *  This function is resolved in this class.  It is assumed that the m_spthermo species thermo
     *  pointer is populated and yields the reference state.
     *
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    //@}

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Set the equation of state parameters
    /*!
     * @internal
     * @param n number of parameters = 1
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
     *  For this phase, the chemical potential is set
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     *
     * eosdata points to the thermo block, and looks like this:
     *
     * @code
     *   <phase id="stoichsolid" >
     *     <thermo model="FixedChemPot">
     *       <chemicalPotential units="J/kmol"> -2.7E7 </chemicalPotential>
     *     </thermo>
     *   </phase>
     * @endcode
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);

    //! Function to set the chemical potential directly
    /*!
     *  @param chemPot  Value of the chemical potential (units J/kmol)
     */
    void setChemicalPotential(doublereal chemPot);

protected:
    //!  Value of the chemical potential of the bath species
    /*!
     *  units are J/kmol
     */
    doublereal chemPot_;
};

}

#endif
