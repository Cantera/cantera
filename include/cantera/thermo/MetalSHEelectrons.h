/**
 * @file MetalSHEelectrons.h
 * Header file for the MetalSHEElectrons class, which represents the
 * electrons in a metal that are consistent with the
 * SHE electrode (see \ref thermoprops and
 * class \link Cantera::MetalSHEelectrons MetalSHEelectrons\endlink)
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_METALSHEELECTRONS_H
#define CT_METALSHEELECTRONS_H

#include "SingleSpeciesTP.h"

namespace Cantera
{

//! Class MetalSHEelectrons represents electrons within a metal, adjacent to an
//! aqueous electrolyte, that are consistent with the SHE reference electrode.
/*!
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * @deprecated To be removed after Cantera 2.4
 *
 * The class is based on the electron having a chemical potential equal to one-
 * half of the entropy of the H2 gas at the system pressure
 *
 * ## Specification of Species Standard State Properties
 *
 * This class inherits from SingleSpeciesTP. It is assumed that the reference
 * state thermodynamics may be obtained by a pointer to a populated species
 * thermodynamic property manager class (see ThermoPhase::m_spthermo). How to
 * relate pressure changes to the reference state thermodynamics is resolved at
 * this level.
 *
 * The enthalpy function is given by the following relation.
 *
 *       \f[
 *            h^o_k(T,P) = h^{ref}_k(T)
 *       \f]
 *
 * The standard state constant-pressure heat capacity is independent of pressure:
 *
 *       \f[
 *            Cp^o_k(T,P) = Cp^{ref}_k(T)
 *       \f]
 *
 * The standard state entropy depends in the following fashion on pressure:
 *
 *       \f[
 *            S^o_k(T,P) = S^{ref}_k(T) -  R \ln(\frac{P}{P_{ref}})
 *       \f]
 *
 * The standard state Gibbs free energy is obtained from the enthalpy and
 * entropy functions:
 *
 *       \f[
 *            \mu^o_k(T,P) =  h^o_k(T,P) - S^o_k(T,P) T
 *       \f]
 *
 *       \f[
 *            \mu^o_k(T,P) =  \mu^{ref}_k(T) + R T \ln( \frac{P}{P_{ref}})
 *       \f]
 *
 * where
 *       \f[
 *            \mu^{ref}_k(T) =   h^{ref}_k(T)   - T S^{ref}_k(T)
 *       \f]
 *
 * The standard state internal energy is obtained from the enthalpy function also
 *
 *       \f[
 *            u^o_k(T,P) = h^o_k(T) - R T
 *       \f]
 *
 * ## Specification of Solution Thermodynamic Properties
 *
 * All solution properties are obtained from the standard state species
 * functions, since there is only one species in the phase.
 *
 * ## %Application within Kinetics Managers
 *
 * The standard concentration is equal to 1.0. This means that the kinetics
 * operator works on an activities basis. Since this is a stoichiometric
 * substance, this means that the concentration of this phase drops out of
 * kinetics expressions since the activity is always equal to one.
 *
 * This is what is expected of electrons. The only effect that this class will
 * have on reactions is in terms of the standard state chemical potential, which
 * is equal to 1/2 of the H2 gas chemical potential, and the voltage assigned to
 * the electron, which is the voltage of the metal.
 *
 * ## Instantiation of the Class
 *
 * The constructor for this phase is located in the default ThermoFactory for
 * %Cantera. A new MetalSHEelectrons object may be created by the following code
 * snippets, where the file metalSHEelectrons.xml exists in a local directory:
 *
 * @code
 *    MetalSHEelectrons *eMetal = new MetalSHEelectrons("metalSHEelectrons.xml", "");
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", iFile + "#MetalSHEelectrons", 0);
 *    MetalSHEelectrons eMetal;
 *    importPhase(*xm, &eMetal);
 * @endcode
 *
 *  @code
 *  ThermoPhase *eMetal = newPhase("MetalSHEelectrons.xml", "MetalSHEelectrons");
 *  @endcode
 *
 * Additionally, this phase may be created without including an XML file with
 * the special command, where the default file is embedded into this object.
 *
 * @code
 *    MetalSHEelectrons *eMetal = new MetalSHEelectrons("MetalSHEelectrons_default.xml", "");
 * @endcode
 *
 * ## XML Example
 *
 * The phase model name for this is called MetalSHEelectrons. It must be
 * supplied as the model attribute of the thermo XML element entry. Within the
 * phase XML block, the density of the phase must be specified though it's not
 * used. An example of an XML file this phase is given below.
 *
 * @code
 * <?xml version="1.0"?>
 * <ctml>
 *   <validate reactions="yes" species="yes"/>
 *
 *   <phase dim="3" id="MetalSHEelectrons">
 *     <elementArray datasrc="elements.xml">
 *         E
 *     </elementArray>
 *     <speciesArray datasrc="#species_Metal_SHEelectrons"> she_electron </speciesArray>
 *     <thermo model="metalSHEelectrons">
 *       <density units="g/cm3">2.165</density>
 *     </thermo>
 *     <transport model="None"/>
 *     <kinetics model="none"/>
 *   </phase>
 *
 *   <!-- species definitions     -->
 *   <speciesData id="species_Metal_SHEelectrons">
 *     <species name="she_electron">
 *       <atomArray> E:1  </atomArray>
 *       <charge> -1 </charge>
 *       <thermo>
 *         <NASA Tmax="1000.0" Tmin="200.0" P0="100000.0">
 *          <floatArray name="coeffs" size="7">
 *           1.172165560E+00,   3.990260375E-03,  -9.739075500E-06,  1.007860470E-08,
 *          -3.688058805E-12, -4.589675865E+02,  3.415051190E-01
 *          </floatArray>
 *         </NASA>
 *         <NASA Tmax="6000.0" Tmin="1000.0" P0="100000.0">
 *            <floatArray name="coeffs" size="7">
 *              1.466432895E+00,  4.133039835E-04, -7.320116750E-08, 7.705017950E-12,
 *             -3.444022160E-16, -4.065327985E+02, -5.121644350E-01
 *           </floatArray>
 *         </NASA>
 *       </thermo>
 *       <density units="g/cm3">2.165</density>
 *     </species>
 *   </speciesData>
 * </ctml>
 * @endcode
 *
 * The model attribute, "MetalSHEelectrons", on the thermo element identifies
 * the phase as being a MetalSHEelectrons object.
 *
 * @ingroup thermoprops
 */
class MetalSHEelectrons : public SingleSpeciesTP
{
public:
    //! Default constructor for the MetalSHEelectrons class
    MetalSHEelectrons();

    //! Construct and initialize a MetalSHEelectrons ThermoPhase object
    //! directly from an ASCII input file
    /*!
     * @param infile name of the input file
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    MetalSHEelectrons(const std::string& infile, const std::string& id = "");

    //! Construct and initialize a MetalSHEelectrons ThermoPhase object
    //! directly from an XML database
    /*!
     *  @param phaseRef XML node pointing to a MetalSHEelectrons description
     *  @param id       Id of the phase.
     */
    MetalSHEelectrons(XML_Node& phaseRef, const std::string& id = "");

    //! @name Mechanical Equation of State
    //! @{

    //! Report the Pressure. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent of pressure.
     * This method simply returns the stored pressure value.
     */
    virtual doublereal pressure() const;

    //! Set the pressure at constant temperature. Units: Pa.
    /*!
     * For an incompressible substance, the density is independent of pressure.
     * Therefore, this method only stores the specified pressure value. It does
     * not modify the density.
     *
     * @param p Pressure (units - Pa)
     */
    virtual void setPressure(doublereal p);

    virtual doublereal isothermalCompressibility() const;
    virtual doublereal thermalExpansionCoeff() const;

    //! @}
    //! @name Activities, Standard States, and Activity Concentrations
    //!
    //! This section is largely handled by parent classes, since there
    //! is only one species. Therefore, the activity is equal to one.
    //! @{

    //! This method returns an array of generalized concentrations
    /*!
     * \f$ C^a_k\f$ are defined such that \f$ a_k = C^a_k / C^0_k, \f$ where
     * \f$ C^0_k \f$ is a standard concentration defined below and \f$ a_k \f$
     * are activities used in the thermodynamic functions. These activity (or
     * generalized) concentrations are used by kinetics manager classes to
     * compute the forward and reverse rates of elementary reactions.
     *
     * For a stoichiometric substance, there is only one species, and the
     * generalized concentration is 1.0.
     *
     * @param c Output array of generalized concentrations. The units depend
     *           upon the implementation of the reaction rate expressions within
     *           the phase.
     */
    virtual void getActivityConcentrations(doublereal* c) const;

    //! Return the standard concentration for the kth species
    /*!
     * The standard concentration \f$ C^0_k \f$ used to normalize the activity
     * (i.e., generalized) concentration. This phase assumes that the kinetics
     * operator works on an dimensionless basis. Thus, the standard
     * concentration is equal to 1.0.
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

    //! Get the array of chemical potentials at unit activity for the species at
    //! their standard states at the current *T* and *P* of the solution.
    /*!
     * For a stoichiometric substance, there is no activity term in the chemical
     * potential expression, and therefore the standard chemical potential and
     * the chemical potential are both equal to the molar Gibbs function.
     *
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P) \f$.
     * The values are evaluated at the current temperature and pressure of the
     * solution
     *
     * @param mu0     Output vector of chemical potentials.
     *                Length: m_kk.
     */
    virtual void getStandardChemPotentials(doublereal* mu0) const;

    //@}
    /// @name  Properties of the Standard State of the Species in the Solution
    //@{

    virtual void getEnthalpy_RT(doublereal* hrt) const;
    virtual void getEntropy_R(doublereal* sr) const;
    virtual void getGibbs_RT(doublereal* grt) const;
    virtual void getCp_R(doublereal* cpr) const;

    //! Returns the vector of nondimensional Internal Energies of the standard
    //! state species at the current *T* and *P* of the solution
    /*!
     * For an incompressible, stoichiometric substance, the molar internal
     * energy is independent of pressure. Since the thermodynamic properties are
     * specified by giving the standard-state enthalpy, the term \f$ P_{ref}
     * \hat v\f$ is subtracted from the specified reference molar enthalpy to
     * compute the standard state molar internal energy.
     *
     * @param urt  output vector of nondimensional standard state
     *             internal energies of the species. Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //@}
    /// @name Thermodynamic Values for the Species Reference States
    //@{

    virtual void getIntEnergy_RT_ref(doublereal* urt) const;
    // @}

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Make the default XML tree
    /*!
     *   @returns a new XML tree containing the default info.
     */
    static XML_Node* makeDefaultXMLTree();

    //! Set the equation of state parameters
    /*!
     * @internal
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
     * For this phase, the density of the phase is specified in this block.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     *
     * eosdata points to the thermo block, and looks like this:
     *
     * @code
     * <phase id="stoichsolid" >
     *   <thermo model="StoichSubstance">
     *     <density units="g/cm3">3.52</density>
     *   </thermo>
     * </phase>
     * @endcode
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);
};

}
#endif
