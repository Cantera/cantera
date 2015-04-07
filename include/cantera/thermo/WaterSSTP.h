/**
 *  @file WaterSSTP.h
 * Declares a %ThermoPhase class consisting of  pure water (see \ref thermoprops
 * and class \link Cantera::WaterSSTP WaterSSTP\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_WATERSSTP_H
#define CT_WATERSSTP_H

#include "SingleSpeciesTP.h"


namespace Cantera
{

class WaterPropsIAPWS;
class WaterProps;
//!  Class for single-component water. This is designed to cover just the
//!  liquid part of water.
/*!
 *  The reference is W. Wagner, A. Prub, "The IAPWS Formulation 1995 for the Thermodynamic
 *  Properties of Ordinary Water Substance for General and Scientific Use,"
 *  J. Phys. Chem. Ref. Dat, 31, 387, 2002.
 *
 * <HR>
 * <H2> Specification of Species Standard %State Properties </H2>
 * <HR>
 *
 *   The offsets used in the steam tables are different than NIST's.
 *   They assume u_liq(TP) = 0.0, s_liq(TP) = 0.0, where TP is the
 *   triple point conditions:
 *
 *      -  u(273.16, rho)    = 0.0
 *      -  s(273.16, rho)    = 0.0
 *      -  psat(273.16)      = 611.655 Pascal
 *      -  rho(273.16, psat) = 999.793 kg m-3
 *
 *   These "steam table" assumptions are used by the WaterPropsIAPWS class.
 *   Therefore, offsets must be calculated to make the thermodynamic
 *   properties calculated within this class to be consistent with
 *   thermo properties within Cantera.
 *
 *   The thermodynamic base state for water is set to the NIST basis here
 *   by specifying constants, #EW_Offset and #SW_Offset, one for energy
 *   quantities and one for entropy quantities. The offsets are
 *   specified so that the following properties hold:
 *
 *   - Delta_Hfo_idealgas(298.15) = -241.826 kJ/gmol
 *   - So_idealgas(298.15, 1bar)  =  188.835  J/gmolK
 *
 *   (From http://webbook.nist.gov)
 *
 *   The "o" here refers to a hypothetical ideal gas state. The way
 *   we achieve this in practice is to evaluate at a very low pressure
 *   and then use the theoretical ideal gas results to scale up to
 *   higher pressures:
 *
 *   Ho(1bar) = H(P0)
 *
 *   So(1bar) = S(P0) + RT ln(1bar/P0)
 *
 * <HR>
 * <H2> %Application within %Kinetics Managers </H2>
 * <HR>
 *
 *   This is unimplemented.
 *
 * <HR>
 * <H2> Instantiation of the Class </H2>
 * <HR>
 *
 * The constructor for this phase is NOT located in the default ThermoFactory
 * for %Cantera. However, a new %WaterSSTP object may be created by
 * the following code snippets, combined with an XML file given in the
 * XML example section.
 *
 * @code
 *      WaterSSTP *w = new WaterSSTP("waterSSTPphase.xml","");
 * @endcode
 *
 * or
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "waterSSTPphase.xml#water", 0);
 *    WaterSSTP *w = new WaterSSTP(*xm);
 * @endcode
 *
 * or by the following call to importPhase():
 *
 * @code
 *    XML_Node *xm = get_XML_NameID("phase", "waterSSTPphase.xml#water", 0);
 *    WaterSSTP water;
 *    importPhase(*xm, &water);
 * @endcode
 *
 * <HR>
 * <H2> XML Example </H2>
 * <HR>
 *
 *   An example of an XML Element named phase setting up a WaterSSTP object with
 *   id "water"   is given below.
 *
 * @code
 * <!-- phase water     -->
 * <phase dim="3" id="water">
 *   <elementArray datasrc="elements.xml">O  H </elementArray>
 *   <speciesArray datasrc="#species_data">H2O</speciesArray>
 *   <state>
 *     <temperature units="K">300.0</temperature>
 *     <pressure units="Pa">101325.0</pressure>
 *   </state>
 *   <thermo model="PureLiquidWater"/>
 *   <kinetics model="none"/>
 * </phase>
 * @endcode
 *
 *  Note the model "PureLiquidWater" indicates the usage of the WaterSSTP object.
 *
 * @ingroup thermoprops
 */
class WaterSSTP : public SingleSpeciesTP
{
public:
    //! Base constructor
    WaterSSTP();

    //! Copy constructor
    WaterSSTP(const WaterSSTP&);

    //! Assignment operator
    WaterSSTP& operator=(const WaterSSTP&);

    //! Full constructor for a water phase
    /*!
     * @param inputFile String name of the input file
     * @param id        string id of the phase name
     */
    explicit WaterSSTP(const std::string& inputFile, const std::string& id = "");

    //! Full constructor for a water phase
    /*!
     * @param phaseRef  XML node referencing the water phase.
     * @param id        string id of the phase name
     */
    explicit WaterSSTP(XML_Node& phaseRef, const std::string& id = "");

    //! Destructor
    virtual ~WaterSSTP();

    //! Duplicator from a ThermoPhase object
    ThermoPhase* duplMyselfAsThermoPhase() const;

    virtual int eosType() const {
        return -1;
    }

    //! @name  Molar Thermodynamic Properties of the Solution
    //! @{

    virtual doublereal cv_mole() const;

    //@}
    /// @name Mechanical Equation of State Properties
    //@{

    virtual doublereal pressure() const;
    virtual void setPressure(doublereal p);

    //! Returns  the isothermal compressibility. Units: 1/Pa.
    /*!
     * The isothermal compressibility is defined as
     * \f[
     * \kappa_T = -\frac{1}{v}\left(\frac{\partial v}{\partial P}\right)_T
     * \f]
     *  or
     * \f[
     * \kappa_T = \frac{1}{\rho}\left(\frac{\partial \rho}{\partial P}\right)_T
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
    virtual doublereal thermalExpansionCoeff() const;

    //! Return the derivative of the volumetric thermal expansion coefficient. Units: 1/K2.
    /*!
     * The thermal expansion coefficient is defined as
     * \f[
     * \beta = \frac{1}{v}\left(\frac{\partial v}{\partial T}\right)_P
     * \f]
     */
    virtual doublereal dthermalExpansionCoeffdT() const;

    //! @}
    //! @name Properties of the Standard State of the Species in the Solution
    //! @{

    //!  Get the gibbs function for the species
    //!  standard states at the current T and P of the solution.
    /*!
     * @param gss Vector of length m_kk, which on return
     *           will contain the
     *           standard state gibbs function for species <I>k</I>.
     */
    virtual void getStandardChemPotentials(doublereal* gss) const;

    //!Get the nondimensional gibbs function for the species
    //! standard states at the current T and P of the solution.
    /*!
     * @param grt Vector of length m_kk, which on return
     *           will contain the nondimensional
     *           standard state gibbs function for species <I>k</I>
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    //! Get the array of nondimensional Enthalpy functions for the standard state species
    //! at the current <I>T</I> and <I>P</I> of the solution.
    /*!
     * @param hrt Vector of length m_kk, which on return
     *            will contain the nondimensional
     *            standard state enthalpy of species <I>k</I>
     */
    void getEnthalpy_RT(doublereal* hrt) const;

    //! Get the nondimensional Entropies for the species
    //! standard states at the current T and P of the solution.
    /*!
     * @param sr Vector of length m_kk, which on return
     *           will contain the nondimensional
     *           standard state entropy for species<I>k</I>
     */
    void getEntropy_R(doublereal* sr) const;

    //!   Get the nondimensional heat capacity at constant pressure
    //!   function for the species standard states at the current T and P of the solution.
    /*!
     * @param cpr Vector of length m_kk, which on return
     *           will contain the nondimensional
     *           constant pressure heat capacity for species <I>k</I>
     */
    virtual void getCp_R(doublereal* cpr) const;

    //! Returns the vector of nondimensional
    //!  internal Energies of the standard state at the current
    //! temperature and pressure of the solution for each species.
    /*!
     * @param urt  Output vector of standard state nondimensional internal energies.
     *             Length: m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //@}
    //! @name Thermodynamic Values for the Species Reference State
    /*!
     *  All functions in this group need to be overrided, because
     *  the  m_spthermo SpeciesThermo function is not adequate for
     *  the real equation of state.
     */
    //@{

    //!  Returns the vector of nondimensional
    //!  enthalpies of the reference state at the current temperature
    //!  of the solution and the reference pressure for the species.
    /*!
     * @param hrt     Output vector containing the nondimensional reference state enthalpies
     *                Length: m_kk.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /*!
      *  Returns the vector of nondimensional
      *  enthalpies of the reference state at the current temperature
      *  of the solution and the reference pressure for the species.
      *
      * @param grt     Output vector containing the nondimensional reference state
      *                Gibbs Free energies.  Length: m_kk.
      */
    virtual void getGibbs_RT_ref(doublereal* grt) const;

    /*!
     *  Returns the vector of the gibbs function of the reference state at the
     *  current temperature of the solution and the reference pressure for the
     *  species. units = J/kmol
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
     * @param cprt   Output vector of nondimensional reference state
     *               heat capacities at constant pressure for the species.
     *               Length: m_kk
     */
    virtual void getCp_R_ref(doublereal* cprt) const;

    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const;
    //! @}

    /// critical temperature
    virtual doublereal critTemperature() const;

    /// critical pressure
    virtual doublereal critPressure() const;

    /// critical density
    virtual doublereal critDensity() const;

    /// saturation pressure
    /*!
     * @param t Temperature (kelvin)
     */
    virtual doublereal satPressure(doublereal t);

    //! Return the fraction of vapor at the current conditions
    /*!
     * Below Tcrit, this routine will always return 0, by definition
     * of the functionality of the routine. Above Tcrit, we query
     * the density to toggle between 0 and 1.
     */
    virtual doublereal vaporFraction() const;

    //! Set the temperature of the phase
    /*!
     * The density and composition of the phase is constant during this
     * operator.
     *
     * @param temp Temperature (Kelvin)
     */
    virtual void setTemperature(const doublereal temp);

    //! Set the density of the phase
    /*!
     * The temperature and composition of the phase is constant during this
     * operator.
     *
     * @param dens value of the density in kg m-3
     */
    virtual void setDensity(const doublereal dens);

    //!Import and initialize a ThermoPhase object  using an XML tree.
    /*!
     * @internal
     *
     *   Here we read extra information about the XML description
     *   of a phase. Regular information about elements and species
     *   and their reference state thermodynamic information
     *   have already been read at this point.
     *   For example, we do not need to call this function for
     *   ideal gas equations of state. This function is called from importPhase()
     *   after the elements and the species are initialized with
     *   default ideal solution level data.
     *
     *   The default implementation in ThermoPhase calls the
     *   virtual function initThermo() and then sets the "state" of the
     *   phase by looking for an XML element named "state", and then
     *   interpreting its contents by calling the virtual function
     *   setStateFromXML().
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

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from function importPhase().
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //! Set equation of state parameter values from XML entries.
    /*!
     * This method is called by function importPhase() in
     * file importCTML.cpp when processing a phase definition in
     * an input file. It should be overloaded in subclasses to set
     * any parameters that are specific to that particular phase
     * model. Note, this method is called before the phase is
     * initialized with elements and/or species.
     *
     * @param eosdata An XML_Node object corresponding to
     *                the "thermo" entry for this phase in the input file.
     */
    virtual void setParametersFromXML(const XML_Node& eosdata);

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterPropsIAPWS* getWater() {
        return m_sub;
    }

    //! Get a pointer to a changeable WaterPropsIAPWS object
    WaterProps* getWaterProps() {
        return m_waterProps;
    }

protected:
    /**
     * @internal
     *        This internal routine must be overwritten because
     *  it is not applicable.
     */
    void _updateThermo() const;

private:
    //! Pointer to the WaterPropsIAPWS that calculates the real properties
    //! of water.
    mutable WaterPropsIAPWS* m_sub;

    //! Pointer to the WaterProps object
    /*!
     *   This class is used to house several approximation
     *   routines for properties of water.
     *
     * This object owns m_waterProps, and the WaterPropsIAPWS object used by
     * WaterProps is m_sub, which is defined above.
     */
    WaterProps* m_waterProps;

    //! Molecular weight of Water -> Cantera assumption
    doublereal m_mw;

    //! Offset constants used to obtain consistency with the NIST database.
    /*!
     *  This is added to all internal energy and enthalpy results.
     *  units = J kmol-1.
     */
    doublereal EW_Offset;

    //! Offset constant used to obtain consistency with NIST convention.
    /*!
     *  This is added to all internal entropy results.
     *  units = J kmol-1 K-1.
     */
    doublereal SW_Offset;

    //! Boolean  is true if object has been properly initialized for calculation
    bool m_ready;

    /**
     *  Since this phase represents a liquid phase, it's an error to
     *  return a gas-phase answer. However, if the below is true, then
     *  a gas-phase answer is allowed. This is used to check the thermodynamic
     *  consistency with ideal-gas thermo functions for example.
     */
    bool m_allowGasPhase;
};

}

#endif
