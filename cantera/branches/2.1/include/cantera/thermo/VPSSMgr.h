/**
 *  @file VPSSMgr.h
 * Declaration file for a virtual base class that manages
 * the calculation of standard state properties for all of the
 * species in a single phase, assuming a variable P and T standard state
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr VPSSMgr\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_VPSSMGR_H
#define CT_VPSSMGR_H

#include "cantera/base/ct_defs.h"
#include "mix_defs.h"
#include "cantera/base/global.h"

namespace Cantera
{

class SpeciesThermoInterpType;
class VPStandardStateTP;
class XML_Node;
class SpeciesThermo;
class PDSS;
/**
 * @defgroup mgrpdssthermocalc Managers for Calculating Standard-State Thermodynamics
 *
 * To compute the thermodynamic properties of multicomponent solutions, it is
 * necessary to know something about the thermodynamic properties of the
 * individual species present in the solution. Exactly what sort of species
 * properties are required depends on the thermodynamic model for the
 * solution. For a gaseous solution (i.e., a gas mixture), the species
 * properties required are usually ideal gas properties at the mixture
 * temperature and at a reference pressure (almost always at 1 bar). For other
 * types of solutions, however, it may not be possible to isolate the species
 * in a "pure" state. For example, the thermodynamic properties of, say, Na+
 * and Cl- in saltwater are not easily determined from data on the properties
 * of solid NaCl, or solid Na metal, or chlorine gas. In this case, the
 * solvation in water is fundamental to the identity of the species, and some
 * other reference state must be used. One common convention for liquid
 * solutions is to use thermodynamic data for the solutes in the limit of
 * infinite dilution within the pure solvent; another convention is to
 * reference all properties to unit molality.
 *
 * In defining these standard states for species in a phase, we make the
 * following definition. A reference state is a standard state of a species in
 * a phase limited to one particular pressure, the reference pressure. The
 * reference state specifies the dependence of all thermodynamic functions as
 * a function of the temperature, in between a minimum temperature and a
 * maximum temperature. The reference state also specifies the molar volume of
 * the species as a function of temperature. The molar volume is a
 * thermodynamic function. A full standard state does the same thing as a
 * reference state, but specifies the thermodynamics functions at all
 * pressures.
 *
 * Class VPSSMgr is the base class for a family of classes that compute
 * properties of all species in a phase in their standard states, for a range
 * of temperatures and pressures.
 *
 * Phases which use the VPSSMGr class must have their respective ThermoPhase
 * objects actually be derivatives of the VPStandardState class. These classes
 * assume that there exists a standard state for each species in the phase,
 * where the Thermodynamic functions are specified as a function of
 * temperature and pressure.  Standard state thermo objects for each species
 * in the phase are all derived from the PDSS virtual base class. Calculators
 * for these standard state thermo , which coordinate the calculation for all
 * of the species in a phase, are all derived from VPSSMgr. In turn, these
 * standard states may employ reference state calculation to aid in their
 * calculations. And the VPSSMgr calculators may also employ SimpleThermo
 * calculators to help in calculating the properties for all of the species in
 * a phase. However, there are some PDSS objects which do not employ reference
 * state calculations. An example of this is a real equation of state for
 * liquid water used within the calculation of brine thermodynamics.
 *
 * Typically calls to calculate standard state thermo properties are virtual
 * calls at the ThermoPhase level. It is left to the child classes of
 * ThermoPhase to specify how these are carried out. Usually, this will
 * involve calling the m_spthermo pointer to a SpeciesThermo object to
 * calculate the reference state thermodynamic properties. Then, the pressure
 * dependence is added in within the child ThermoPhase object to complete the
 * specification of the standard state. The VPStandardStateTP class, however,
 * redefines the calls to the calculation of standard state properties to use
 * VPSSMgr class calls.  A listing of these classes and important pointers are
 * supplied below.
 *
 *     - ThermoPhase
 *          - \link Cantera::ThermoPhase::m_spthermo m_spthermo\endlink
 *                 This is a pointer to a %SpeciesThermo manager class that
 *                 handles the reference %state Thermodynamic calculations.
 *     - VPStandardStateTP (inherits from %ThermoPhase)
 *          - \link Cantera::ThermoPhase::m_spthermo m_spthermo\endlink
 *                 %SpeciesThermo manager handling reference %state Thermodynamic calculations.
 *                  may or may not be used by the VPSSMgr class. For species
 *                  which don't have a reference state class defined, a default
 *                  class, called STITbyPDSS which is installed into the SpeciesThermo
 *                  class, actually calculates reference state
 *                  thermo by calling a PDSS object.
 *          - \link Cantera::VPStandardStateTP::m_VPSS_ptr m_VPSS_ptr\endlink
 *                  This is a pointer to a %VPSSMgr class which handles the
 *                  standard %state thermo calculations. It may
 *                  or may not use the pointer, m_spthermo, in its calculations.
 *
 *   The following classes inherit from VPSSMgr. Each of these classes
 *   handle multiple species and by definition all of the species in a phase.
 *   It is a requirement that a VPSSMgr object handles all of the
 *   species in a phase.
 *
 *   - VPSSMgr_IdealGas
 *      - standardState model = "IdealGas"
 *      - This model assumes that all species in the phase obey the
 *        ideal gas law for their pressure dependence. The manager
 *        uses a SpeciesThermo object to handle the calculation of the
 *        reference state.
 *   - VPSSMgr_ConstVol
 *      - standardState model = "ConstVol"
 *      - This model assumes that all species in the phase obey the
 *        constant partial molar volume pressure dependence.
 *        The manager uses a SpeciesThermo object to handle the
 *        calculation of the reference state.
 *   - VPSSMgr_Water_ConstVol
 *      - standardState model = "Water_ConstVol"
 *      - This model assumes that all species but one in the phase obey the
 *        constant partial molar volume pressure dependence.
 *        The manager uses a SpeciesThermo object to handle the
 *        calculation of the reference state for those species.
 *        Species 0 is assumed to be water, and a real equation
 *        of state is used to model the T, P behavior.
 *   - VPSSMgr_Water_HKFT
 *      - standardState model = "Water_HKFT"
 *      - This model assumes that all species but one in the phase obey the
 *        HKFT equation of state.
 *        Species 0 is assumed to be water, and a real equation
 *        of state is used to model the T, P behavior.
 *   - VPSSMgr_General
 *      - standardState model = "General"
 *      - This model is completely general. Nothing is assumed at this
 *        level. Calls consist of loops to PDSS property evaluations.
 *
 *  The choice of which VPSSMgr object to be used is implicitly made by
 *  %Cantera by querying the XML data file for compatibility.
 *  However, each of these VPSSMgr objects may be explicitly requested in the XML file
 *  by adding in the following XML node into the thermo section of the
 *  phase XML Node. For example, the code example listed below
 *  explicitly requests that the VPSSMgr_IdealGas
 *  object be used to handle the standard state thermodynamics calculations.
 *
 *  @code
 *    <phase id="Silane_Pyrolysis" dim="3">
 *       . . .
 *       <thermo model="VPIdealGas">
 *          <standardState model="IdealGas"\>
 *       <\thermo>
 *       . . .
 *    <\phase>
 *  @endcode
 *
 *  If it turns out that the VPSSMgr_IdealGas class can not handle the standard
 *  state calculation, then %Cantera will fail during the instantiation phase
 *  printing out an informative error message.
 *
 *  In the source code listing above, the thermo model, VPIdealGas ,was requested. The
 *  thermo model specifies the type of ThermoPhase object to use. In this case
 *  the object IdealSolnGasVPSS (with the ideal gas suboption) is used. IdealSolnGasVPSS
 *  inherits from VPStandardStateTP, so that it actually has a VPSSMgr pointer
 *  to be specified. Note, in addition to the IdealGas entry to the model
 *  parameter in standardState node, we could have also specified the "General"
 *  option. The general option will always work. An example of this
 *  usage is listed below.
 *
 *  @code
 *    <phase id="Silane_Pyrolysis" dim="3">
 *       . . .
 *       <thermo model="VPIdealGas">
 *          <standardState model="General"\>
 *       <\thermo>
 *       . . .
 *    <\phase>
 *  @endcode
 *
 *  The "General" option will cause the VPSSMgr_General %VPSSMgr class to be
 *  used. In this manager, the calculations are all handled at the PDSS object
 *  level. This is completely general, but, may be significantly slower.
 *
 * @ingroup thermoprops
 */

//! Virtual base class for the classes that manage the calculation
//! of standard state properties for all the species in a phase.
/*!
 *  This class defines the interface which all subclasses must implement.
 *
 * Class VPSSMgr is the base class for a family of classes that compute
 * properties of a set of species in their standard state at a range of
 * temperatures and pressures.
 *
 *  If #m_useTmpRefStateStorage is set to true, then the following internal
 *  arrays, containing information about the reference arrays,
 *  are calculated and kept up to date at every call.
 *
 *  - #m_h0_RT
 *  - #m_g0_RT
 *  - #m_s0_R
 *  - #m_cp0_R
 *
 *  The virtual function #_updateRefStateThermo() is supplied to do this
 *  and may be reimplemented in child routines. A default implementation
 *  based on the speciesThermo class is supplied in this base class.
 *  #_updateStandardStateThermo() is called whenever a reference state
 *   property is needed.
 *
 *  When  #m_useTmpStandardStateStorage is true, then the following
 *  internal arrays, containing information on the standard state properties
 *  are calculated and kept up to date.
 *
 *  -  #m_hss_RT;
 *  -  #m_cpss_R;
 *  -  #m_gss_RT;
 *  -  #m_sss_R;
 *  -  #m_Vss
 *
 *  The virtual function #_updateStandardStateThermo() is supplied to do this
 *  and must be reimplemented in child routines,
 *  when  #m_useTmpStandardStateStorage is true.
 *  It may be optionally reimplemented in child routines if
 *  #m_useTmpStandardStateStorage is false.
 *  #_updateStandardStateThermo() is called whenever a standard state property is needed.
 *
 *  This class is usually used for nearly incompressible phases. For those phases, it
 *  makes sense to change the equation of state independent variable from
 *  density to pressure.
 *
 */
class VPSSMgr
{
public:
    //! Constructor
    /*!
     * @param vptp_ptr Pointer to the Variable pressure %ThermoPhase object
     *                 This object must have already been malloced.
     * @param spth     Pointer to the optional SpeciesThermo object
     *                 that will handle the calculation of the reference
     *                 state thermodynamic coefficients.
     */
    VPSSMgr(VPStandardStateTP* vptp_ptr, SpeciesThermo* spth = 0);

    //! Destructor
    virtual ~VPSSMgr();

    //! Copy Constructor
    VPSSMgr(const VPSSMgr& right);

    //! Assignment operator
    VPSSMgr& operator=(const VPSSMgr& right);

    //! Duplication routine for objects which derive from VPSSMgr
    /*!
     *  This function can be used to duplicate objects derived from VPSSMgr
     *  even if the application only has a pointer to VPSSMgr to work with.
     */
    virtual VPSSMgr* duplMyselfAsVPSSMgr() const;

    //! @name  Properties of the Standard State of the Species in the Solution
    //! @{

    //!Get the array of chemical potentials at unit activity.
    /*!
     * These are the standard state chemical potentials \f$ \mu^0_k(T,P)
     * \f$. The values are evaluated at the current temperature and pressure.
     *
     * @param mu   Output vector of standard state chemical potentials.
     *             length = m_kk. units are J / kmol.
     */
    virtual void getStandardChemPotentials(doublereal* mu) const;

    /**
     * Get the nondimensional Gibbs functions for the species at their
     * standard states of solution at the current T and P of the solution.
     *
     * @param grt    Output vector of nondimensional standard state
     *               Gibbs free energies. length = m_kk.
     */
    virtual void getGibbs_RT(doublereal* grt) const;

    /**
     * Get the nondimensional Enthalpy functions for the species at their
     * standard states at the current *T* and *P* of the solution.
     *
     * @param hrt     Output vector of standard state enthalpies.
     *                length = m_kk. units are unitless.
     */
    virtual void getEnthalpy_RT(doublereal* hrt) const;

    //! Return a reference to a vector of the molar enthalpies of the
    //! species in their standard states
    const vector_fp& enthalpy_RT() const {
        return m_hss_RT;
    }

    /**
     * Get the array of nondimensional Enthalpy functions for the standard
     * state species at the current *T* and *P* of the solution.
     *
     * @param sr     Output vector of nondimensional standard state
     *               entropies. length = m_kk.
     */
    virtual void getEntropy_R(doublereal* sr) const;

    //! Return a reference to a vector of the entropies of the species
    const vector_fp& entropy_R() const {
        return m_sss_R;
    }

    //! Returns the vector of nondimensional internal Energies of the standard
    //! state at the current temperature and pressure of the solution for each
    //! species.
    /*!
     * The internal energy is calculated from the enthalpy from the
     * following formula:
     *
     * \f[
     *  u^{ss}_k(T,P) = h^{ss}_k(T)  - P * V^{ss}_k
     * \f]
     *
     * @param urt    Output vector of nondimensional standard state
     *               internal energies. length = m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    //! Get the nondimensional Heat Capacities at constant pressure for the
    //! standard state of the species at the current T and P.
    /*!
     * This is redefined here to call the internal function,  _updateStandardStateThermo(),
     * which calculates all standard state properties at the same time.
     *
     * @param cpr    Output vector containing the the nondimensional Heat
     *               Capacities at constant pressure for the standard state of
     *               the species. Length: m_kk.
     */
    virtual void getCp_R(doublereal* cpr) const;

    //! Return a reference to a vector of the constant pressure
    //! heat capacities of the species
    const vector_fp& cp_R() const {
        return m_cpss_R;
    }

    //! Get the molar volumes of each species in their standard states at the
    //! current *T* and *P* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * This is redefined here to call the internal function,
     *  _updateStandardStateThermo(), which calculates all standard state
     *  properties at the same time.
     *
     * @param vol Output vector of species volumes. length = m_kk.
     *            units =  m^3 / kmol
     */
    virtual void getStandardVolumes(doublereal* vol) const;
    virtual const vector_fp& getStandardVolumes() const;

    //! Return a reference to a vector of the species standard molar volumes
    //! @deprecated Use getStandardVolumes()
    const vector_fp& standardVolumes() const {
        warn_deprecated("VPSSMgr::standardVolumes");
        return m_Vss;
    }

public:
    //@}
    /*! @name Thermodynamic Values for the Species Reference States
     *  There are also temporary variables for holding the species reference-
     *  state values of Cp, H, S, and V at the last temperature and reference
     *  pressure called. These functions are not recalculated if a new call is
     *  made using the previous temperature. All calculations are done within
     *  the routine _updateRefStateThermo().
     */
    //@{

    /*!
     *  Returns the vector of nondimensional enthalpies of the reference state
     *  at the current temperature of the solution and the reference pressure
     *  for the species.
     *
     * @param hrt Output vector contains the nondimensional enthalpies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /*!
     *  Returns the vector of nondimensional Gibbs free energies of the
     *  reference state at the current temperature of the solution and the
     *  reference pressure for the species.
     *
     * @param grt Output vector contains the nondimensional Gibbs free energies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const ;


    //! Return a reference to the vector of Gibbs free energies of the species
    const vector_fp& Gibbs_RT_ref() const {
        return m_g0_RT;
    }

    /*!
     *  Returns the vector of the gibbs function of the reference state at the
     *  current temperature of the solution and the reference pressure for the
     *  species. units = J/kmol
     *
     * @param g   Output vector contain the Gibbs free energies
     *            of the reference state of the species
     *            length = m_kk, units = J/kmol.
     */
    virtual void getGibbs_ref(doublereal* g) const ;

    /*!
     *  Returns the vector of nondimensional entropies of the reference state
     *  at the current temperature of the solution and the reference pressure
     *  for the species.
     *
     * @param er  Output vector contain the nondimensional entropies
     *            of the species in their reference states
     *            length: m_kk, units: dimensionless.
     */
    virtual void getEntropy_R_ref(doublereal* er) const ;

    /*!
     *  Returns the vector of nondimensional constant pressure heat capacities
     *  of the reference state at the current temperature of the solution and
     *  reference pressure for the species.
     *
     * @param cpr  Output vector contains the nondimensional heat capacities
     *             of the species in their reference states
     *             length: m_kk, units: dimensionless.
     */
    virtual void getCp_R_ref(doublereal* cpr) const ;

    //!  Get the molar volumes of the species reference states at the current
    //!  *T* and *P_ref* of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const ;

    //@}
    /*! @name Setting the Internal State of the System
     *  All calls to change the internal state of the system's T and P
     *  are done through these routines
     *      - setState_TP()
     *      - setState_T()
     *      - setState_P()
     *
     *  These routine in turn call the following underlying virtual functions
     *
     *   - _updateRefStateThermo()
     *   - _updateStandardStateThermo()
     *
     *  An important point to note is that between calls the assumption
     *  that the underlying PDSS objects will retain their set Temperatures
     *  and Pressure CAN NOT BE MADE. For efficiency reasons, we may twiddle
     *  these to get derivatives.
     */
    //@{

    //! Set the temperature (K) and pressure (Pa)
    /*!
     *  This sets the temperature and pressure and triggers
     *  calculation of underlying quantities
     *
     * @param T    Temperature (K)
     * @param P    Pressure (Pa)
     */
    virtual void setState_TP(doublereal T, doublereal P);

    //! Set the temperature (K)
    /*!
     * @param T    Temperature (K)
     */
    virtual void setState_T(doublereal T);

    //! Set the  pressure (Pa)
    /*!
     * @param P    Pressure (Pa)
     */
    virtual void setState_P(doublereal P);

    //! Return the temperature stored in the object
    doublereal temperature() const {
        return m_tlast;
    }

    //! Return the pressure stored in the object
    doublereal pressure() const {
        return m_plast;
    }

    //! Return the pointer to the reference-state Thermo calculator
    //! SpeciesThermo object.
    SpeciesThermo* SpeciesThermoMgr() {
        return m_spthermo;
    }

    //! Updates the internal standard state thermodynamic vectors at the
    //! current T and P of the solution.
    /*!
     * If you are to peek internally inside the object, you need to
     * call these functions after setState functions in order to be sure
     * that the vectors are current.
     */
    virtual void updateStandardStateThermo();

    //! Updates the internal reference state thermodynamic vectors at the
    //! current T of the solution and the reference pressure.
    /*!
     * If you are to peek internally inside the object, you need to
     * call these functions after setState functions in order to be sure
     * that the vectors are current.
     */
    virtual void updateRefStateThermo() const;

protected:

    //! Updates the standard state thermodynamic functions at the
    //! current T and P of the solution.
    /*!
     * @internal
     *
     * If m_useTmpStandardStateStorage is true, this function must be called
     * for every call to functions in this class. It checks to see whether the
     * temperature or pressure has changed and thus the ss thermodynamics
     * functions for all of the species must be recalculated.
     *
     * This function is responsible for updating the following internal members,
     * when  m_useTmpStandardStateStorage is true.
     *
     *  -  m_hss_RT;
     *  -  m_cpss_R;
     *  -  m_gss_RT;
     *  -  m_sss_R;
     *  -  m_Vss
     *
     *  If m_useTmpStandardStateStorage is not true, this function may be
     *  required to be called by child classes to update internal member data.
     *
     *  Note, the base class implementation will throw an error. It must be
     *  reimplemented in derived classes.
     *
     *  Underscore updates never check for the state of the system
     *  They just do the calculation.
     */
    virtual void _updateStandardStateThermo();

    //! Updates the reference state thermodynamic functions at the
    //! current T of the solution and the reference pressure
    /*!
     *  Underscore updates never check for the state of the system
     *  They just do the calculation.
     */
    virtual void _updateRefStateThermo() const;

public:
    //@}
    //! @name Utility Methods - Reports on various quantities
    /*!
     * The following methods are used in the process of reporting
     * various states and attributes
     */
    //@{

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     * @param index  Species index
     */
    virtual PDSS_enumType reportPDSSType(int index = -1) const ;

    //! This utility function reports the type of manager
    //! for the calculation of ss properties
    /*!
     *  @return Returns an enum type called VPSSMgr_enumType, which is a list
     *          of the known VPSSMgr objects
     */
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;

    //! Minimum temperature.
    /*!
     * If no argument is supplied, this method returns the minimum temperature
     * for which \e all parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the minimum temperature for
     * species k in the phase.
     *
     * @param k    Species index
     */
    virtual doublereal minTemp(size_t k=npos) const ;

    //! Maximum temperature.
    /*!
     * If no argument is supplied, this method returns the maximum temperature
     * for which \e all parameterizations are valid. If an integer index k is
     * supplied, then the value returned is the maximum temperature for
     * parameterization k.
     *
     * @param k  Species Index
     */
    virtual doublereal maxTemp(size_t k=npos) const;

    //! The reference-state pressure for the standard state
    /*!
     * Returns the reference state pressure in Pascals for species k. If k is
     * left out of the argument list, it returns the reference state pressure
     * for the first species. Note that some SpeciesThermo implementations,
     * such as those for ideal gases, require that all species in the same
     * phase have the same reference state pressures.
     *
     * @param k Species index. Default is -1, which returns
     *          the generic answer.
     */
    virtual doublereal refPressure(size_t k=npos) const ;

    //@}
    /*! @name Initialization Methods - For Internal use
     * The following methods are used in the process of constructing the phase
     * and setting its parameters from a specification in an input file. They
     * are not normally used in application programs. To see how they are
     * used, see files importCTML.cpp and ThermoFactory.cpp.
     */
    //@{

    //! @internal Initialize the object
    /*!
     * This method is provided to allow subclasses to perform any
     * initialization required after all species have been added. For example,
     * it might be used to resize internal work arrays that must have an entry
     * for each species.  The base class implementation does nothing, and
     * subclasses that do not require initialization do not need to overload
     * this method.  When importing a CTML phase description, this method is
     * called just prior to returning from function importPhase().
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //! Initialize the lengths within the object
    /*!
     *  Note this function is not virtual
     */
    void initLengths();

    //! Finalize the thermo after all species have been entered
    /*!
     *  This function is the LAST initialization routine to be called. It's
     *  called after createInstallPDSS() has been called for each species in
     *  the phase, and after initThermo() has been called. It's called via an
     *  inner-to-outer onion shell like manner.
     *
     *  In this routine, we currently calculate the reference pressure,
     *  the minimum and maximum temperature for the applicability
     *  of the thermo formulation.
     *
     *  @param phaseNode   Reference to the phaseNode XML node.
     *  @param id          ID of the phase.
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Install specific content for species k in the reference-state
    //! thermodynamic SpeciesManager object
    /*!
     * This occurs before matrices are sized appropriately.
     *
     * @param k           Species index in the phase
     * @param speciesNode XML Node corresponding to the species
     * @param phaseNode_ptr Pointer to the XML Node corresponding
     *                      to the phase which owns the species
     */
    void installSTSpecies(size_t k,  const XML_Node& speciesNode,
                          const XML_Node* phaseNode_ptr);

    //! Install specific content for species k in the standard-state
    //! thermodynamic calculator and also create/return a PDSS object
    //! for that species.
    /*!
     * This occurs before matrices are sized appropriately.
     *
     * @param k           Species index in the phase
     * @param speciesNode XML Node corresponding to the species
     * @param phaseNode_ptr Pointer to the XML Node corresponding
     *                      to the phase which owns the species
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);

    //! Initialize the internal shallow pointers in this object
    /*!
     * There are a bunch of internal shallow pointers that point to the owning
     * VPStandardStateTP and SpeciesThermo objects. This function reinitializes
     * them. This function is called like an onion.
     *
     *  @param vp_ptr   Pointer to the VPStandardStateTP standard state
     *  @param sp_ptr   Pointer to the SpeciesThermo standard state
     */
    virtual void initAllPtrs(VPStandardStateTP* vp_ptr, SpeciesThermo* sp_ptr);

protected:
    //! Number of species in the phase
    size_t m_kk;

    //! Variable pressure ThermoPhase object
    VPStandardStateTP* m_vptp_ptr;

    //!  Pointer to reference state thermo calculator
    /*!
     * Note, this can have a value of 0
     */
    SpeciesThermo* m_spthermo;

    //! The last temperature at which the standard state thermodynamic
    //! properties were calculated at.
    mutable doublereal    m_tlast;

    //! The last pressure at which the Standard State thermodynamic
    //! properties were calculated at.
    mutable doublereal    m_plast;

    /*!
     * Reference pressure (Pa) must be the same for all species
     * - defaults to 1 atm.
     */
    mutable doublereal m_p0;

    //! minimum temperature for the standard state calculations
    doublereal m_minTemp;

    //! maximum temperature for the standard state calculations
    doublereal m_maxTemp;

    /*!
     * boolean indicating whether temporary reference state storage is used
     * -> default is false
     */
    bool m_useTmpRefStateStorage;

    /*!
     * Vector containing the species reference enthalpies at T = m_tlast
     * and P = p_ref.
     */
    mutable vector_fp      m_h0_RT;

    /**
     * Vector containing the species reference constant pressure
     * heat capacities at T = m_tlast    and P = p_ref.
     */
    mutable vector_fp      m_cp0_R;

    /**
     * Vector containing the species reference Gibbs functions
     * at T = m_tlast  and P = p_ref.
     */
    mutable vector_fp      m_g0_RT;

    /**
     * Vector containing the species reference entropies
     * at T = m_tlast and P = p_ref.
     */
    mutable vector_fp      m_s0_R;

    //! Vector containing the species reference molar volumes
    mutable vector_fp      m_V0;

    /*!
     * boolean indicating whether temporary standard state storage is used
     * -> default is false
     */
    bool m_useTmpStandardStateStorage;

    /**
     * Vector containing the species Standard State enthalpies at T = m_tlast
     * and P = m_plast.
     */
    mutable vector_fp      m_hss_RT;

    /**
     * Vector containing the species Standard State constant pressure
     * heat capacities at T = m_tlast and P = m_plast.
     */
    mutable vector_fp      m_cpss_R;

    /**
     * Vector containing the species Standard State Gibbs functions
     * at T = m_tlast and P = m_plast.
     */
    mutable vector_fp      m_gss_RT;

    /**
     * Vector containing the species Standard State entropies
     * at T = m_tlast and P = m_plast.
     */
    mutable vector_fp      m_sss_R;

    /**
     * Vector containing the species standard state volumes
     * at T = m_tlast and P = m_plast
     */
    mutable vector_fp      m_Vss;

    //! species reference enthalpies - used by individual PDSS objects
    /*!
     * Vector containing the species reference enthalpies at T = m_tlast
     * and P = p_ref.
     */
    mutable vector_fp      mPDSS_h0_RT;

    //! species reference heat capacities - used by individual PDSS objects
    /**
     * Vector containing the species reference constant pressure
     * heat capacities at T = m_tlast    and P = p_ref.
     */
    mutable vector_fp      mPDSS_cp0_R;

    //! species reference gibbs free energies - used by individual PDSS objects
    /**
     * Vector containing the species reference Gibbs functions
     * at T = m_tlast  and P = p_ref.
     */
    mutable vector_fp      mPDSS_g0_RT;

    //! species reference entropies - used by individual PDSS objects
    /**
     * Vector containing the species reference entropies
     * at T = m_tlast and P = p_ref.
     */
    mutable vector_fp      mPDSS_s0_R;

    //! species reference state molar Volumes - used by individual PDSS objects
    /**
     * Vector containing the rf molar volumes
     * at T = m_tlast and P = p_ref.
     */
    mutable vector_fp      mPDSS_V0;

    //! species standard state enthalpies - used by individual PDSS objects
    /*!
     * Vector containing the species standard state enthalpies at T = m_tlast
     * and P = p_ref.
     */
    mutable vector_fp      mPDSS_hss_RT;

    //! species standard state heat capacities - used by individual PDSS objects
    /**
     * Vector containing the species standard state constant pressure
     * heat capacities at T = m_tlast    and P = p_ref.
     */
    mutable vector_fp      mPDSS_cpss_R;

    //! species standard state gibbs free energies - used by individual PDSS objects
    /**
     * Vector containing the species standard state Gibbs functions
     * at T = m_tlast  and P = p_ref.
     */
    mutable vector_fp      mPDSS_gss_RT;

    //! species standard state entropies - used by individual PDSS objects
    /**
     * Vector containing the species standard state entropies
     * at T = m_tlast and P = p_ref.
     */
    mutable vector_fp      mPDSS_sss_R;

    //! species standard state molar Volumes - used by individual PDSS objects
    /**
     * Vector containing the ss molar volumes
     * at T = m_tlast and P = p_ref.
     */
    mutable vector_fp      mPDSS_Vss;

    friend class PDSS;
private:
    //! Error message to indicate an unimplemented feature
    /*!
     * @param msg  Error message string
     */
    void err(const std::string& msg) const;
};
//@}
}

#endif
