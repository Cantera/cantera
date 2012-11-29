/**
 *  @file VPSSMgr_Water_HKFT.h
 * Declaration file for a derived class that handles the calculation
 * of standard state thermo properties for real water and
 *  a set of species which have the HKFT equation of state
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr_Water_HKFT VPSSMgr_Water_HKFT\endlink).
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_VPSSMGR_WATER_HKFT_H
#define CT_VPSSMGR_WATER_HKFT_H

#include "VPSSMgr.h"

namespace Cantera
{

class SpeciesThermoInterpType;
class VPStandardStateTP;
class SpeciesThermo;
class PDSS;
class PDSS_Water;

//! Virtual base class for the species thermo manager classes.
/*!
 *  This class defines the interface which all subclasses must implement.
 *
 * Class %VPSSSpeciesThermo is the base class
 * for a family of classes that compute properties of a set of
 * species in their reference state at a range of temperatures.
 * Note, the pressure dependence of the reference state is not
 * handled by this particular species standard state model.
 *
 * @ingroup mgrpdssthermocalc
 */
class VPSSMgr_Water_HKFT : public VPSSMgr
{

public:


    //! Constructor
    /*!
     * @param vptp_ptr Pointer to the Variable pressure %ThermoPhase object
     *                 This object must have already been malloced.
     *
     * @param spth     Pointer to the optional SpeciesThermo object
     *                 that will handle the calculation of the reference
     *                 state thermodynamic coefficients.
     */
    VPSSMgr_Water_HKFT(VPStandardStateTP* vptp_ptr,
                       SpeciesThermo* spth);

    //! Destructor
    virtual ~VPSSMgr_Water_HKFT();

    //! Copy Constructor for the %SpeciesThermo object.
    /*!
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    VPSSMgr_Water_HKFT(const VPSSMgr_Water_HKFT& right);

    //! Assignment operator for the %SpeciesThermo object
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    VPSSMgr_Water_HKFT& operator=(const VPSSMgr_Water_HKFT& right);

    //! Duplication routine for objects which inherit from
    //! %VPSSSpeciesThermo
    /*!
     *  This virtual routine can be used to duplicate %VPSSSpeciesThermo  objects
     *  inherited from %VPSSSpeciesThermo even if the application only has
     *  a pointer to %VPSSSpeciesThermo to work with.
     */
    virtual VPSSMgr* duplMyselfAsVPSSMgr() const;

    /*!
     * @name  Properties of the Standard State of the Species in the Solution
     *
     *  Within VPStandardStateTP, these properties are calculated via a common routine,
     *  _updateStandardStateThermo(),
     *  which must be overloaded in inherited objects.
     *  The values are cached within this object, and are not recalculated unless
     *  the temperature or pressure changes.
     */
    //@{


    //@}
    /// @name Thermodynamic Values for the Species Reference States (VPStandardStateTP)
    /*!
     *  There are also temporary
     *  variables for holding the species reference-state values of Cp, H, S, and V at the
     *  last temperature and reference pressure called. These functions are not recalculated
     *  if a new call is made using the previous temperature.
     *  All calculations are done within the routine  _updateRefStateThermo().
     */
    //@{

    /*!
     *  Returns the vector of nondimensional
     *  enthalpies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param hrt Output vector contains the nondimensional enthalpies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;

    /*!
     *  Returns the vector of nondimensional
     *  Gibbs free energies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param grt Output vector contains the nondimensional Gibbs free energies
     *            of the reference state of the species
     *            length = m_kk, units = dimensionless.
     */
    virtual void getGibbs_RT_ref(doublereal* grt) const ;

    /*!
      *  Returns the vector of the
      *  gibbs function of the reference state at the current temperature
      *  of the solution and the reference pressure for the species.
      *  units = J/kmol
      *
      * @param g   Output vector contain the Gibbs free energies
      *            of the reference state of the species
      *            length = m_kk, units = J/kmol.
      */
    virtual void getGibbs_ref(doublereal* g) const ;


    /*!
     *  Returns the vector of nondimensional
     *  entropies of the reference state at the current temperature
     *  of the solution and the reference pressure for the species.
     *
     * @param er  Output vector contain the nondimensional entropies
     *            of the species in their reference states
     *            length: m_kk, units: dimensionless.
     */
    virtual void getEntropy_R_ref(doublereal* er) const ;

    /*!
     *  Returns the vector of nondimensional
     *  constant pressure heat capacities of the reference state
     *  at the current temperature of the solution
     *  and reference pressure for the species.
     *
     * @param cpr  Output vector contains the nondimensional heat capacities
     *             of the species in their reference states
     *             length: m_kk, units: dimensionless.
     */
    virtual void getCp_R_ref(doublereal* cpr) const ;

    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const ;

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

    //@}
    /// @name Setting the Internal State of the System
    /*!
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

    //! Updates the internal reference state thermodynamic vectors at the
    //! current T of the solution and the reference pressure.
    /*!
     *  This is called to make sure that the internal thermodynamic members
     *  are up-to-date. It checks against an internal value of m_tempRef
     *  to see whether the values are current.
     */
    virtual void updateRefStateThermo() const;

private:

    //! Updates the reference state thermodynamic functions at the current T
    //! and a calculated Pref that is safe.
    /*!
     *
     * This function is responsible for updating the following internal members
     *
     *  -  m_h0_RT;
     *  -  m_cp0_R;
     *  -  m_g0_RT;
     *  -  m_s0_R;
     *  -  m_V0
     *
     *  It always does the calculation. No checking is ever done to see
     *  if the calculation is necessary.
     *
     *  m_p0 is calculated within this routine given the value of the temperature.
     *  This is necessary because we are using a real equation of state for
     *  water.
     *
     *  The state of the system is left at (m_tlast, m_plast) at the end
     *  of the routine.
     */
    virtual void _updateRefStateThermo() const;

    //! Updates the standard state thermodynamic functions at the current T and P of the solution.
    /*!
     * @internal
     *
     * If m_useTmpStandardStateStorage is true,
     * this function must be called for every call to functions in this
     * class. It checks to see whether the temperature or pressure has changed and
     * thus the ss thermodynamics functions for all of the species
     * must be recalculated.
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
     *  Note, this will throw an error. It must be reimplemented in derived classes.
     *
     */
    virtual void _updateStandardStateThermo();


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
     *
     * @param index  Species index
     */
    virtual PDSS_enumType reportPDSSType(int index = -1) const ;


    //! This utility function reports the type of manager
    //! for the calculation of ss properties
    /*!
     *
     *
     */
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;

    //@}
    //! @name Initialization Methods - For Internal use (VPStandardState)
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see files importCTML.cpp and
     * ThermoFactory.cpp.
     */
    //@{

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
     *
     * @see importCTML.cpp
     */
    virtual void initThermo();

    //! Finalize the thermo after all species have been entered
    /*!
     *  This function is the LAST initialization routine to be
     *  called. It's called after createInstallPDSS() has been
     *  called for each species in the phase, and after initThermo()
     *  has been called.
     *  It's called via an inner-to-outer onion shell like manner.
     *
     *
     *  @param phaseNode   Reference to the phaseNode XML node.
     *  @param id          ID of the phase.
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //! Install specific content for species k in the standard-state
    //! thermodynamic calculator and also create/return a PDSS object
    //! for that species.
    /*!
     * This occurs before matrices are sized appropriately.
     *
     * @param k             Species index in the phase
     * @param speciesNode   XML Node corresponding to the species
     * @param phaseNode_ptr Pointer to the XML Node corresponding
     *                      to the phase which owns the species
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);

    //@}

private:

    //! Shallow pointer to the water object
    PDSS_Water* m_waterSS;

    //! Last reference temperature calculated
    /*!
     * Reference state calculations are totally separated from
     * standard state calculations.
     */
    mutable doublereal m_tlastRef;
};
//@}
}

#endif

