/**
 *  @file VPSSMgr_ConstVol.h
 *  Declarations for a derived class for the calculation of multiple-species thermodynamic
 *  property managers for variable temperature and pressure standard
 *  states assuming constant volume (see class
 *  \link Cantera::VPSSMgr_ConstVol VPSSMgr_ConstVol \endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_VPSSMGR_CONSTVOL_H
#define CT_VPSSMGR_CONSTVOL_H

#include "cantera/base/ct_defs.h"
#include "VPSSMgr.h"

namespace Cantera
{

class SpeciesThermoInterpType;
class PDSS;

//! Constant Molar Volume e VPSS species thermo manager class
/*!
 *  The calculation of multiple-species thermodynamic
 *  property managers for variable temperature and pressure standard
 *  states assuming a constant partial molar volume assumption.
 *
 *  @ingroup mgrpdssthermocalc
 */
class VPSSMgr_ConstVol : public VPSSMgr
{

public:

    //! Constructor
    /*!
     *  @param vp_ptr Pointer to the owning VPStandardStateTP  object
     *                for the phase. It's a requirement that this be
     *                already malloced.
     *  @param spth   Pointer to the SpeciesThermo object for the
     *                phase. It's a requirement that this be already
     *                malloced.
     */
    VPSSMgr_ConstVol(VPStandardStateTP* vp_ptr, SpeciesThermo* spth);

    //! Destructor
    virtual ~VPSSMgr_ConstVol();

    //! Copy Constructor
    /*!
     * @param right    Reference to %VPSSMgr_ConstVol object to be copied into the
     *                 current one.
     */
    VPSSMgr_ConstVol(const VPSSMgr_ConstVol& right);

    //! Assignment operator for the %VPSSMgr_ConstVol object
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %VPSSMgr_ConstVol object to be copied into the
     *                 current one.
     */
    VPSSMgr_ConstVol& operator=(const VPSSMgr_ConstVol& right);

    //! Duplicator routine for the VPSSMgr base class
    /*!
     *  This virtual routine can be used to duplicate %VPSSMgr objects
     *  inherited from %VPSSMgr even if the application only has
     *  a pointer to %VPSSMgr to work with.
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

protected:

    //! Updates the standard state thermodynamic functions at the current
    //! T and P of the solution.
    /*!
     * @internal
     *
     * If m_useTmpStandardStateStorage is true,
     * this function must be called whenever the temperature or pressure
     * has changed.
     *
     * This function is responsible for updating the following internal members,
     * when m_useTmpStandardStateStorage is true.
     *
     *  -  m_hss_RT;
     *  -  m_cpss_R;
     *  -  m_gss_RT;
     *  -  m_sss_R;
     *  -  m_Vss
     *
     *  If m_useTmpStandardStateStorage is not true, this function may be
     *  required to be called every time information is requested from
     *  this object.
     */
    virtual void _updateStandardStateThermo();

    //@}

    /// @name Thermodynamic Values for the Species Reference States
    /*!
     *  There are also temporary
     *  variables for holding the species reference-state values of Cp, H, S, and V at the
     *  last temperature and reference pressure called. These functions are not recalculated
     *  if a new call is made using the previous temperature.
     *  All calculations are done within the routine  _updateRefStateThermo().
     *  _updateRefStateThermo() is defined in the parent object.
     */
    //@{

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


    //!  Get the molar volumes of the species reference states at the current
    //!  <I>T</I> and <I>P_ref</I> of the solution.
    /*!
     * units = m^3 / kmol
     *
     * @param vol     Output vector containing the standard state volumes.
     *                Length: m_kk.
     */
    virtual void getStandardVolumes_ref(doublereal* vol) const ;

    //@}

    //! @name Initialization Methods - For Internal use
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally seen by application programs
     */
    //@{

public:
    //! Initialize the VPSSMgr object
    /*!
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species. It is called after createInstallPDSS() and
     * before initThermoXML().
     *
     * @internal
     */
    virtual void initThermo();

    //! Initialize the thermo for this standard state thermo calculator
    /*!
     *  This task is done last, after createInstallPDSS() and after
     *  initThermo().
     *
     *  @param phaseNode   Reference to the phase node in the XML tree
     *  @param id          string name of the phase
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //!  Create and install a constant volume pressure dependent
    //!  standard state for one species within this object
    /*!
     *  This function sets up the internal data within this object for
     *  handling the calculation of the standard state for the species.
     *
     *  -   It registers the species with the SpeciesThermo object for the
     *      containing VPStandardStateTP phase.
     *  -   It grabs the molar volume property and installs its value within
     *      this object.
     *  -   It also creates a PDSS object, which basically contains a
     *      duplication of some of this information and returns a pointer to
     *      the new object.
     *  .
     *
     *  @param k Species index within the phase
     *  @param speciesNode Reference to the species node in the XML tree
     *  @param phaseNode_ptr Pointer to the phase node in the XML tree
     *
     *  @return Returns a pointer to the a newly malloced PDSS object
     *          containing the parameterization
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);
    //@}

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
     */
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;

};
//@}
}

#endif
