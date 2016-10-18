/**
 *  @file VPSSMgr_ConstVol.h
 *  Declarations for a derived class for the calculation of multiple-species thermodynamic
 *  property managers for variable temperature and pressure standard
 *  states assuming constant volume (see class
 *  \link Cantera::VPSSMgr_ConstVol VPSSMgr_ConstVol \endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_VPSSMGR_CONSTVOL_H
#define CT_VPSSMGR_CONSTVOL_H

#include "VPSSMgr.h"

namespace Cantera
{
//! Constant Molar Volume e VPSS species thermo manager class
/*!
 * The calculation of multiple-species thermodynamic property managers for
 * variable temperature and pressure standard states assuming a constant partial
 * molar volume assumption.
 *
 *  @ingroup mgrpdssthermocalc
 */
class VPSSMgr_ConstVol : public VPSSMgr
{
public:
    //! Constructor
    /*!
     *  @param vp_ptr Pointer to the owning VPStandardStateTP object for the
     *                phase.
     *  @param spth   Pointer to the MultiSpeciesThermo object for the phase.
     */
    VPSSMgr_ConstVol(VPStandardStateTP* vp_ptr, MultiSpeciesThermo* spth);

    VPSSMgr_ConstVol(const VPSSMgr_ConstVol& right);
    VPSSMgr_ConstVol& operator=(const VPSSMgr_ConstVol& right);
    virtual VPSSMgr* duplMyselfAsVPSSMgr() const;

    /*!
     * @name  Properties of the Standard State of the Species in the Solution
     *
     * Within VPStandardStateTP, these properties are calculated via a common
     * routine, _updateStandardStateThermo(), which must be overloaded in
     * inherited objects. The values are cached within this object, and are
     * not recalculated unless the temperature or pressure changes.
     */
    //@{

protected:
    virtual void _updateStandardStateThermo();

    //@}
    /*! @name Thermodynamic Values for the Species Reference States
     *
     *  There are also temporary variables for holding the species reference-
     *  state values of Cp, H, S, and V at the last temperature and reference
     *  pressure called. These functions are not recalculated if a new call is
     *  made using the previous temperature. All calculations are done within
     *  the routine _updateRefStateThermo(). _updateRefStateThermo() is
     *  defined in the parent object.
     */
    //@{

    virtual void getGibbs_RT_ref(doublereal* grt) const;
    virtual void getStandardVolumes_ref(doublereal* vol) const;

    //@}
    /*! @name Initialization Methods - For Internal use
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally seen by application programs
     */
    //@{

public:
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //!  Create and install a constant volume pressure dependent
    //!  standard state for one species within this object
    /*!
     *  This function sets up the internal data within this object for
     *  handling the calculation of the standard state for the species.
     *
     *  -   It registers the species with the MultiSpeciesThermo object for the
     *      containing VPStandardStateTP phase.
     *  -   It grabs the molar volume property and installs its value within
     *      this object.
     *  -   It also creates a PDSS object, which basically contains a
     *      duplication of some of this information and returns a pointer to
     *      the new object.
     *
     *  @param k Species index within the phase
     *  @param speciesNode Reference to the species node in the XML tree
     *  @param phaseNode_ptr Pointer to the phase node in the XML tree
     *  @return A pointer to the a newly created PDSS object containing the
     *          parameterization
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);
    //@}

    virtual PDSS_enumType reportPDSSType(int index = -1) const;
    virtual VPSSMgr_enumType reportVPSSMgrType() const;
};

}

#endif
