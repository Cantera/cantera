/**
 *  @file VPSSMgr_General.h
 *  Declaration file for a derived class that handles the calculation
 * of standard state thermo properties for
 * a set of species belonging to a single phase in a completely general
 * but slow way (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr_General VPSSMgr_General\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_VPSSMGR_GENERAL_H
#define CT_VPSSMGR_GENERAL_H

#include "VPSSMgr.h"

namespace Cantera
{
//!  Class that handles the calculation of standard state thermo properties for
//!  a set of species belonging to a single phase in a completely general
//!  but slow way.
/*!
 *   This class manages the calculation of standard state thermo properties
 *   for a set of species belonging to a single phase in a completely general
 *   but slow way. The way this does this is to call the underlying PDSS
 *   routines one at a time for every species.
 *
 * @ingroup mgrpdssthermocalc
 */
class VPSSMgr_General : public VPSSMgr
{
public:
    //! Constructor
    /*!
     *  @param vp_ptr Pointer to the owning VPStandardStateTP object for the
     *                phase.
     *  @param spth   Pointer to the MultiSpeciesThermo object for the phase.
     */
    VPSSMgr_General(VPStandardStateTP* vp_ptr,
                    MultiSpeciesThermo* spth);

    VPSSMgr_General(const VPSSMgr_General& right);
    VPSSMgr_General& operator=(const VPSSMgr_General& right);
    virtual VPSSMgr* duplMyselfAsVPSSMgr() const;

protected:
    /*!
     * @name  Properties of the Standard State of the Species in the Solution
     *
     * Within VPStandardStateTP, these properties are calculated via a common
     * routine, _updateStandardStateThermo(), which must be overloaded in
     * inherited objects. The values are cached within this object, and are
     * not recalculated unless the temperature or pressure changes.
     */
    //@{
    virtual void _updateStandardStateThermo();
    virtual void _updateRefStateThermo() const;
    //@}

    /*! @name Thermodynamic Values for the Species Reference States
     * There are also temporary variables for holding the species reference-
     * state values of Cp, H, S, and V at the last temperature and reference
     * pressure called. These functions are not recalculated if a new call is
     * made using the previous temperature. All calculations are done within
     * the routine _updateRefStateThermo().
     */
    //@{
    virtual void getGibbs_ref(doublereal* g) const;
    //@}

    /*! @name Initialization Methods - For Internal use
     * The following methods are used in the process of constructing the phase
     * and setting its parameters from a specification in an input file. They
     * are not normally used in application programs. To see how they are
     * used, see importPhase().
     */
    //@{
    virtual void initThermo();
    //@}

private:
    //! Local factory routine for the creation of PDSS objects
    /*!
     * This routine is specific to the VPSSMgr_General object. It will create
     * a PDSS object for species k, by searching and querying for the
     * "standardState" XML node in the standard state description of the
     * species. If this XML node doesn't exist, it will assume that the
     * standard state is an ideal gas. It decides on the attribute, "model",
     * what PDSS object to create.
     *
     * @param speciesNode XML node for the standard state of the species
     * @param k  Species number
     * @param phaseNode_ptr   pointer to the phase XML node
     * @param doST  output variable indicating whether the
     *           instantiation has resulted in a MultiSpeciesThermo object
     *           being created and registered with the MultiSpeciesThermo
     *           manager class.
     * @returns the pointer to a newly created PDSS object
     */
    PDSS* returnPDSS_ptr(size_t k, const XML_Node& speciesNode,
                         const XML_Node* const phaseNode_ptr, bool& doST);

public:
    //! Factory routine for the creation of PDSS objects that are
    //! then internally registered with this VPSSMgr object
    /*!
     * This function sets up the internal data within this object for handling
     * the calculation of the standard state for the species.
     *
     * This routine will create a PDSS object for species k, by searching and
     * querying for the "standardState" XML node in the standard state
     * description of the species. It will then store the object's pointer in
     * a vector of pointers, and it will own the object.
     *
     * @param k  Species number
     * @param speciesNode XML node for the standard state of the species
     * @param phaseNode_ptr   pointer to the phase XML node
     * @return the pointer to the newly created PDSS object
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                     const XML_Node* const phaseNode_ptr);

    virtual PDSS_enumType reportPDSSType(int index = -1) const;
    virtual VPSSMgr_enumType reportVPSSMgrType() const;
    virtual void initAllPtrs(VPStandardStateTP* vp_ptr, MultiSpeciesThermo* sp_ptr);

private:
    //! Shallow pointers containing the PDSS objects for the species
    //! in this phase. This object doesn't own these pointers.
    std::vector<PDSS*> m_PDSS_ptrs;
};

}

#endif
