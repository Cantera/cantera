/**
 *  @file VPSSMgr_IdealGas.h
 * Declaration file for a derived class that handles the calculation
 * of standard state thermo properties for
 *  a set of species which have an Ideal Gas dependence
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr_IdealGas VPSSMgr_IdealGas\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_VPSSMGR_IDEALGAS_H
#define CT_VPSSMGR_IDEALGAS_H

#include "VPSSMgr.h"

namespace Cantera
{
//! A VPSSMgr where all species in the phase obey an ideal gas equation of state
/**
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 */
class VPSSMgr_IdealGas : public VPSSMgr
{
public:
    //! Basic constructor that initializes the object
    /*!
     * @param vp_ptr Pointer to the owning ThermoPhase
     * @param spth   Species thermo pointer.
     */
    VPSSMgr_IdealGas(VPStandardStateTP* vp_ptr, MultiSpeciesThermo* spth);

    VPSSMgr_IdealGas(const VPSSMgr_IdealGas& right);
    VPSSMgr_IdealGas& operator=(const VPSSMgr_IdealGas& right);
    virtual VPSSMgr* duplMyselfAsVPSSMgr() const;

    /*! @name  Properties of the Standard State of the Species in the Solution
     * Within VPStandardStateTP, these properties are calculated via a common
     * routine, _updateStandardStateThermo(), which must be overloaded in
     * inherited objects. The values are cached within this object, and are
     * not recalculated unless the temperature or pressure changes.
     */
    //@{
    virtual void getIntEnergy_RT(doublereal* urt) const;
    virtual void getStandardVolumes(doublereal* vol) const;
    //@}

protected:
    virtual void _updateStandardStateThermo();

public:
    //! Create and install an ideal gas standard state manager for one species
    //! within this object
    /*!
     *  This function sets up the internal data within this object for
     *  handling the calculation of the standard state for the species.
     *
     *  -   It registers the species with the MultiSpeciesThermo object for the
     *      containing VPStandardStateTP phase.
     *  -   It also creates a PDSS object, which basically contains a
     *      duplication of some of this information and returns a pointer to
     *      the new object.
     *  .
     *  @param k             Species index within the phase
     *  @param speciesNode   Reference to the species node in the XML tree
     *  @param phaseNode_ptr Pointer to the phase node in the XML tree
     *  @return a pointer to the a newly created PDSS object containing the
     *          parameterization
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);

    virtual PDSS_enumType reportPDSSType(int index = -1) const;
    virtual VPSSMgr_enumType reportVPSSMgrType() const;
};
}

#endif
