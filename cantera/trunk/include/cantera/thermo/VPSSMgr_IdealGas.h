/**
 *  @file VPSSMgr_IdealGas.h
 * Declaration file for a derived class that handles the calculation
 * of standard state thermo properties for
 *  a set of species which have an Ideal Gas dependence
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr_IdealGas VPSSMgr_IdealGas\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_VPSSMGR_IDEALGAS_H
#define CT_VPSSMGR_IDEALGAS_H

#include "VPSSMgr.h"

namespace Cantera
{
//! A VPSSMgr where all species in the phase obey an ideal gas equation of state
class VPSSMgr_IdealGas : public VPSSMgr
{
public:
    //! Basic constructor that initializes the object
    /*!
     * @param vp_ptr Pointer to the owning ThermoPhase
     * @param spth   Species thermo pointer.
     */
    VPSSMgr_IdealGas(VPStandardStateTP* vp_ptr, SpeciesThermo* spth);

    //! Copy Constructor
    VPSSMgr_IdealGas(const VPSSMgr_IdealGas& right);

    //! Assignment operator
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

    /*! @name Initialization Methods - For Internal use
     * The following methods are used in the process of constructing the phase
     * and setting its parameters from a specification in an input file. They
     * are not normally used in application programs. To see how they are
     * used, see importPhase().
     */
    //@{
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);
    //@}

    //! Create and install an ideal gas standard state manager for one species
    //! within this object
    /*!
     *  This function sets up the internal data within this object for
     *  handling the calculation of the standard state for the species.
     *
     *  -   It registers the species with the SpeciesThermo object for the
     *      containing VPStandardStateTP phase.
     *  -   It also creates a PDSS object, which basically contains a
     *      duplication of some of this information and returns a pointer to
     *      the new object.
     *  .
     *  @param k             Species index within the phase
     *  @param speciesNode   Reference to the species node in the XML tree
     *  @param phaseNode_ptr Pointer to the phase node in the XML tree
     *  @return Returns a pointer to the a newly malloced PDSS object
     *          containing the parameterization
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);

    virtual PDSS_enumType reportPDSSType(int index = -1) const ;
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;
};
}

#endif
