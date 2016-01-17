/**
 *  @file VPSSMgr_Water_ConstVol.h
 * Declaration file for a derived class that handles the calculation
 * of standard state thermo properties for real water and
 *  a set of species which have a constant molar volume pressure
 * dependence
 * (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr_ConstVol VPSSMgr_ConstVol\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_VPSSMGR_WATER_CONSTVOL_H
#define CT_VPSSMGR_WATER_CONSTVOL_H

#include "VPSSMgr.h"

namespace Cantera
{
class PDSS_Water;

//! Handles the calculation of standard state thermo properties for real water
//! and a set of species which have a constant molar volume pressure
//! dependence.
class VPSSMgr_Water_ConstVol : public VPSSMgr
{
public:
    //! Base Constructor
    /*!
     * Initialize the object.
     *
     *  @param vp_ptr   Pointer to the VPStandardStateTP standard state
     *  @param sp_ptr   Pointer to the SpeciesThermo standard state
     */
    VPSSMgr_Water_ConstVol(VPStandardStateTP* vp_ptr, SpeciesThermo* sp_ptr);

    //! Copy Constructor
    VPSSMgr_Water_ConstVol(const VPSSMgr_Water_ConstVol& right);

    //! Assignment operator
    VPSSMgr_Water_ConstVol& operator=(const VPSSMgr_Water_ConstVol& right);

    virtual VPSSMgr* duplMyselfAsVPSSMgr() const;

private:
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

public:

    /*! @name Thermodynamic Values for the Species Reference States
     * There are also temporary variables for holding the species reference-
     * state values of Cp, H, S, and V at the last temperature and reference
     * pressure called. These functions are not recalculated if a new call is
     * made using the previous temperature. All calculations are done within
     * the routine  _updateRefStateThermo().
     */
    //@{

    virtual void getEnthalpy_RT_ref(doublereal* hrt) const;
    virtual void getGibbs_RT_ref(doublereal* grt) const ;
    virtual void getGibbs_ref(doublereal* g) const ;
    virtual void getEntropy_R_ref(doublereal* er) const ;
    virtual void getCp_R_ref(doublereal* cpr) const ;
    virtual void getStandardVolumes_ref(doublereal* vol) const ;

    //! @}
    /*! @name Initialization Methods - For Internal use
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see importPhase().
     */
    //@{

    virtual void initThermo();
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);
    //@}

    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);

    virtual PDSS_enumType reportPDSSType(int index = -1) const ;
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;
    virtual void initAllPtrs(VPStandardStateTP* vp_ptr, SpeciesThermo* sp_ptr);

private:
    //! Pointer to the Water PDSS object.
    /*!
     * This is a shallow copy. The water PDSS object is owned by the VPStandardStateTP
     * object.
     */
    PDSS_Water* m_waterSS;
};
}

#endif
