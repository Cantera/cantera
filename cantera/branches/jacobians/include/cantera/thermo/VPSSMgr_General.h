/**
 *  @file VPSSMgr_General.h
 *  Declaration file for a derived class that handles the calculation
 * of standard state thermo properties for
 * a set of species belonging to a single phase in a completely general
 * but slow way (see \ref mgrpdssthermocalc and
 * class \link Cantera::VPSSMgr_General VPSSMgr_General\endlink).
 */
/*
 * Copyright (2007) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_VPSSMGR_GENERAL_H
#define CT_VPSSMGR_GENERAL_H

#include "cantera/base/ct_defs.h"
#include "VPSSMgr.h"

namespace Cantera
{

template<typename ValAndDerivType> class SpeciesThermoInterpType;
class VPStandardStateTP;
template<typename ValAndDerivType> class SpeciesThermo;
class PDSS;


//!  Class that handles the calculation of standard state thermo properties for
//!  a set of species belonging to a single phase in a completely general
//!  but slow way.
/*!
 *   This class manages the calculation of standard state thermo properties for
 *   a set of species belonging to a single phase in a completely general
 *   but slow way.
 *   The way this does this is to call the underlying PDSS routines one at a
 *   time for every species.
 *
 * @ingroup mgrpdssthermocalc
 */
class VPSSMgr_General : public VPSSMgr
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
    VPSSMgr_General(VPStandardStateTP* vp_ptr,
                    SpeciesThermo<doublereal> * spth);

    //! Destructor
    virtual ~VPSSMgr_General();

    //! Copy Constructor for the %SpeciesThermo object.
    /*!
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    VPSSMgr_General(const VPSSMgr_General& right);

    //! Assignment operator for the %SpeciesThermo object
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    VPSSMgr_General& operator=(const VPSSMgr_General& right);

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


protected:

    //! Internally updates the standard state thermodynamic functions at the current
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

    //! Finalize the thermo objects after all species have been entered
    /*!
     *  This function is the LAST initialization routine to be
     *  called. It's called after createInstallPDSS() has been
     *  called for each species in the phase, and after initThermo()
     *  has been called.
     *  It's called via an inner-to-outer onion-shell like manner.
     *
     *  Currently, this routine passed control to the parent class
     *  without doing anything.
     *
     *  @param phaseNode   Reference to the phaseNode XML node.
     *  @param id          ID of the phase.
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

private:
    //! Local factory routine for the creation of PDSS objects
    /*!
     *   This routine is specific to the VPSSMgr_General object.
     *   It will create a PDSS object for species k, by searching
     *   and querying for the "standardState" XML node in the standard
     *   state description of the species. If this XML node doesn't
     *   exist, it will assume that the standard state is an ideal
     *   gas.
     *   It decides on the attribute, "model", what PDSS object
     *   to create.
     *
     *   @param k  Species number
     *   @param speciesNode XML node for the standard state of the species
     *   @param phaseNode_ptr   pointer to the phase XML node
     *   @param doST  output variable indicating whether the
     *             instantiation has resulted in a SpeciesThermo object
     *             being created and registered with the SpeciesThermo
     *             manager class.
     *
     *   @return  Returns the pointer to a malloced PDSS object
     */
    PDSS* returnPDSS_ptr(size_t k, const XML_Node& speciesNode,
                         const XML_Node* const phaseNode_ptr, bool& doST);

public:

    //! Factory routine for the creation of PDSS objects that are
    //! then internally registered with this VPSSMgr object
    /*!
     *  This function sets up the internal data within this object for
     *  handling the calculation of the standard state for the species.
     *
     *   This routine
     *   will create a PDSS object for species k, by searching
     *   and querying for the "standardState" XML node in the standard
     *   state description of the species.
     *   It will then store the object's pointer in a vector of pointers,
     *   and it will own the object.
     *
     *   @param k  Species number
     *   @param speciesNode XML node for the standard state of the species
     *   @param phaseNode_ptr   pointer to the phase XML node
     *
     *   @return  Returns the pointer to the malloced PDSS object
     */
    virtual  PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                     const XML_Node* const phaseNode_ptr);

    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param index  Species index
     */
    virtual PDSS_enumType reportPDSSType(int index = -1) const ;


    //! This utility function reports the type of manager
    //! for the calculation of the standard state properties
    /*!
     *
     */
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;

    //! Initialize the internal shallow pointers in this object
    /*!
     * There are a bunch of internal shallow pointers that point to the owning
     * VPStandardStateTP and SpeciesThermo objects. This function reinitializes
     * them. This function is called like an onion.
     *
     *  @param vp_ptr   Pointer to the VPStandardStateTP standard state
     *  @param sp_ptr   Pointer to the SpeciesThermo standard state
     */
    virtual void initAllPtrs(VPStandardStateTP* vp_ptr, SpeciesThermo<doublereal> * sp_ptr);

private:

    //! Shallow pointers containing the PDSS objects for the species
    //! in this phase.
    /*!
     * This object doesn't own these pointers.
     */
    std::vector<PDSS*> m_PDSS_ptrs;
};
//@}
}

#endif

