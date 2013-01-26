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

#include "cantera/base/ct_defs.h"
#include "PDSS.h"
#include "VPSSMgr.h"

namespace Cantera
{

template<typename ValAndDerivType> class SpeciesThermoInterpType;
class VPStandardStateTP;
template<typename ValAndDerivType> class SpeciesThermo;


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
class VPSSMgr_IdealGas : public VPSSMgr
{

public:


    //! Basic constructor that initializes the object
    /*!
     * @param vp_ptr Pointer to the owning ThermoPhase
     * @param spth   Species thermo pointer.
     */
    VPSSMgr_IdealGas(VPStandardStateTP* vp_ptr, SpeciesThermo<doublereal> * spth);

    //! Destructor
    virtual ~VPSSMgr_IdealGas();

    //! Copy Constructor for the %SpeciesThermo object.
    /*!
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    VPSSMgr_IdealGas(const VPSSMgr_IdealGas& right);

    //! Assignment operator for the %SpeciesThermo object
    /*!
     *  This is NOT a virtual function.
     *
     * @param right    Reference to %SpeciesThermo object to be copied into the
     *                 current one.
     */
    VPSSMgr_IdealGas& operator=(const VPSSMgr_IdealGas& right);

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

    /**
     *  Returns the vector of nondimensional
     *  internal Energies of the standard state at the current temperature
     *  and pressure of the solution for each species.
     * \f[
     *  u^{ss}_k(T,P) = h^{ss}_k(T)  - P * V^{ss}_k
     * \f]
     *
     * @param urt    Output vector of nondimensional standard state
     *               internal energies. length = m_kk.
     */
    virtual void getIntEnergy_RT(doublereal* urt) const;

    /**
     * Get the molar volumes of each species in their standard
     * states at the current
     * <I>T</I> and <I>P</I> of the solution.
     * units = m^3 / kmol
     *
     * This is redefined here to call the internal function,  _updateStandardStateThermo(),
     * which calculates all standard state properties at the same time.
     *
     * @param vol Output vector of species volumes. length = m_kk.
     *            units =  m^3 / kmol
     */
    virtual void getStandardVolumes(doublereal* vol) const;

protected:

    //! Updates the standard state thermodynamic functions at the current
    //! T and P of the solution.
    /*!
     * @internal
     *
     * If m_useTmpStandardStateStorage is true,
     * this function must be called every time the temperature or pressure
     * has changed.
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
     *  required to be called everytime this class is invoked.
     *
     */
    virtual void _updateStandardStateThermo();

public:

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



    //! @name Initialization Methods - For Internal use (VPStandardState)
    /*!
     * The following methods are used in the process of constructing
     * the phase and setting its parameters from a specification in an
     * input file. They are not normally used in application programs.
     * To see how they are used, see files importCTML.cpp and
     * ThermoFactory.cpp.
     *
     */
    //@{

    //! Initialize the thermo for this standard state thermo calculator
    /*!
     *  This task is done last, after createInstallPDSS() and after
     *  initThermo().
     *
     *  @param phaseNode   Reference to the phase node in the XML tree
     *  @param id          string name of the phase
     */
    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id);

    //!  Create and install an ideal gas standard state manager
    //!  for one species within this object
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
     *
     *  @return Returns a pointer to the a newly malloced PDSS object
     *          containing the parameterization
     */
    virtual PDSS* createInstallPDSS(size_t k, const XML_Node& speciesNode,
                                    const XML_Node* const phaseNode_ptr);


    //! This utility function reports the type of parameterization
    //! used for the species with index number index.
    /*!
     *
     * @param index  Species index
     */
    virtual PDSS_enumType reportPDSSType(int index = -1) const ;


    //! This utility function reports the type of manager
    //! for the calculation of standard state properties
    /*!
     *
     */
    virtual VPSSMgr_enumType reportVPSSMgrType() const ;

};
//@}
}

#endif

