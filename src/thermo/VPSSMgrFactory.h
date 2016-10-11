/**
 *  @file VPSSMgrFactory.h
 *    Header for factory to build instances of classes that manage the
 *    standard-state thermodynamic properties of a set of species
 *    (see \ref mgrpdssthermocalc and class \link Cantera::VPSSMgrFactory VPSSMgrFactory\endlink);
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef VPSSSPECIESTHERMO_FACTORY_H
#define VPSSSPECIESTHERMO_FACTORY_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/thermo/VPSSMgr.h"

namespace Cantera
{
//! Throw a named error for an unknown or missing vpss species thermo model.
/*!
 * @ingroup mgrpdssthermocalc
 */
class UnknownVPSSMgrModel: public CanteraError
{
public:
    //! Constructor
    /*!
     * @param proc                    Function name error occurred.
     * @param VPSSMgrModel            Unrecognized species thermo calculator name
     * @deprecated Unused. To be removed after Cantera 2.3.
     */
    UnknownVPSSMgrModel(const std::string& proc,
                        const std::string& VPSSMgrModel) :
        CanteraError(proc, "Specified VPSSMgr model "
                     + VPSSMgrModel +
                     " does not match any known type.") {}
};

//! Factory to build instances of classes that manage the
//! standard-state thermodynamic properties of a set of species.
/*!
 *  This class is responsible for making the decision concerning which
 *  derivative of VPSSMgr object to use. The VPSSMgr object is used to calculate
 *  thermodynamic functions for the standard state. It queries the database of
 *  species to understand what the requirements are for the submodels for all of
 *  the species in the phase. Then, it picks the derived VPSSMgr object to use
 *  and passes it back to the calling routine. It doesn't load any data into the
 *  derived VPSSMgr object.
 *
 *  Making the choice of VPSSMgr types is the only thing this class does.
 *
 * This class is implemented as a singleton -- one in which only one instance is
 * needed.  The recommended way to access the factory is to call this static
 * method, which instantiates the class if it is the first call, but otherwise
 * simply returns the pointer to the existing instance.
 *
 * @ingroup mgrpdssthermocalc
 */
class VPSSMgrFactory : public Factory<VPSSMgr, VPStandardStateTP*, MultiSpeciesThermo*>
{
public:
    //! Static method to return an instance of this class
    /*!
     * This class is implemented as a singleton -- one in which only one
     * instance is needed.  The recommended way to access the factory is to call
     * this static method, which instantiates the class if it is the first call,
     * but otherwise simply returns the pointer to the existing instance.
     */
    static VPSSMgrFactory* factory() {
        std::unique_lock<std::mutex> lock(vpss_species_thermo_mutex);
        if (!s_factory) {
            s_factory = new VPSSMgrFactory;
        }
        return s_factory;
    }

    //! Delete static instance of this class
    /*!
     * If it is necessary to explicitly delete the factory before the process
     * terminates (for example, when checking for memory leaks) then this method
     * can be called to delete it.
     */
    void deleteFactory();

    //! String conversion to an enumType
    /*!
     *  This routine is a string conversion. The string is obtained from the
     *  standardState model attribute and converted to a VPSSMgr_enumType type.
     *
     * @param ssModel  String representing the VPSSMGr object
     * @deprecated Unused. To be removed after Cantera 2.3.
     */
    virtual VPSSMgr_enumType
    VPSSMgr_StringConversion(const std::string& ssModel) const;

    //! Create a new species variable pressure standard state calculator
    /*!
     * @param type    The enumerated type of the standard state calculator
     * @param vp_ptr  Variable pressure standard state ThermoPhase object
     *                that will be the owner.
     */
    virtual VPSSMgr* newVPSSMgr(VPSSMgr_enumType type, VPStandardStateTP* vp_ptr);

    //! Create a new species property manager for a group of species
    /*!
     * This routine will look through species nodes. It will discover what each
     * species needs for its species property managers. Then, it will create and
     * return the proper species property manager to use.
     *
     * @param vp_ptr         Variable pressure standard state ThermoPhase object
     *                       that will be the owner.
     * @param phaseNode_ptr  Pointer to the ThermoPhase phase XML Node
     * @param spDataNodeList Vector of XML_Nodes, each of which is a species XML
     *                       Node. There are m_kk of these.
     * @returns              a pointer to a newly created species
     *                       property manager object.
     */
    virtual VPSSMgr* newVPSSMgr(VPStandardStateTP* vp_ptr,
                                XML_Node* phaseNode_ptr,
                                std::vector<XML_Node*> & spDataNodeList);

private:
    //! pointer to the sole instance of this class
    static VPSSMgrFactory* s_factory;

    //! Decl of the static mutex variable that locks the
    //! VPSSMgr factory singleton
    static std::mutex vpss_species_thermo_mutex;

    //! Constructor. This is made private, so that only the static
    //! method factory() can instantiate the class.
    VPSSMgrFactory();
};

////////////////////// Convenience functions ////////////////////
//  These functions allow using a different factory class that
//  derives from VPSSMgrFactory.
//////////////////////////////////////////////////////////////////

//! Create a new species thermo manager instance, by specifying
//! the type and (optionally) a pointer to the factory to use to create it.
/*!
 * This utility program will look through species nodes. It will discover what
 * each species needs for its species property managers. Then, it will create
 * and return the proper species property manager to use.
 *
 *  These functions allow using a different factory class that
 *  derives from VPSSMgrFactory.
 *
 * @param type         Species thermo type.
 * @param vp_ptr       Variable pressure standard state ThermoPhase object
 *                     that will be the owner.
 * @param f            Pointer to a VPSSMgrFactory. optional parameter.
 *                     Defaults to NULL.
 * @deprecated Unused. To be removed after Cantera 2.3.
 */
VPSSMgr* newVPSSMgr(VPSSMgr_enumType type,
                    VPStandardStateTP* vp_ptr, VPSSMgrFactory* f=0);

//! Function to return VPSSMgr manager
/*!
 * This utility program will look through species nodes. It will discover what
 * each species needs for its species property managers. Then, it will alloc
 * and return the proper species property manager to use.
 *
 *  These functions allow using a different factory class that
 *  derives from VPSSMgrFactory.
 *
 * @param vp_ptr         Variable pressure standard state ThermoPhase object
 *                       that will be the owner.
 * @param phaseNode_ptr  Pointer to the ThermoPhase phase XML Node
 * @param spDataNodeList This vector contains a list
 *                       of species XML nodes that will be in the phase
 * @param f              Pointer to a VPSSMgrFactory. optional
 *                       parameter. Defaults to NULL.
 * @deprecated The `VPSSMgrFactory*` argument to this function is deprecated and
 *     will be removed after Cantera 2.3.
 */
VPSSMgr* newVPSSMgr(VPStandardStateTP* vp_ptr,
                    XML_Node* phaseNode_ptr,
                    std::vector<XML_Node*> & spDataNodeList,
                    VPSSMgrFactory* f=0);
}

#endif
