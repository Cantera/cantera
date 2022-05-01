//! @file Interface.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_INTERFACE_H
#define CT_INTERFACE_H

#include "cantera/base/Solution.h"
#include "cantera/thermo/SurfPhase.h"
#include "cantera/kinetics/InterfaceKinetics.h"

namespace Cantera
{

//! A container class holding managers for all pieces defining an interface
class Interface : public Solution
{
private:
    Interface();

public:
    ~Interface() {}
    Interface(const Interface&) = delete;
    Interface& operator=(const Interface&) = delete;

    //! Create an empty Interface object
    static shared_ptr<Interface> create() {
        return shared_ptr<Interface>(new Interface());
    }

    //! Set the reacting phase thermo object
    void setThermo(shared_ptr<ThermoPhase> thermo) override;

    //! Set the Kinetics object
    void setKinetics(shared_ptr<Kinetics> kinetics) override;

    //! Get the surface phase thermo object
    shared_ptr<SurfPhase> thermo() {
        return m_surf;
    }

    //! Get the surface phase Kinetics object
    shared_ptr<InterfaceKinetics> kinetics() {
        return m_surfkin;
    }

protected:
    shared_ptr<SurfPhase> m_surf; //!< Surface phase ThermoPhase manager
    shared_ptr<InterfaceKinetics> m_surfkin;  //!< Kinetics manager
};

//! Create and initialize a new Interface from an input file
/*!
 * This constructor wraps newPhase() and newKinetics()
 *
 * @param infile name of the input file
 * @param name   name of the surface phase in the file.
 *               If this is blank, the first phase in the file is used.
 * @param adjacent vector containing names of adjacent phases that participate in this
 *                 phases kinetics. If empty, adjacent phases will be instantiated based
 *                 on the phase definition.
 * @returns an initialized Interface object.
 */
shared_ptr<Interface> newInterface(const std::string& infile,
    const std::string& name="", const std::vector<std::string>& adjacent={});


//! Create and initialize a new Interface from an input file
/*!
 * This constructor wraps newPhase() and newKinetics()
 *
 * @param infile name of the input file
 * @param name   name of the phase in the file. If this is the empty string, the first
 *               phase in the file is used.
 * @param adjacent vector containing adjacent Solution objects. If empty, adjacent
 *                 phases will be instantiated based on the phase definition.
 * @returns an initialized Interface object.
 */
shared_ptr<Interface> newInterface(const std::string& infile,
    const std::string& name, const std::vector<shared_ptr<Solution>>& adjacent);

//! Create and initialize a new Interface from AnyMap objects
/*!
 * This constructor wraps newPhase() and newKinetics()
 *
 * @param phaseNode the node containing the phase definition (that is, thermo model,
 *     list of species, and initial state)
 * @param rootNode the root node of the tree containing the phase definition, which
 *     will be used as the default location from which to read species definitions.
 * @param adjacent vector containing adjacent Solution objects. If empty, adjacent
 *                 phases will be instantiated based on the phase definition.
 * @returns an initialized Interface object.
 */
shared_ptr<Interface> newInterface(AnyMap& phaseNode, const AnyMap& rootNode=AnyMap(),
    const std::vector<shared_ptr<Solution>>& adjacent={});
}

#endif
