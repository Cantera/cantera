//! @file ReactorNode.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_REACTORNODE_H
#define CT_REACTORNODE_H

#include "cantera/base/global.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

class Solution;

/**
 * Base class for reactor nodes.
 * @ingroup reactorGroup
 */
class ReactorNode
{
protected:
    //! Instantiate a ReactorNode object with Solution contents.
    //! @param sol  Solution object to be set.
    //! @param name  Name of the reactor.
    //! @since New in %Cantera 3.1.
    ReactorNode(shared_ptr<Solution> sol, const string& name="(none)");

public:
    ReactorNode() = default;
    virtual ~ReactorNode();
    ReactorNode(const ReactorNode&) = delete;
    ReactorNode& operator=(const ReactorNode&) = delete;

    //! String indicating the reactor model implemented. Usually
    //! corresponds to the name of the derived class.
    virtual string type() const {
        return "ReactorNode";
    }

    //! Return the name of this reactor
    string name() const {
        return m_name;
    }

    //! Set the name of this reactor
    void setName(const string& name) {
        m_name = name;
    }

    //! Set the default name of a reactor. Returns `false` if it was previously set.
    void setDefaultName(map<string, int>& counts);

    //! Set the state of the Phase object associated with this reactor to the
    //! reactor's current state.
    virtual void restoreState();

    //! Set the state of the reactor to correspond to the state of the
    //! associated ThermoPhase object. This is the inverse of restoreState().
    //! Calling this will trigger integrator reinitialization.
    virtual void syncState();

protected:
    //! Composite thermo/kinetics handler.
    shared_ptr<Solution> m_solution;

    string m_name;  //!< Reactor name.
    bool m_defaultNameSet = false;  //!< `true` if default name has been previously set.
};
}

#endif
