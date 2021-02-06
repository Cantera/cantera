//! @file ReactorFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef REACTOR_FACTORY_H
#define REACTOR_FACTORY_H

#include "cantera/zeroD/ReactorBase.h"
#include "cantera/base/FactoryBase.h"

namespace Cantera
{

class ReactorFactory : public Factory<ReactorBase>
{
public:
    static ReactorFactory* factory() {
        std::unique_lock<std::mutex> lock(reactor_mutex);
        if (!s_factory) {
            s_factory = new ReactorFactory;
        }
        return s_factory;
    }

    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(reactor_mutex);
        delete s_factory;
        s_factory = 0;
    }

    //! Create a new reactor by type identifier.
    /*!
     * @param n the type to be created.
     */
    virtual ReactorBase* newReactor(int n);

    //! Create a new reactor by type name.
    /*!
     * @param reactorType the type to be created.
     */
    virtual ReactorBase* newReactor(const std::string& reactorType);

    //! Register a new reactor type identifier.
    /*!
     * @param name the name of the reactor type.
     * @param type the type identifier of the reactor.
     * Integer type identifiers are used by clib and matlab interfaces.
     *
     * @deprecated To be removed after Cantera 2.5.
     */
    void reg_type(const std::string& name, const int type) {
        m_types[type] = name;
    }

protected:
    //! Map containing reactor type identifier / reactor type name pairs.
    //! @deprecated To be removed after Cantera 2.5.
    std::unordered_map<int, std::string> m_types;

private:
    static ReactorFactory* s_factory;
    static std::mutex reactor_mutex;
    ReactorFactory();
};

//! Create a Reactor object of the specified type
inline ReactorBase* newReactor(const std::string& model)
{
    return ReactorFactory::factory()->newReactor(model);
}

}

#endif
