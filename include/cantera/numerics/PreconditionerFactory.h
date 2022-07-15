//! @file PreconditionerFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef PRECONDITIONER_FACTORY_H
#define PRECONDITIONER_FACTORY_H

#include "cantera/base/FactoryBase.h"

namespace Cantera
{

class PreconditionerBase;

//! Factory class to create preconditioner objects
class PreconditionerFactory : public Factory<PreconditionerBase>
{
public:
    static PreconditionerFactory* factory() {
        std::unique_lock<std::mutex> lock(precon_mutex);
        if (!s_factory) {
            s_factory = new PreconditionerFactory;
        }
        return s_factory;
    };

    //! Delete preconditioner factory
    virtual void deleteFactory() {
        std::unique_lock<std::mutex> lock(precon_mutex);
        delete s_factory;
        s_factory = 0;
    };

private:
    static PreconditionerFactory* s_factory;
    static std::mutex precon_mutex;
    PreconditionerFactory();
};

//! Create a Preconditioner object of the specified type
inline std::shared_ptr<PreconditionerBase> newPreconditioner(const std::string& precon)
{
    return std::shared_ptr<PreconditionerBase>(PreconditionerFactory::factory()->create(precon));
};

}

#endif
