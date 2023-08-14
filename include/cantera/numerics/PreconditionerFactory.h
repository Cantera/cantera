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
    static PreconditionerFactory* factory();

    //! Delete preconditioner factory
    void deleteFactory() override;

private:
    static PreconditionerFactory* s_factory;
    static std::mutex precon_mutex;
    PreconditionerFactory();
};

//! Create a Preconditioner object of the specified type
shared_ptr<PreconditionerBase> newPreconditioner(const string& precon);

}

#endif
