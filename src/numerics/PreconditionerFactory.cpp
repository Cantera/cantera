//! @file PreconditionerFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/PreconditionerFactory.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/oneD/MultiJac.h"

namespace Cantera
{

PreconditionerFactory* PreconditionerFactory::factory() {
    std::unique_lock<std::mutex> lock(precon_mutex);
    if (!s_factory) {
        s_factory = new PreconditionerFactory;
    }
    return s_factory;
};

//! Delete preconditioner factory
void PreconditionerFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(precon_mutex);
    delete s_factory;
    s_factory = 0;
};

PreconditionerFactory* PreconditionerFactory::s_factory = 0;
std::mutex PreconditionerFactory::precon_mutex;

PreconditionerFactory::PreconditionerFactory()
{
    reg("Adaptive", []() { return new AdaptivePreconditioner(); });
    reg("banded-direct", []() { return new MultiJac(); });
}

shared_ptr<PreconditionerBase> newPreconditioner(const string& precon)
{
    return shared_ptr<PreconditionerBase>(PreconditionerFactory::factory()->create(precon));
};

}
