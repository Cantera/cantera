//! @file SystemJacobianFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/SystemJacobianFactory.h"
#include "cantera/numerics/AdaptivePreconditioner.h"
#include "cantera/oneD/MultiJac.h"

namespace Cantera
{

SystemJacobianFactory* SystemJacobianFactory::factory() {
    std::unique_lock<std::mutex> lock(jac_mutex);
    if (!s_factory) {
        s_factory = new SystemJacobianFactory;
    }
    return s_factory;
};

void SystemJacobianFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(jac_mutex);
    delete s_factory;
    s_factory = 0;
};

SystemJacobianFactory* SystemJacobianFactory::s_factory = 0;
std::mutex SystemJacobianFactory::jac_mutex;

SystemJacobianFactory::SystemJacobianFactory()
{
    reg("Adaptive", []() { return new AdaptivePreconditioner(); });
    reg("banded-direct", []() { return new MultiJac(); });
}

shared_ptr<SystemJacobian> newSystemJacobian(const string& type)
{
    return shared_ptr<SystemJacobian>(SystemJacobianFactory::factory()->create(type));
};

}
