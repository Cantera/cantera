//! @file SystemJacobianFactory.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef SYSTEMJACOBIANFACTORY_FACTORY_H
#define SYSTEMJACOBIANFACTORY_FACTORY_H

#include "cantera/base/FactoryBase.h"
#include "cantera/base/global.h"

namespace Cantera
{

class SystemJacobian;

//! Factory class to create Jacobian objects for use by linear solvers
class SystemJacobianFactory : public Factory<SystemJacobian>
{
public:
    static SystemJacobianFactory* factory();
    void deleteFactory() override;

private:
    static SystemJacobianFactory* s_factory;
    static std::mutex jac_mutex;
    SystemJacobianFactory();
};

//! Create a SystemJacobian object of the specified type
shared_ptr<SystemJacobian> newSystemJacobian(const string& type);

//! @deprecated To be removed after %Cantera 3.2. Renamed to newSystemJacobian()
inline shared_ptr<SystemJacobian> newPreconditioner(const string& type) {
    warn_deprecated("newPreconditioner",
                    "To be removed after Cantera 3.2. Renamed to newSystemJacobian");
    return newSystemJacobian(type);
}

}

#endif
