//! @file PreconditionerFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/numerics/PreconditionerFactory.h"
#include "cantera/numerics/AdaptivePreconditioner.h"


using namespace std;
namespace Cantera
{

PreconditionerFactory* PreconditionerFactory::s_factory = 0;
std::mutex PreconditionerFactory::precon_mutex;

PreconditionerFactory::PreconditionerFactory()
{
    reg("Adaptive", []() { return new AdaptivePreconditioner(); });
}

}
