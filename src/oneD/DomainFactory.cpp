//! @file DomainFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/oneD/DomainFactory.h"
#include "cantera/oneD/Boundary1D.h"
#include "cantera/oneD/Flow1D.h"
#include "cantera/oneD/IonFlow.h"
#include "cantera/transport/Transport.h"

namespace Cantera
{

DomainFactory* DomainFactory::s_factory = 0;
std::mutex DomainFactory::domain_mutex;

DomainFactory::DomainFactory()
{
    reg("inlet", [](shared_ptr<Solution> solution, const string& id) {
        return new Inlet1D(solution, id);
    });
    reg("empty", [](shared_ptr<Solution> solution, const string& id) {
        return new Empty1D(solution, id);
    });
    reg("symmetry-plane", [](shared_ptr<Solution> solution, const string& id) {
        return new Symm1D(solution, id);
    });
    reg("outlet", [](shared_ptr<Solution> solution, const string& id) {
        return new Outlet1D(solution, id);
    });
    reg("outlet-reservoir", [](shared_ptr<Solution> solution, const string& id) {
        return new OutletRes1D(solution, id);
    });
    reg("surface", [](shared_ptr<Solution> solution, const string& id) {
        return new Surf1D(solution, id);
    });
    reg("reacting-surface", [](shared_ptr<Solution> solution, const string& id) {
        return new ReactingSurf1D(solution, id);
    });
    reg("gas-flow", [](shared_ptr<Solution> solution, const string& id) {
        return new Flow1D(solution, id);
    });
    reg("ion-flow", [](shared_ptr<Solution> solution, const string& id) {
        return new IonFlow(solution, id);
    });
    reg("free-flow", [](shared_ptr<Solution> solution, const string& id) {
        Flow1D* ret;
        if (solution->transport()->transportModel() == "ionized-gas") {
            ret = new IonFlow(solution, id);
        } else {
            ret = new Flow1D(solution, id);
        }
        ret->setFreeFlow();
        return ret;
    });
    reg("axisymmetric-flow", [](shared_ptr<Solution> solution, const string& id) {
        Flow1D* ret;
        if (solution->transport()->transportModel() == "ionized-gas") {
            ret = new IonFlow(solution, id);
        } else {
            ret = new Flow1D(solution, id);
        }
        ret->setAxisymmetricFlow();
        return ret;
    });
    reg("unstrained-flow", [](shared_ptr<Solution> solution, const string& id) {
        Flow1D* ret;
        if (solution->transport()->transportModel() == "ionized-gas") {
            ret = new IonFlow(solution, id);
        } else {
            ret = new Flow1D(solution, id);
        }
        ret->setUnstrainedFlow();
        return ret;
    });
}

DomainFactory* DomainFactory::factory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    if (!s_factory) {
        s_factory = new DomainFactory;
    }
    return s_factory;
}

void DomainFactory::deleteFactory()
{
    std::unique_lock<std::mutex> lock(domain_mutex);
    delete s_factory;
    s_factory = 0;
}

}
