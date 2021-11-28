//! @file TransportFactory.cpp Implementation file for class TransportFactory.

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

// known transport models
#include "cantera/transport/MultiTransport.h"
#include "cantera/transport/MixTransport.h"
#include "cantera/transport/UnityLewisTransport.h"
#include "cantera/transport/IonGasTransport.h"
#include "cantera/transport/WaterTransport.h"
#include "cantera/transport/DustyGasTransport.h"
#include "cantera/transport/HighPressureGasTransport.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

using namespace std;

namespace Cantera
{
TransportFactory* TransportFactory::s_factory = 0;

// declaration of static storage for the mutex
std::mutex TransportFactory::transport_mutex;

//////////////////// class TransportFactory methods //////////////

TransportFactory::TransportFactory()
{
    reg("none", []() { return new Transport(); });
    addAlias("none", "Transport");
    addAlias("none", "None");
    addAlias("none", "");
    reg("unity-Lewis-number", []() { return new UnityLewisTransport(); });
    addAlias("unity-Lewis-number", "UnityLewis");
    reg("mixture-averaged", []() { return new MixTransport(); });
    addAlias("mixture-averaged", "Mix");
    reg("mixture-averaged-CK", []() { return new MixTransport(); });
    addAlias("mixture-averaged-CK", "CK_Mix");
    reg("multicomponent", []() { return new MultiTransport(); });
    addAlias("multicomponent", "Multi");
    reg("multicomponent-CK", []() { return new MultiTransport(); });
    addAlias("multicomponent-CK", "CK_Multi");
    reg("ionized-gas", []() { return new IonGasTransport(); });
    addAlias("ionized-gas", "Ion");
    reg("water", []() { return new WaterTransport(); });
    addAlias("water", "Water");
    reg("high-pressure", []() { return new HighPressureGasTransport(); });
    addAlias("high-pressure", "HighP");
    m_CK_mode["CK_Mix"] = m_CK_mode["mixture-averaged-CK"] = true;
    m_CK_mode["CK_Multi"] = m_CK_mode["multicomponent-CK"] = true;
}

void TransportFactory::deleteFactory()
{
    std::unique_lock<std::mutex> transportLock(transport_mutex);
    delete s_factory;
    s_factory = 0;
}

Transport* TransportFactory::newTransport(const std::string& transportModel,
        ThermoPhase* phase, int log_level)
{
    if (transportModel != "DustyGas" && canonicalize(transportModel) == "none") {
        return create(transportModel);
    }
    if (!phase) {
        throw CanteraError("TransportFactory::newTransport",
            "Valid phase definition required for initialization of "
            "new '{}' object", transportModel);
    }

    vector_fp state;
    Transport* tr = 0;
    phase->saveState(state);

    if (transportModel == "DustyGas") {
        tr = new DustyGasTransport;
        Transport* gastr = new MultiTransport;
        gastr->init(phase, 0, log_level);
        DustyGasTransport* dtr = (DustyGasTransport*)tr;
        dtr->initialize(phase, gastr);
    } else {
        tr = create(transportModel);
        int mode = m_CK_mode[transportModel] ? CK_Mode : 0;
        tr->init(phase, mode, log_level);
    }
    phase->restoreState(state);
    return tr;
}

Transport* TransportFactory::newTransport(ThermoPhase* phase, int log_level)
{
    std::string transportModel = "None";
    XML_Node& phaseNode = phase->xml();
    AnyMap& input = phase->input();
    if (input.hasKey("transport")) {
        transportModel = input["transport"].asString();
    } else if (phaseNode.hasChild("transport")) {
        transportModel = phaseNode.child("transport").attrib("model");
    }
    return newTransport(transportModel, phase,log_level);
}

Transport* newTransportMgr(const std::string& model, ThermoPhase* thermo, int log_level)
{
    TransportFactory* f = TransportFactory::factory();
    return f->newTransport(model, thermo, log_level);
}

Transport* newDefaultTransportMgr(ThermoPhase* thermo, int loglevel)
{
    return TransportFactory::factory()->newTransport(thermo, loglevel);
}

}
