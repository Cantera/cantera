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
#include "cantera/base/stringUtils.h"
#include "cantera/base/utilities.h"

namespace Cantera
{
TransportFactory* TransportFactory::s_factory = 0;

// declaration of static storage for the mutex
std::mutex TransportFactory::transport_mutex;

//////////////////// class TransportFactory methods //////////////

TransportFactory::TransportFactory()
{
    reg("none", []() { return new Transport(); });
    addDeprecatedAlias("none", "Transport");
    addDeprecatedAlias("none", "None");
    addDeprecatedAlias("none", "");
    reg("unity-Lewis-number", []() { return new UnityLewisTransport(); });
    addDeprecatedAlias("unity-Lewis-number", "UnityLewis");
    reg("mixture-averaged", []() { return new MixTransport(); });
    addDeprecatedAlias("mixture-averaged", "Mix");
    reg("mixture-averaged-CK", []() { return new MixTransport(); });
    addDeprecatedAlias("mixture-averaged-CK", "CK_Mix");
    reg("multicomponent", []() { return new MultiTransport(); });
    addDeprecatedAlias("multicomponent", "Multi");
    reg("multicomponent-CK", []() { return new MultiTransport(); });
    addDeprecatedAlias("multicomponent-CK", "CK_Multi");
    reg("ionized-gas", []() { return new IonGasTransport(); });
    addDeprecatedAlias("ionized-gas", "Ion");
    reg("water", []() { return new WaterTransport(); });
    addDeprecatedAlias("water", "Water");
    reg("high-pressure", []() { return new HighPressureGasTransport(); });
    addDeprecatedAlias("high-pressure", "HighP");
    m_CK_mode["CK_Mix"] = m_CK_mode["mixture-averaged-CK"] = true;
    m_CK_mode["CK_Multi"] = m_CK_mode["multicomponent-CK"] = true;
}

TransportFactory* TransportFactory::factory() {
    std::unique_lock<std::mutex> transportLock(transport_mutex);
    if (!s_factory) {
        s_factory = new TransportFactory();
    }
    return s_factory;
}

void TransportFactory::deleteFactory()
{
    std::unique_lock<std::mutex> transportLock(transport_mutex);
    delete s_factory;
    s_factory = 0;
}

Transport* TransportFactory::newTransport(const string& transportModel,
        ThermoPhase* phase, int log_level)
{
    if (log_level != -7) {
        warn_deprecated("TransportFactory::newTransport", "The log_level parameter "
            "is deprecated and will be removed after Cantera 3.1.");
    }
    if (transportModel != "DustyGas" && canonicalize(transportModel) == "none") {
        return create(transportModel);
    }
    if (!phase) {
        throw CanteraError("TransportFactory::newTransport",
            "Valid phase definition required for initialization of "
            "new '{}' object", transportModel);
    }

    vector<double> state;
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
    if (log_level != -7) {
        warn_deprecated("TransportFactory::newTransport", "The log_level parameter "
            "is deprecated and will be removed after Cantera 3.1.");
    }
    string transportModel = "none";
    AnyMap& input = phase->input();
    if (input.hasKey("transport")) {
        transportModel = input["transport"].asString();
    }
    return newTransport(transportModel, phase,log_level);
}

shared_ptr<Transport> newTransport(shared_ptr<ThermoPhase> thermo, const string& model)
{
    Transport* tr;
    if (model == "default") {
        tr = TransportFactory::factory()->newTransport(thermo.get());
    } else {
        tr = TransportFactory::factory()->newTransport(model, thermo.get());
    }
    return shared_ptr<Transport>(tr);
}

}
