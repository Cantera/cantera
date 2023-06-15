//! @file PDSSFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/thermo/PDSSFactory.h"
#include "cantera/thermo/PDSS_IdealGas.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/PDSS_ConstVol.h"
#include "cantera/thermo/PDSS_SSVol.h"
#include "cantera/thermo/PDSS_HKFT.h"
#include "cantera/thermo/PDSS_IonsFromNeutral.h"

namespace Cantera
{

PDSSFactory* PDSSFactory::s_factory = 0;
std::mutex PDSSFactory::thermo_mutex;

PDSSFactory::PDSSFactory()
{
    reg("ideal-gas", []() { return new PDSS_IdealGas(); });
    reg("constant-volume", []() { return new PDSS_ConstVol(); });
    addDeprecatedAlias("constant-volume", "constant_incompressible");
    addDeprecatedAlias("constant-volume", "constant-incompressible");
    reg("liquid-water-IAPWS95", []() { return new PDSS_Water(); });
    addDeprecatedAlias("liquid-water-IAPWS95", "waterPDSS");
    addDeprecatedAlias("liquid-water-IAPWS95", "waterIAPWS");
    addDeprecatedAlias("liquid-water-IAPWS95", "water");
    reg("ions-from-neutral-molecule", []() { return new PDSS_IonsFromNeutral(); });
    addDeprecatedAlias("ions-from-neutral-molecule", "IonFromNeutral");
    addDeprecatedAlias("ions-from-neutral-molecule", "ions-from-neutral");
    reg("molar-volume-temperature-polynomial", []() { return new PDSS_SSVol(); });
    addDeprecatedAlias("molar-volume-temperature-polynomial", "temperature_polynomial");
    reg("density-temperature-polynomial", []() { return new PDSS_SSVol(); });
    addDeprecatedAlias("density-temperature-polynomial", "density_temperature_polynomial");
    reg("HKFT", []() { return new PDSS_HKFT(); });
}

PDSSFactory* PDSSFactory::factory() {
    std::unique_lock<std::mutex> lock(thermo_mutex);
    if (!s_factory) {
        s_factory = new PDSSFactory;
    }
    return s_factory;
}

void PDSSFactory::deleteFactory() {
    std::unique_lock<std::mutex> lock(thermo_mutex);
    delete s_factory;
    s_factory = 0;
}

PDSS* PDSSFactory::newPDSS(const std::string& model)
{
    return create(model);
}

PDSS* newPDSS(const std::string& model)
{
    return PDSSFactory::factory()->newPDSS(model);
}

}
