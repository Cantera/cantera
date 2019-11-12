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
    reg("constant-incompressible", []() { return new PDSS_ConstVol(); });
    addAlias("constant-incompressible", "constant_incompressible");
    addAlias("constant-incompressible", "constant-volume");
    reg("water", []() { return new PDSS_Water(); });
    addAlias("water", "waterPDSS");
    addAlias("water", "waterIAPWS");
    addAlias("water", "liquid-water-IAPWS95");
    reg("ions-from-neutral", []() { return new PDSS_IonsFromNeutral(); });
    addAlias("ions-from-neutral", "IonFromNeutral");
    addAlias("ions-from-neutral", "ions-from-neutral-molecule");
    reg("temperature_polynomial", []() { return new PDSS_SSVol(); });
    addAlias("temperature_polynomial", "molar-volume-temperature-polynomial");
    reg("density_temperature_polynomial", []() { return new PDSS_SSVol(); });
    addAlias("density_temperature_polynomial", "density-temperature-polynomial");
    reg("HKFT", []() { return new PDSS_HKFT(); });
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
