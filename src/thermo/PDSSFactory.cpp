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
    addAlias("constant-volume", "constant_incompressible");
    addAlias("constant-volume", "constant-incompressible");
    reg("liquid-water-IAPWS95", []() { return new PDSS_Water(); });
    addAlias("liquid-water-IAPWS95", "waterPDSS");
    addAlias("liquid-water-IAPWS95", "waterIAPWS");
    addAlias("liquid-water-IAPWS95", "water");
    reg("ions-from-neutral-molecule", []() { return new PDSS_IonsFromNeutral(); });
    addAlias("ions-from-neutral-molecule", "IonFromNeutral");
    addAlias("ions-from-neutral-molecule", "ions-from-neutral");
    reg("molar-volume-temperature-polynomial", []() { return new PDSS_SSVol(); });
    addAlias("molar-volume-temperature-polynomial", "temperature_polynomial");
    reg("density-temperature-polynomial", []() { return new PDSS_SSVol(); });
    addAlias("density-temperature-polynomial", "density_temperature_polynomial");
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
