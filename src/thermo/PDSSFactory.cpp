//! @file PDSSFactory.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

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
    m_synonyms["constant_incompressible"] = "constant-incompressible";
    reg("water", []() { return new PDSS_Water(); });
    m_synonyms["waterPDSS"] = m_synonyms["waterIAPWS"] = "water";
    reg("ions-from-neutral", []() { return new PDSS_IonsFromNeutral(); });
    m_synonyms["IonFromNeutral"] = "ions-from-neutral";
    reg("temperature_polynomial", []() { return new PDSS_SSVol(); });
    m_synonyms["temperature-polynomial"] = "temperature_polynomial";
    m_synonyms["density_temperature_polynomial"] = "temperature_polynomial";
    m_synonyms["density-temperature-polynomial"] = "temperature_polynomial";
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
