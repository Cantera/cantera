/**
 *  @file FalloffFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/kinetics/reaction_defs.h"

namespace Cantera
{

FalloffFactory* FalloffFactory::s_factory = 0;
std::mutex FalloffFactory::falloff_mutex;

FalloffFactory::FalloffFactory()
{
    reg("Lindemann", []() { return new Falloff(); });
    addAlias("Lindemann", "Simple");
    reg("Troe", []() { return new Troe(); });
    reg("SRI", []() { return new SRI(); });
}

Falloff* FalloffFactory::newFalloff(int type, const vector_fp& c)
{
    warn_deprecated("FalloffFactory::newFalloff",
        "Instantiation using magic numbers is deprecated; use string "
        "identifier instead. To be removed after Cantera 2.5.");
    static const std::unordered_map<int, std::string> types {
        {SIMPLE_FALLOFF, "Simple"},
        {TROE_FALLOFF, "Troe"},
        {SRI_FALLOFF, "SRI"}
    };

    Falloff* f = create(types.at(type));
    f->init(c);
    return f;
}

Falloff* FalloffFactory::newFalloff(const std::string& type, const vector_fp& c)
{
    Falloff* f = create(type);
    f->init(c);
    return f;
}

shared_ptr<Falloff> newFalloff(int type, const vector_fp& c)
{
    shared_ptr<Falloff> f(FalloffFactory::factory()->newFalloff(type, c));
    return f;
}

shared_ptr<Falloff> newFalloff(const std::string& type, const vector_fp& c)
{
    shared_ptr<Falloff> f(FalloffFactory::factory()->newFalloff(type, c));
    return f;
}

}
