/**
 *  @file FalloffFactory.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/kinetics/reaction_defs.h"

namespace Cantera
{

FalloffFactory* FalloffFactory::s_factory = 0;
std::mutex FalloffFactory::falloff_mutex;

FalloffFactory::FalloffFactory()
{
    reg("Simple", []() { return new Falloff(); });
    reg("Troe", []() { return new Troe(); });
    reg("SRI", []() { return new SRI(); });
}

Falloff* FalloffFactory::newFalloff(int type, const vector_fp& c)
{
    static const std::unordered_map<int, std::string> types {
        {SIMPLE_FALLOFF, "Simple"},
        {TROE_FALLOFF, "Troe"},
        {SRI_FALLOFF, "SRI"}
    };

    Falloff* f = create(types.at(type));
    f->init(c);
    return f;
}

shared_ptr<Falloff> newFalloff(int type, const vector_fp& c)
{
    shared_ptr<Falloff> f(FalloffFactory::factory()->newFalloff(type, c));
    return f;
}

}
