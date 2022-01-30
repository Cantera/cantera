/**
 *  @file FalloffFactory.cpp
 *
 *  @deprecated  Deprecated in Cantera 2.6 and removed thereafter. Replaced by
 *      FalloffRate objects managed by MultiRate evaluators.
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
    reg("Lindemann", []() { return new Lindemann(); });
    addAlias("Lindemann", "Simple");
    reg("Troe", []() { return new Troe(); });
    reg("SRI", []() { return new SRI(); });
    reg("Tsang", []() { return new Tsang(); });
}

Falloff* FalloffFactory::newFalloff(const std::string& type, const vector_fp& c)
{
    Falloff* f = create(type);
    f->init(c);
    return f;
}

shared_ptr<Falloff> newFalloff(const std::string& type, const vector_fp& c)
{
    shared_ptr<Falloff> f(FalloffFactory::factory()->newFalloff(type, c));
    return f;
}

}
