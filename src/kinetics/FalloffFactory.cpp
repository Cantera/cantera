/**
 *  @file FalloffFactory.cpp
 */
// Copyright 2001  California Institute of Technology

#include "cantera/kinetics/FalloffFactory.h"
#include "cantera/kinetics/reaction_defs.h"

namespace Cantera
{

FalloffFactory* FalloffFactory::s_factory = 0;
mutex_t FalloffFactory::falloff_mutex;

Falloff* FalloffFactory::newFalloff(int type, const vector_fp& c)
{
    Falloff* f;
    switch (type) {
    case SIMPLE_FALLOFF:
        f = new Falloff();
        break;
    case TROE_FALLOFF:
        f = new Troe();
        break;
    case SRI_FALLOFF:
        f = new SRI();
        break;
    default:
        return 0;
    }
    f->init(c);
    return f;
}

shared_ptr<Falloff> newFalloff(int type, const vector_fp& c)
{
    shared_ptr<Falloff> f(FalloffFactory::factory()->newFalloff(type, c));
    return f;
}

}
