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
    case SRI3_FALLOFF:
        f = new SRI3();
        break;
    case SRI5_FALLOFF:
        f = new SRI5();
        break;
    default:
        return 0;
    }
    f->init(c);
    return f;
}

}
