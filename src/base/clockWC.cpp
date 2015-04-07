/**
 * @file clockWC.cpp
 *    Definitions for a simple class that implements an Ansi C wall clock timer
 *     (see \ref Cantera::clockWC).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include <time.h>
#include "cantera/base/clockWC.h"

namespace Cantera
{
clockWC::clockWC() :
    last_num_ticks(clock()),
    clock_rollovers(0u),
    start_ticks(0),
    inv_clocks_per_sec(1./(double)CLOCKS_PER_SEC),
    clock_width((double)(1L<<((int)sizeof(clock_t)*8-2))*4./(double)CLOCKS_PER_SEC)
{
    start_ticks = last_num_ticks;
}

double clockWC::start()
{
    start_ticks = last_num_ticks = clock();
    clock_rollovers = 0u;
    return 0.0;
}

double clockWC::secondsWC()
{
    clock_t num_ticks = clock();
    if (num_ticks < last_num_ticks) {
        clock_rollovers++;
    }
    double value = (num_ticks - start_ticks) * inv_clocks_per_sec;
    if (clock_rollovers) {
        value += clock_rollovers * clock_width;
    }
    last_num_ticks = num_ticks;
    return value;
}
}
