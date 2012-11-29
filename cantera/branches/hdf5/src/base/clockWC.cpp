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

/*
 * Reinitialize the tick counters within the object
 */
double clockWC::start()
{
    start_ticks = last_num_ticks = clock();
    clock_rollovers = 0u;
    return 0.0;
}

/*
 *    Returns system cpu and wall clock time in seconds. This
 *    is a strictly Ansi C timer, since clock() is defined as an
 *    Ansi C function. On some machines clock() returns type
 *    unsigned long (HP) and on others (SUN) it returns type long.
 *       An attempt to recover the actual time for clocks which have
 *    rolled over is made also. However, it only works if this
 *    function is called fairly regularily during
 *    the solution procedure.
 *
 *    clock() -> returns the time in microseconds. Division by
 *               the macro CLOCKS_PER_SEC recovers the time in seconds.
 */
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
    return(value);
}
}
