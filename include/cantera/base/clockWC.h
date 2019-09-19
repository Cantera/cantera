/**
 * @file clockWC.h
 *    Declarations for a simple class that implements an Ansi C wall clock timer
 *   (see \ref Cantera::clockWC).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CLOCKWC_H
#define CT_CLOCKWC_H

#include <time.h>
namespace Cantera
{

//! The class provides the wall clock timer in seconds
/*!
 * This routine relies on the ANSI C routine, clock(), for its basic operation.
 * Therefore, it should be fairly portable.
 *
 * The clock will rollover if the calculation is long enough. The wraparound
 * time is roughly 72 minutes for a 32 bit system. This object senses that by
 * seeing if the raw tick counter is has decreased from the last time. If it
 * senses a wraparound has occurred, it increments an internal counter to
 * account for this. Therefore, for long calculations, this object must be
 * called at regular intervals for the seconds timer to be accurate.
 *
 * An example of how to use the timer is given below. timeToDoCalcs contains the
 * wall clock time calculated for the operation.
 *
 * @code
 * clockWC wc;
 * do_hefty_calculations_atLeastgreaterThanAMillisecond();
 * double timeToDoCalcs = wc.secondsWC();
 * @endcode
 *
 * In general, the process to be timed must take more than a millisecond for
 * this clock to enough of a significant resolution to be accurate.
 *
 * @ingroup globalUtilFuncs
 *
 */
class clockWC
{
public:
    //! Constructor
    /*!
     * This also serves to initialize the ticks within the object
     */
    clockWC();

    //! Resets the internal counters and returns the wall clock time in seconds
    double start();

    //! Returns the wall clock time in seconds since the last reset.
    /*!
     * Returns system cpu and wall clock time in seconds. This is a strictly
     * Ansi C timer, since clock() is defined as an Ansi C function. On some
     * machines clock() returns type unsigned long (HP) and on others (SUN)
     * it returns type long. An attempt to recover the actual time for clocks
     * which have rolled over is made also. However, it only works if this
     * function is called fairly regularily during the solution procedure.
     */
    double secondsWC();

private:
    //! Counters the value of the number of ticks from the last call.
    clock_t last_num_ticks;

    //! Number of clock rollovers since the last initialization
    /*!
     * The clock will rollover if the calculation is long enough. This object
     * senses that by seeing if the raw tick counter is has decreased from the
     * last time.
     */
    unsigned int clock_rollovers;

    //! Counter containing the value of the number of ticks from the first call
    //! (or the reset call).
    clock_t start_ticks;

    //! internal constant containing clock ticks per second
    const double inv_clocks_per_sec;

    //! internal constant containing the total number of ticks per rollover.
    const double clock_width;
};
}
#endif
