#ifndef CT_KERNEL_EQUIL_H
#define CT_KERNEL_EQUIL_H

//#include "ChemEquil.h"
#include "MultiPhase.h"

namespace Cantera {

    //-----------------------------------------------------------
    //              convenience functions
    //-----------------------------------------------------------

    /**
     * Return variable is equal to the number of subroutine attempts
     * it took to equilibrate the system.
     */
    int equilibrate(thermo_t& s, const char* XY,
        int solver = -1, doublereal rtol = 1.0e-9, int maxsteps = 1000, 
        int maxiter = 100, int loglevel = -99);

    doublereal equilibrate(MultiPhase& s, const char* XY,
        doublereal tol = 1.0e-9, int maxsteps = 1000, int maxiter = 100,
        int loglevel = -99);

}

#endif
