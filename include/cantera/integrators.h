/**
 * @file integrators.h
 * ODE integrators. Currently, the only integrator is CVODE.
 * @deprecated To be removed after Cantera 2.3. Include relevant headers directly.
 */
#ifndef CT_INTEG_H_INCL
#define CT_INTEG_H_INCL

#pragma message "Deprecated. integrators.h will be removed after Cantera 2.3. Include relevant headers directly."

#include "numerics/Integrator.h"
#include "numerics/DAE_Solver.h"
#include "numerics/IDA_Solver.h"

#endif
