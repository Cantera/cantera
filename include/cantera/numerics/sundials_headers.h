#ifndef CT_SUNDIALS_HEADERS
#define CT_SUNDIALS_HEADERS

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#include "idas/idas.h"

#if CT_SUNDIALS_USE_LAPACK
    #include "sunlinsol/sunlinsol_lapackdense.h"
    #include "sunlinsol/sunlinsol_lapackband.h"
#else
    #include "sunlinsol/sunlinsol_dense.h"
    #include "sunlinsol/sunlinsol_band.h"
#endif
#include "sunlinsol/sunlinsol_spgmr.h"
#include "cvodes/cvodes_diag.h"

#if SUNDIALS_VERSION_MAJOR < 7
    #include "cvodes/cvodes_direct.h"
    #include "idas/idas_direct.h"
    #include "idas/idas_spils.h"
    #include "cvodes/cvodes_spils.h"
#endif

#if SUNDIALS_VERSION_MAJOR < 6
    typedef realtype sunrealtype;
    typedef booleantype sunbooleantype;
#endif

#define CV_SS 1
#define IDA_SS 1
#define CV_SV 2
#define IDA_SV 2

typedef long int sd_size_t;

#endif
