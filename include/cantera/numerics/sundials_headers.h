#ifndef CT_SUNDIALS_HEADERS
#define CT_SUNDIALS_HEADERS

#include "sundials/sundials_types.h"
#include "sundials/sundials_math.h"
#include "sundials/sundials_nvector.h"
#include "nvector/nvector_serial.h"
#include "cvodes/cvodes.h"
#include "idas/idas.h"
#if CT_SUNDIALS_VERSION >= 30
    #if CT_SUNDIALS_USE_LAPACK
        #include "sunlinsol/sunlinsol_lapackdense.h"
        #include "sunlinsol/sunlinsol_lapackband.h"
    #else
        #include "sunlinsol/sunlinsol_dense.h"
        #include "sunlinsol/sunlinsol_band.h"
    #endif
    #include "sunlinsol/sunlinsol_spgmr.h"
    #include "cvodes/cvodes_direct.h"
    #include "cvodes/cvodes_diag.h"
    #include "cvodes/cvodes_spils.h"
    #include "idas/idas_direct.h"
    #include "idas/idas_spils.h"
#else
    #if CT_SUNDIALS_USE_LAPACK
        #include "cvodes/cvodes_lapack.h"
        #include "idas/idas_lapack.h"
    #else
        #include "cvodes/cvodes_dense.h"
        #include "cvodes/cvodes_band.h"
        #include "idas/idas_dense.h"
        #include "idas/idas_band.h"
    #endif
    #include "cvodes/cvodes_diag.h"
    #include "cvodes/cvodes_spgmr.h"
    #include "idas/idas_spgmr.h"
#endif

#define CV_SS 1
#define IDA_SS 1
#define CV_SV 2
#define IDA_SV 2

#if CT_SUNDIALS_VERSION < 25
typedef int sd_size_t;
#else
typedef long int sd_size_t;
#endif

#endif
