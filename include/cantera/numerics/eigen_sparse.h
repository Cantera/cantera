#include "cantera/base/ct_defs.h"
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#else
#include "cantera/ext/Eigen/Sparse"
#include "cantera/ext/Eigen/SparseQR"
#endif

namespace Cantera {

}
