#include "cantera/base/ct_defs.h"
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>
#else
#include "cantera/ext/Eigen/Sparse"
#include "cantera/ext/Eigen/SparseQR"
#include "cantera/ext/Eigen/SparseLU"
#include<cantera/ext/Eigen/IterativeLinearSolvers>
#endif

namespace Cantera {

}
