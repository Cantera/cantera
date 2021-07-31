#include "cantera/base/ct_defs.h"
#if CT_USE_SYSTEM_EIGEN
#include <Eigen/Dense>
#else
#include "cantera/ext/Eigen/Dense"
#endif

namespace Cantera {
    typedef Eigen::Map<Eigen::MatrixXd> MappedMatrix;
    typedef Eigen::Map<Eigen::VectorXd> MappedVector;
    typedef Eigen::Map<const Eigen::VectorXd> ConstMappedVector;
    typedef Eigen::Map<Eigen::RowVectorXd> MappedRowVector;
    typedef Eigen::Map<const Eigen::RowVectorXd> ConstMappedRowVector;
}
