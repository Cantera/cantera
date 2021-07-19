//! @file IdealGasConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASCONSTPRESSMOLE_REACTOR_H
#define CT_IDEALGASCONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/*!
 * IdealGasConstPressureMoleReactor is a class for ideal gas constant-pressure reactors
 * which use a state of moles.
 */
class IdealGasConstPressureMoleReactor : public MoleReactor
{
public:
    IdealGasConstPressureMoleReactor() {}

    virtual std::string type() const {
        return "IdealGasConstPressureMoleReactor";
    };

    virtual size_t componentIndex(const std::string& nm) const;

    virtual std::string componentName(size_t k);

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(double* y);

    virtual void initialize(double t0 = 0.0);

    virtual void eval(double t, double* LHS, double* RHS);

    virtual void updateState(double* y);

    virtual Eigen::SparseMatrix<double> jacobian(double t, double* y);

protected:
    vector_fp m_hk; //!< Species molar enthalpies

    const int m_sidx = 1;
};

}

#endif
