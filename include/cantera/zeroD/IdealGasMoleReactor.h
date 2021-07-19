//! @file IdealGasMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASMOLE_REACTOR_H
#define CT_IDEALGASMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/*!
 * IdealGasMoleReactor is a class for ideal gas constant-volume reactors which use a
 * state of moles.
 */
class IdealGasMoleReactor : public MoleReactor
{
public:
    IdealGasMoleReactor() {}

    virtual std::string type() const {
        return "IdealGasMoleReactor";
    }

    virtual size_t componentIndex(const std::string& nm) const;

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(double* y);

    virtual void initialize(double t0 = 0.0);

    virtual void eval(double t, double* LHS, double* RHS);

    virtual void updateState(double* y);

    virtual Eigen::SparseMatrix<double> jacobian(double t, double* y);

protected:
    vector_fp m_uk; //!< Species molar internal energies
};

}

#endif
