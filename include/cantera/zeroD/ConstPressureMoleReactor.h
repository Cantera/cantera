//! @file ConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTPRESSMOLE_REACTOR_H
#define CT_CONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/**
 * ConstPressureMoleReactor is a class for constant-pressure reactors
 * which use a state of moles.
 * @since New in %Cantera 3.0
 * @ingroup reactorGroup
 */
class ConstPressureMoleReactor : public MoleReactor
{
public:
    using MoleReactor::MoleReactor; // inherit constructors

    string type() const override {
        return "ConstPressureMoleReactor";
    };

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void eval(double t, double* LHS, double* RHS) override;

    vector<size_t> steadyConstraints() const override {
        throw CanteraError("ConstPressureMoleReactor::steadyConstraints",
            "ConstPressureMoleReactor is not currently compatible with the steady-state"
            " solver.\nSee https://github.com/Cantera/enhancements/issues/234");
    }

    void updateState(double* y) override;

    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    double upperBound(size_t k) const override;
    double lowerBound(size_t k) const override;
    void resetBadValues(double* y) override;

protected:
    const size_t m_sidx = 1;
};

}

#endif
