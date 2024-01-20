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

    size_t componentIndex(const string& nm) const override;

    string componentName(size_t k) override;

    void getState(double* y) override;

    void initialize(double t0=0.0) override;

    void eval(double t, double* LHS, double* RHS) override;

    void updateState(double* y) override;

protected:
    const size_t m_sidx = 1;
};

}

#endif
