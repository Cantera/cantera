//! @file ConstPressureMoleReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTPRESSMOLE_REACTOR_H
#define CT_CONSTPRESSMOLE_REACTOR_H

#include "cantera/zeroD/MoleReactor.h"

namespace Cantera
{

/*!
 * ConstPressureMoleReactor is a class for constant-pressure reactors
 * which use a state of moles.
 * @since New in Cantera 3.0
 */
class ConstPressureMoleReactor : public MoleReactor
{
public:
    ConstPressureMoleReactor() {}

    virtual std::string type() const {
        return "ConstPressureMoleReactor";
    };

    virtual size_t componentIndex(const std::string& nm) const;

    virtual std::string componentName(size_t k);

    virtual void getState(double* y);

    virtual void initialize(double t0 = 0.0);

    virtual void eval(double t, double* LHS, double* RHS);

    virtual void updateState(double* y);

protected:
    const size_t m_sidx = 1;
};

}

#endif
