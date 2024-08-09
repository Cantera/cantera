//! @file StFlow.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STFLOW_H
#define CT_STFLOW_H

#include "Flow1D.h"

namespace Cantera
{

/**
 *  This class represents 1D flow domains that satisfy the one-dimensional
 *  similarity solution for chemically-reacting, axisymmetric flows.
 *
 *  @deprecated To be removed after %Cantera 3.1; replaced by Flow1D.
 *  @ingroup flowGroup
 */
class StFlow : public Flow1D
{
public:
    //! Create a new flow domain.
    //! @param ph Object representing the gas phase. This object will be used
    //!     to evaluate all thermodynamic, kinetic, and transport properties.
    //! @param nsp Number of species.
    //! @param points Initial number of grid points
    StFlow(ThermoPhase* ph = 0, size_t nsp = 1, size_t points = 1);

    //! Delegating constructor
    StFlow(shared_ptr<ThermoPhase> th, size_t nsp = 1, size_t points = 1);

    //! Create a new flow domain.
    //! @param sol  Solution object used to evaluate all thermodynamic, kinetic, and
    //!     transport properties
    //! @param id  name of flow domain
    //! @param points  initial number of grid points
    StFlow(shared_ptr<Solution> sol, const string& id="", size_t points=1);

    void eval(size_t j, double* x, double* r, integer* mask, double rdt) override;

    //! Evaluate all residual components at the right boundary.
    virtual void evalRightBoundary(double* x, double* res, int* diag, double rdt);

    void evalContinuity(size_t j, double* x, double* r, int* diag, double rdt) override;

protected:
    double wdot(size_t k, size_t j) const {
        return m_wdot(k,j);
    }

    //! Write the net production rates at point `j` into array `m_wdot`
    void getWdot(double* x, size_t j) {
        setGas(x,j);
        m_kin->getNetProductionRates(&m_wdot(0,j));
    }

    //! Evaluate the residual function. This function is called in eval
    //! after updateProperties is called.
    virtual void evalResidual(double* x, double* rsd, int* diag,
                              double rdt, size_t jmin, size_t jmax);
};

}

#endif
