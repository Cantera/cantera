//! @file ConstPressureReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_CONSTP_REACTOR_H
#define CT_CONSTP_REACTOR_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class ConstPressureReactor is a class for constant-pressure reactors. The
 * reactor may have an arbitrary number of inlets and outlets, each of which
 * may be connected to a "flow device" such as a mass flow controller, a
 * pressure regulator, etc. Additional reactors may be connected to the other
 * end of the flow device, allowing construction of arbitrary reactor
 * networks.
 * @ingroup reactorGroup
 */
class ConstPressureReactor : public Reactor
{
public:
    ConstPressureReactor(shared_ptr<Solution> sol, const string& name="(none)");
    ConstPressureReactor(shared_ptr<Solution> sol, bool clone,
                         const string& name="(none)");

    string type() const override {
        return "ConstPressureReactor";
    }

    void getState(span<double> y) override;
    void eval(double t, span<double> LHS, span<double> RHS) override;
    void evalSteady(double t, span<double> LHS, span<double> RHS) override;
    vector<size_t> initializeSteady() override;

    void updateState(span<const double> y) override;

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass", "enthalpy",
    //! the name of a homogeneous phase species, or the name of a surface
    //! species.
    size_t componentIndex(const string& nm) const override;
    string componentName(size_t k) override;
    double upperBound(size_t k) const override;
    double lowerBound(size_t k) const override;
    void resetBadValues(span<double> y) override;

protected:
    double m_initialMass; //!< Initial mass [kg]; used for steady-state calculations
};

}

#endif
