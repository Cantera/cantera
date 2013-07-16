/**
 *  @file IdealGasReactor.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_IDEALGASREACTOR_H
#define CT_IDEALGASREACTOR_H

#include "Reactor.h"
#include "cantera/kinetics/Kinetics.h"

namespace Cantera
{

/**
 * Class IdealGasReactor is a class for stirred reactors that is specifically
 * optimized for ideal gases. In this formulation, temperature replaces the
 * total internal energy as a state variable.
 */
class IdealGasReactor : public Reactor
{
public:
    IdealGasReactor() {}

    virtual int type() const {
        return IdealGasReactorType;
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getInitialConditions(doublereal t0, size_t leny,
                                      doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    virtual void evalEqs(doublereal t, doublereal* y,
                         doublereal* ydot, doublereal* params);

    virtual void updateState(doublereal* y);

    virtual size_t componentIndex(const std::string& nm) const;

protected:
    vector_fp m_uk; //!< Species molar internal energies
};

}

#endif
