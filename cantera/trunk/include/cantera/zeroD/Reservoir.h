/**
 * @file Reservoir.h
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_RESERVOIR_H
#define CT_RESERVOIR_H

#include "ReactorBase.h"

namespace Cantera
{

class Reservoir : public ReactorBase
{
public:
    Reservoir() {}
    virtual int type() const {
        return ReservoirType;
    }
    virtual void initialize(doublereal t0 = 0.0) {}

    void insert(Cantera::ThermoPhase& contents) {
        setThermoMgr(contents);
    }
};

}

#endif
