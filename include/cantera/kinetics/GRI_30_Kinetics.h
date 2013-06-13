/**
 * @file GRI_30_Kinetics.h
 *
 * @ingroup chemkinetics
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_GRI30_KINETICS_H
#define CT_GRI30_KINETICS_H

#include "GasKinetics.h"

namespace Cantera
{
const int cGRI_30_Kinetics = cGasKinetics + 1;

/**
 *  Kinetics manager implementing reaction mechanism GRI-Mech 3.0
 *  @deprecated
 */
class GRI_30_Kinetics : public GasKinetics
{
public:
    /// Default constructor.
    GRI_30_Kinetics(thermo_t* th=0);

    virtual int type() const {
        return cGRI_30_Kinetics;
    }

    virtual void getNetProductionRates(doublereal* net) {
        gri30_updateROP();
        get_wdot(&m_ropnet[0], net);
    }

private:
    void gri30_update_rates_T();
    void gri30_updateROP();

    /**
     * Update the equilibrium constants in molar units.
    */
    void gri30_updateKc();

    void get_wdot(const doublereal* rop, doublereal* wdot);
    void update_kc(const doublereal* grt, doublereal  c0, doublereal* rkc);
    void update_rates(doublereal t, doublereal  tlog, doublereal* rf);
    void eval_ropnet(const doublereal* c, const doublereal* rf, const doublereal* rkc, doublereal* r);
};
}

#endif
