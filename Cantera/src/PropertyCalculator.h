/**
 *  @file PropertyCalculator.h 
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_PROP_CALC_H
#define CT_PROP_CALC_H

#include "ct_defs.h"

namespace Cantera {

    template<class M>
    class PropertyCalculator {
    public:
        virtual doublereal value(const M& s) =0;
    };

    template<class M>
    class EnthalpyCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.enthalpy_mass();
        }
    };

    template<class M>
    class EntropyCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.entropy_mass();
        }
    };

    template<class M>    
    class TemperatureCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.temperature();
        }
    };

    template<class M>    
    class PressureCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.pressure();
        }
    };

    template<class M>    
    class DensityCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.density();
        }
    };

    template<class M>    
    class IntEnergyCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.intEnergy_mass();
        }
    };
}

#endif

