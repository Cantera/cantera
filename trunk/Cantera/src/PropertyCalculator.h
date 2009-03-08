/**
 *  @file PropertyCalculator.h 
 *
 * $Author: dggoodwin $
 * $Revision: 1.3 $
 * $Date: 2006/07/11 15:34:51 $
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_PROP_CALC_H
#define CT_PROP_CALC_H

#include "ct_defs.h"

namespace Cantera {

    /// Classes used by ChemEquil. These classes are used only by the
    /// ChemEquil equilibrium solver. Each one returns a particular
    /// property of the object supplied as the argument.
    ///
    template<class M>
    class PropertyCalculator {
    public:
        virtual ~PropertyCalculator(){}
        virtual doublereal value(const M& s) =0;
        virtual string symbol() =0;
    };

    template<class M>
    class EnthalpyCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.enthalpy_mass();
        }
        virtual string symbol() { return "H"; }
    };

    template<class M>
    class EntropyCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.entropy_mass();
        }
        virtual string symbol() { return "S"; }
    };

    template<class M>    
    class TemperatureCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.temperature();
        }
        virtual string symbol() { return "T"; }
    };

    template<class M>    
    class PressureCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.pressure();
        }
        virtual string symbol() { return "P"; }
    };

    template<class M>    
    class DensityCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.density();
        }
        virtual string symbol() { return "V"; }
    };

    template<class M>    
    class IntEnergyCalculator : public PropertyCalculator<M> {
    public:
        virtual doublereal value(const M& s) {
            return s.intEnergy_mass();
        }
        virtual string symbol() { return "U"; }
    };
}

#endif

