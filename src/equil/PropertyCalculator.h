/**
 *  @file PropertyCalculator.h
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_PROP_CALC_H
#define CT_PROP_CALC_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

/// Classes used by ChemEquil. These classes are used only by the
/// ChemEquil equilibrium solver. Each one returns a particular
/// property of the object supplied as the argument.
///
template<class M>
class PropertyCalculator
{
public:
    virtual ~PropertyCalculator() {}
    virtual doublereal value(const M& s) =0;
    virtual std::string symbol() =0;
};

template<class M>
class EnthalpyCalculator : public PropertyCalculator<M>
{
public:
    virtual doublereal value(const M& s) {
        return s.enthalpy_mass();
    }
    virtual std::string symbol() {
        return "H";
    }
};

template<class M>
class EntropyCalculator : public PropertyCalculator<M>
{
public:
    virtual doublereal value(const M& s) {
        return s.entropy_mass();
    }
    virtual std::string symbol() {
        return "S";
    }
};

template<class M>
class TemperatureCalculator : public PropertyCalculator<M>
{
public:
    virtual doublereal value(const M& s) {
        return s.temperature();
    }
    virtual std::string symbol() {
        return "T";
    }
};

template<class M>
class PressureCalculator : public PropertyCalculator<M>
{
public:
    virtual doublereal value(const M& s) {
        return s.pressure();
    }
    virtual std::string symbol() {
        return "P";
    }
};

template<class M>
class DensityCalculator : public PropertyCalculator<M>
{
public:
    virtual doublereal value(const M& s) {
        return s.density();
    }
    virtual std::string symbol() {
        return "V";
    }
};

template<class M>
class IntEnergyCalculator : public PropertyCalculator<M>
{
public:
    virtual doublereal value(const M& s) {
        return s.intEnergy_mass();
    }
    virtual std::string symbol() {
        return "U";
    }
};
}

#endif

