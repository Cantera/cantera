//! @file FlowDevice.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/ReactorBase.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/numerics/Func1.h"

namespace Cantera
{

bool FlowDevice::install(ReactorBase& in, ReactorBase& out)
{
    if (m_in || m_out) {
        throw CanteraError("FlowDevice::install", "Already installed");
    }
    m_in =  &in;
    m_out = &out;
    m_in->addOutlet(*this);
    m_out->addInlet(*this);

    // construct adapters between inlet and outlet species
    const ThermoPhase& mixin = m_in->contents();
    const ThermoPhase& mixout = m_out->contents();

    m_nspin = mixin.nSpecies();
    m_nspout = mixout.nSpecies();
    string nm;
    size_t ki, ko;
    for (ki = 0; ki < m_nspin; ki++) {
        nm = mixin.speciesName(ki);
        ko = mixout.speciesIndex(nm);
        m_in2out.push_back(ko);
    }
    for (ko = 0; ko < m_nspout; ko++) {
        nm = mixout.speciesName(ko);
        ki = mixin.speciesIndex(nm);
        m_out2in.push_back(ki);
    }
    return true;
}

void FlowDevice::setPressureFunction(Func1* f)
{
    m_pfunc = f;
}

double FlowDevice::evalPressureFunction()
{
    double delta_P = in().pressure() - out().pressure();
    if (m_pfunc) {
        return m_pfunc->eval(delta_P);
    }
    return delta_P;
}

void FlowDevice::setTimeFunction(Func1* g)
{
    m_tfunc = g;
}

double FlowDevice::evalTimeFunction()
{
    if (m_tfunc) {
        return m_tfunc->eval(m_time);
    }
    return 1.;
}

double FlowDevice::outletSpeciesMassFlowRate(size_t k)
{
    if (k >= m_nspout) {
        return 0.0;
    }
    size_t ki = m_out2in[k];
    if (ki == npos) {
        return 0.0;
    }
    return m_mdot * m_in->massFraction(ki);
}

double FlowDevice::enthalpy_mass()
{
    return m_in->enthalpy_mass();
}

}
