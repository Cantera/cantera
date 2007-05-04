
#include "FlowDevice.h"
#include "ReactorBase.h"
#include "Func1.h"

using namespace std;

namespace CanteraZeroD {

    bool FlowDevice::install(ReactorBase& in, ReactorBase& out) {
        if (m_in || m_out) return false;
        m_in =  &in;
        m_out = &out;
        m_in->addOutlet(*this);
        m_out->addInlet(*this);
        
        // construct adapters between inlet and outlet species
        phase_t* mixin = &m_in->contents();
        phase_t* mixout = &m_out->contents();
        if (mixin == 0 || mixout == 0) return false;
        
        m_nspin = mixin->nSpecies();
        m_nspout = mixout->nSpecies();
        string nm;
        int ki, ko;
        for (ki = 0; ki < m_nspin; ki++) {
            nm = mixin->speciesName(ki);
            ko = mixout->speciesIndex(nm);
            m_in2out.push_back(ko);
        }
        for (ko = 0; ko < m_nspout; ko++) {
            nm = mixout->speciesName(ko);
            ki = mixin->speciesIndex(nm);
            m_out2in.push_back(ki);
        }
        return true; 
    }

    void FlowDevice::setFunction(Func1* f) {
        m_func = f;
    }


    /**
     * Mass flow rate of outlet species k.  Returns zero if this
     * species is not present in the upstream mixture.
     */
    doublereal FlowDevice::outletSpeciesMassFlowRate(int k) {
        if (k < 0 || k >= m_nspout) return 0.0;
        int ki = m_out2in[k];
        if (ki < 0) return 0.0;
        return m_mdot * m_in->massFraction(ki);
    }

    doublereal FlowDevice::enthalpy_mass() { return m_in->enthalpy_mass(); }

}        
