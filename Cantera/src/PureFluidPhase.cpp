
//#include "ct_defs.h"
#include "xml.h"
#include "PureFluidPhase.h"

namespace Cantera {

    void PureFluid::
    initThermo() {
        if (m_sub) delete m_sub;
        m_sub = tpx::GetSub(m_subflag);
        if (m_sub == 0) {
            throw CanteraError("PureFluid::initThermo",
                "could not create new substance object.");
        }
        m_mw = m_sub->MolWt();
        m_weight[0] = m_mw;
        setMolecularWeight(0,m_mw);
        double one = 1.0;
        setMoleFractions(&one);
        double cp0_R, h0_RT, s0_R, T0, p;
        T0 = 298.15;
        if (T0 < m_sub->Tcrit()) {
            m_sub->Set(tpx::TX, T0, 1.0);
            p = 0.01*m_sub->P();
        }
        else {
            p = 0.001*m_sub->Pcrit();
        }
        m_sub->Set(tpx::TP, T0, p);
        
        m_spthermo->update_one(0, T0, &cp0_R, &h0_RT, &s0_R);
        double s_R = s0_R - log(p/refPressure());
        m_sub->setStdState(h0_RT*GasConstant*298.15/m_mw,
            s_R*GasConstant/m_mw, T0, p);
        if (m_verbose) {
            writelog("PureFluid::initThermo: initialized phase "
                +id()+"\n");
        }
    }

    void PureFluid::
    setParametersFromXML(const XML_Node& eosdata) {
        eosdata.require("model","PureFluid");
        m_subflag = atoi(eosdata["fluid_type"].c_str());
        if (m_subflag < 0) 
            throw CanteraError("PureFluid::setParametersFromXML",
                "missing or negative substance flag");
    }

}
