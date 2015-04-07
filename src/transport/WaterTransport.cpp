//! @file WaterTransport.cpp
#include "cantera/base/ct_defs.h"
#include "cantera/thermo/WaterPropsIAPWS.h"
#include "cantera/transport/TransportBase.h"
#include "cantera/numerics/DenseMatrix.h"
#include "cantera/transport/LiquidTransportParams.h"
#include "cantera/thermo/VPStandardStateTP.h"

#include "cantera/transport/WaterTransport.h"
#include "cantera/thermo/PDSS_Water.h"
#include "cantera/thermo/WaterSSTP.h"
#include "cantera/thermo/WaterProps.h"

using namespace std;

namespace Cantera
{

WaterTransport::WaterTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim)
{
    initTP();
}

WaterTransport::WaterTransport(const WaterTransport& right) :
    Transport(right.m_thermo, right.m_nDim)
{
    *this = right;
}

WaterTransport&  WaterTransport::operator=(const  WaterTransport& right)
{
    if (&right != this) {
        return *this;
    }
    Transport::operator=(right);

    // All pointers in this routine are shallow pointers. Therefore, it's
    // ok just to reinitialize them
    initTP();

    return *this;
}

Transport*   WaterTransport::duplMyselfAsTransport() const
{
    return new WaterTransport(*this);
}

void WaterTransport::initTP()
{
    // The expectation is that we have a VPStandardStateTP derived object
    VPStandardStateTP* vpthermo = dynamic_cast<VPStandardStateTP*>(m_thermo);
    if (!vpthermo) {

        WaterSSTP* wsstp = dynamic_cast<WaterSSTP*>(m_thermo);
        if (!wsstp) {
            throw CanteraError("WaterTransport::initTP()",
                               "Expectation is that ThermoPhase be a VPStandardStateTP");
        } else {

            m_sub = wsstp->getWater();
            AssertTrace(m_sub != 0);
            // Get a pointer to a changeable WaterProps object
            m_waterProps = wsstp->getWaterProps();
            AssertTrace(m_waterProps != 0);
        }
    } else {
        m_waterPDSS = dynamic_cast<PDSS_Water*>(vpthermo->providePDSS(0));
        if (!m_waterPDSS) {
            throw CanteraError("WaterTransport::initTP()",
                               "Expectation is that first species be water with a PDSS_Water object");
        }
        // Get a pointer to a changeable WaterPropsIAPWS object
        m_sub = m_waterPDSS->getWater();
        AssertTrace(m_sub != 0);
        // Get a pointer to a changeable WaterProps object
        m_waterProps = m_waterPDSS->getWaterProps();
        AssertTrace(m_waterProps != 0);
    }
}

doublereal WaterTransport::viscosity()
{
    return m_waterProps->viscosityWater();
}

doublereal WaterTransport::thermalConductivity()
{
    return m_waterProps->thermalConductivityWater();
}

}
