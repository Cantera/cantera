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

#include <iostream>
using namespace std;


namespace Cantera
{

//! default constructor
WaterTransport::WaterTransport(thermo_t* thermo, int ndim) :
    Transport(thermo, ndim)
{
    initTP();
}

//  Copy Constructor for the %WaterThermo object.
/*
 *    @param right  ThermoPhase to be copied
 */
WaterTransport::WaterTransport(const WaterTransport& right) :
    Transport(right.m_thermo, right.m_nDim)
{
    *this = right;
}

// Assignment operator
/*
 *
 * @param right    Reference to %WaterTransport object to be copied into the
 *                 current one.
 */
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

// Duplication routine for objects which inherit from %Transport
/*
 *  This virtual routine can be used to duplicate %Transport objects
 *  inherited from %Transport even if the application only has
 *  a pointer to %Transport to work with.
 *
 *  These routines are basically wrappers around the derived copy
 *  constructor.
 */
Transport*   WaterTransport::duplMyselfAsTransport() const
{
    return new WaterTransport(*this);
}


// virtual destructor
WaterTransport::~WaterTransport()
{
}

// Routine to do some common initializations at the start of using
// this routine.
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

// Returns the viscosity of water at the current conditions
// (kg/m/s)
/*
 *  This function calculates the value of the viscosity of pure
 *  water at the current T and P.
 *
 *  The formulas used are from the paper
 *     J. V. Sengers, J. T. R. Watson, "Improved International
 *     Formulations for the Viscosity and Thermal Conductivity of
 *     Water Substance", J. Phys. Chem. Ref. Data, 15, 1291 (1986).
 *
 *  The formulation is accurate for all temperatures and pressures,
 *  for steam and for water, even near the critical point.
 *  Pressures above 500 MPa and temperature above 900 C are suspect.
 */
doublereal WaterTransport::viscosity()
{
    return m_waterProps->viscosityWater();
}

// Returns the thermal conductivity of water at the current conditions
// (W/m/K)
/*
 *  This function calculates the value of the thermal conductivity of
 *  water at the current T and P.
 *
 *  The formulas used are from the paper
 *     J. V. Sengers, J. T. R. Watson, "Improved International
 *     Formulations for the Viscosity and Thermal Conductivity of
 *     Water Substance", J. Phys. Chem. Ref. Data, 15, 1291 (1986).
 *
 *  The formulation is accurate for all temperatures and pressures,
 *  for steam and for water, even near the critical point.
 *  Pressures above 500 MPa and temperature above 900 C are suspect.
 */
doublereal WaterTransport::thermalConductivity()
{
    return m_waterProps->thermalConductivityWater();
}

}
