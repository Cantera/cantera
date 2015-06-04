/**
 *  @file WaterTransport.h Header file defining class WaterTransport
 */
#ifndef CT_WATERTRAN_H
#define CT_WATERTRAN_H

#include "LiquidTransportParams.h"
#include "cantera/thermo/WaterPropsIAPWS.h"

namespace Cantera
{
//! @{
const int LVISC_CONSTANT     = 0;
const int LVISC_WILKES       = 1;
const int LVISC_MIXTUREAVG   = 2;

const int LDIFF_MIXDIFF_UNCORRECTED     = 0;
const int LDIFF_MIXDIFF_FLUXCORRECTED  = 1;
const int LDIFF_MULTICOMP_STEFANMAXWELL  = 2;
//! @}


class WaterProps;
class PDSS_Water;

//! Transport Parameters for pure water
//! @ingroup tranprops
class WaterTransport : public Transport
{
public:
    //! default constructor
    /*!
     *  @param thermo   ThermoPhase object that represents the phase.
     *                  Defaults to zero
     *
     *  @param ndim     Number of dimensions of the flux expressions.
     *                  Defaults to a value of one.
     */
    WaterTransport(thermo_t* thermo = 0, int ndim = 1);

    WaterTransport(const WaterTransport& right);
    WaterTransport&  operator=(const  WaterTransport& right);
    virtual Transport* duplMyselfAsTransport() const;

    virtual int model() const {
        return cWaterTransport;
    }

    //! Returns the viscosity of water at the current conditions
    //! (kg/m/s)
    /*!
     *  This function calculates the value of the viscosity of pure
     *  water at the current T and P.
     *
     *  The formulas used are from the paper
     *
     *     J. V. Sengers, J. T. R. Watson, "Improved International
     *     Formulations for the Viscosity and Thermal Conductivity of
     *     Water Substance", J. Phys. Chem. Ref. Data, 15, 1291 (1986).
     *
     *  The formulation is accurate for all temperatures and pressures,
     *  for steam and for water, even near the critical point.
     *  Pressures above 500 MPa and temperature above 900 C are suspect.
     */
    virtual doublereal viscosity();

    virtual doublereal bulkViscosity() {
        return 0.0;
    }

    //! Returns the thermal conductivity of water at the current conditions
    //! (W/m/K)
    /*!
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
    virtual doublereal thermalConductivity();

private:

    //! Routine to do some common initializations at the start of using
    //! this routine.
    void initTP();

    //! Pointer to the WaterPropsIAPWS object, which does the actual calculations
    //! for the real equation of state
    /*!
     * This object owns m_sub
     */
    mutable WaterPropsIAPWS* m_sub;

    //! Pointer to the WaterProps object
    /*!
     *   This class is used to house several approximation
     *   routines for properties of water.
     *
     * This object owns m_waterProps, and the WaterPropsIAPWS object used by
     * WaterProps is m_sub, which is defined above.
     */
    WaterProps* m_waterProps;

    //! Pressure dependent standard state object for water
    /*!
     *  We assume that species 0 is water, with a PDSS_Water object.
     */
    PDSS_Water* m_waterPDSS;
};
}
#endif
