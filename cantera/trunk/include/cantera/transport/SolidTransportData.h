/**
 *  @file SolidTransportData.h
 *  Header file defining class SolidTransportData
 */

#ifndef CT_SOLIDTRANSPORTDATA_H
#define CT_SOLIDTRANSPORTDATA_H

#include "cantera/transport/TransportBase.h"
#include "cantera/transport/TransportParams.h"
#include "cantera/base/FactoryBase.h"
#include "cantera/transport/LTPspecies.h"

namespace Cantera
{

//! Class SolidTransportData holds transport parameters for a
//! specific solid-phase species.
/*!
 * A SolidTransportData object is created for a solid phase
 * (not for each species as happens for the analogous LiquidTransportData).
 *
 * This class is mainly used to collect transport properties from the parse
 * phase in the TranportFactory and transfer them to the Transport class.
 * Transport properties are expressed by subclasses of LTPspecies. Note that
 * we use the liquid phase species model for the solid phases. That is, for
 * the time being at least, we ignore mixing models for solid phases and just
 * specify a transport property at the level that we specify the transport
 * property for a species in the liquid phase. One may need to be careful
 * about deleting pointers to LTPspecies objects created in the
 * TransportFactory.
 *
 * All of the pointers in this class are shallow pointers. Therefore, this is
 * a passthrough class, which keeps track of pointer ownership by zeroing
 * pointers as we go. Yes, Yes, yes, this is not good.
 */
class SolidTransportData : public TransportParams
{
public:
    SolidTransportData();
    SolidTransportData(const SolidTransportData& right);
    SolidTransportData& operator=(const SolidTransportData& right);
    ~SolidTransportData();

    //! A SolidTransportData object is instantiated for each species.
    //! This is the species name for which this object is instantiated.
    std::string speciesName;

    //! Model type for the ionic conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* ionConductivity;

    //! Model type for the thermal conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* thermalConductivity;

    //! Model type for the electrical conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* electConductivity;

    //! Model type for the defectDiffusivity -- or more like a defect diffusivity in the context of the solid phase.
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* defectDiffusivity;

    //! Model type for the defectActivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* defectActivity;

protected:
    //protected members of SolidTransportData are analogous to those found in TransportParams

    //! Local storage of the number of species
    //    int nsp_;

    //!  Pointer to the ThermoPhase object
    //    thermo_t* thermo;

    //! Local storage of the molecular weights of the species
    /*!
     *  Length is nsp_ and units are kg kmol-1.
     */
    //    vector_fp mw;

    //! Maximum temperatures for parameter fits
    //    doublereal tmax;

    //! Minimum temperatures for parameter fits
    //    doublereal tmin;

    //! Pointer to the xml tree describing the implementation of transport for this object
    //    XML_Writer* xml;

    //! Log level
    //    int log_level;

};

}
#endif
