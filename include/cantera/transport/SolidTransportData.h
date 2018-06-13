/**
 *  @file SolidTransportData.h
 *  Header file defining class SolidTransportData
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_SOLIDTRANSPORTDATA_H
#define CT_SOLIDTRANSPORTDATA_H

#include "cantera/transport/TransportParams.h"
#include "cantera/transport/LTPspecies.h"

namespace Cantera
{

//! Class SolidTransportData holds transport parameters for a specific solid-
//! phase species.
/*!
 * @deprecated To be removed after Cantera 2.4
 *
 * A SolidTransportData object is created for a solid phase
 * (not for each species as happens for the analogous LiquidTransportData).
 *
 * This class is mainly used to collect transport properties from the parse
 * phase in the TranportFactory and transfer them to the Transport class.
 * Transport properties are expressed by subclasses of LTPspecies. Note that we
 * use the liquid phase species model for the solid phases. That is, for the
 * time being at least, we ignore mixing models for solid phases and just
 * specify a transport property at the level that we specify the transport
 * property for a species in the liquid phase. One may need to be careful about
 * deleting pointers to LTPspecies objects created in the TransportFactory.
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
    LTPspecies* ionConductivity;

    //! Model type for the thermal conductivity
    LTPspecies* thermalConductivity;

    //! Model type for the electrical conductivity
    LTPspecies* electConductivity;

    //! Model type for the defectDiffusivity -- or more like a defect diffusivity in the context of the solid phase.
    LTPspecies* defectDiffusivity;

    //! Model type for the defectActivity
    LTPspecies* defectActivity;
};

}
#endif
