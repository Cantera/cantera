/**
 *  @file LiquidTransportData.h
 *  Header file defining class LiquidTransportData
 */
#ifndef CT_LIQUIDTRANSPORTDATA_H
#define CT_LIQUIDTRANSPORTDATA_H

#include "TransportBase.h"
#include "cantera/base/FactoryBase.h"
#include "LTPspecies.h"

namespace Cantera
{

//! Class LiquidTransportData holds transport parameters for a
//! specific liquid-phase species.
/*!
 * A LiquidTransportData object is created for each species.
 *
 * This class is mainly used to collect transport properties from the parse
 * phase in the TransportFactory and transfer them to the Transport class.
 * Transport properties are expressed by subclasses of LTPspecies. One may
 * need to be careful about deleting pointers to LTPspecies objects created in
 * the TransportFactory.
 *
 *  All of the pointers in this class are shallow pointers. Therefore, this
 *  is a passthrough class, which keeps track of pointer ownership by zeroing
 *  pointers as we go. Yes, Yes, yes, this is not good.
 */
class LiquidTransportData
{
public:
    LiquidTransportData();
    LiquidTransportData(const LiquidTransportData& right);
    LiquidTransportData& operator=(const LiquidTransportData& right);
    ~LiquidTransportData();

    //! A LiquidTransportData object is instantiated for each species.
    //! This is the species name for which this object is instantiated.
    std::string speciesName;

    //! Model type for the hydroradius
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* hydroRadius;

    //! Model type for the viscosity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* viscosity;

    //! Model type for the ionic conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* ionConductivity;

    //! Model type for the mobility ratio
    /*!
     *  shallow pointers that should be zero during destructor
     */
    std::vector<LTPspecies*> mobilityRatio;

    //! Model type for the self diffusion coefficients
    /*!
     *  shallow pointers that should be zero during destructor
     */
    std::vector<LTPspecies*> selfDiffusion;

    //! Model type for the thermal conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* thermalCond;

    //! Model type for the electrical conductivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* electCond;

    //! Model type for the speciesDiffusivity
    /*!
     *  shallow pointer that should be zero during destructor
     */
    LTPspecies* speciesDiffusivity;
};

}
#endif
