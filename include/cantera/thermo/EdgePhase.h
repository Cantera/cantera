/**
 *  @file EdgePhase.h Declarations for the EdgePhase ThermoPhase object, which
 *       models the interface between two surfaces (see \ref thermoprops and
 *       \link Cantera::EdgePhase EdgePhase\endlink).
 */

//  Copyright 2002 California Institute of Technology

#ifndef CT_EDGEPHASE_H
#define CT_EDGEPHASE_H

#include "SurfPhase.h"

namespace Cantera
{

//! A thermodynamic phase representing a one dimensional edge between two
//! surfaces.
/*!
 * This thermodynamic function is largely a wrapper around the SurfPhase
 * thermodynamic object.
 *
 * All of the equations and formulations carry through from SurfPhase to this
 * EdgePhase object. It should be noted however, that dimensional object with
 * length dimensions, have their dimensions reduced by one.
 *
 * @ingroup thermoprops
 */
class EdgePhase : public SurfPhase
{
public:
    //! Constructor
    /*!
     * @param n0  Surface site density (kmol m-1).
     */
    EdgePhase(doublereal n0=1.0);

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
    EdgePhase(const EdgePhase& right);

    //! Assignment Operator
    /*!
     * @param right Object to be copied
     */
    EdgePhase& operator=(const EdgePhase& right);

    //! Duplicator from a ThermoPhase object
    ThermoPhase* duplMyselfAsThermoPhase() const;

    //! returns the equation of state type
    virtual int eosType() const {
        return cEdge;
    }

    //! Set the Equation-of-State parameters by reading an XML Node Input
    /*!
     * The Equation-of-State data consists of one item, the site density.
     *
     * @param thermoData   Reference to an XML_Node named thermo
     *                     containing the equation-of-state data. The
     *                     XML_Node is within the phase XML_Node describing
     *                     the EdgePhase object.
     *
     * An example of the contents of the thermoData XML_Node is provided
     * below. The units attribute is used to supply the units of the
     * site density in any convenient form. Internally it is changed
     * into MKS form.
     *
     * @code
     *    <thermo model="Edge">
     *       <site_density units="mol/cm"> 3e-15 </site_density>
     *    </thermo>
     * @endcode
     */
    virtual void setParametersFromXML(const XML_Node& thermoData);
};
}

#endif
