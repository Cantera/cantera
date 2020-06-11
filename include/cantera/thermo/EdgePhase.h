/**
 *  @file EdgePhase.h Declarations for the EdgePhase ThermoPhase object, which
 *       models the interface between two surfaces (see \ref thermoprops and
 *       \link Cantera::EdgePhase EdgePhase\endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

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

    virtual std::string type() const {
        return "Edge";
    }

    //! Set the Equation-of-State parameters by reading an XML Node Input
    /*!
     * The Equation-of-State data consists of one item, the site density.
     *
     * @param thermoData   Reference to an XML_Node named thermo containing the
     *                     equation-of-state data. The XML_Node is within the
     *                     phase XML_Node describing the EdgePhase object.
     *
     * An example of the contents of the thermoData XML_Node is provided below.
     * The units attribute is used to supply the units of the site density in
     * any convenient form. Internally it is changed into MKS form.
     *
     * @code
     *    <thermo model="Edge">
     *       <site_density units="mol/cm"> 3e-15 </site_density>
     *    </thermo>
     * @endcode
     *
     * @deprecated The XML input format is deprecated and will be removed in
     *     Cantera 3.0.
     */
    virtual void setParametersFromXML(const XML_Node& thermoData);
};
}

#endif
