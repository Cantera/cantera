/**
 *  @file EdgePhase.h Declarations for the EdgePhase ThermoPhase object, which
 *       models the interface between two surfaces (see @ref thermoprops and
 *       @link Cantera::EdgePhase EdgePhase@endlink).
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
 * EdgePhase object. It should be noted that dimensional quantities with
 * dimensions including a length have their dimensions reduced by one.
 *
 * @ingroup thermoprops
 */
class EdgePhase : public SurfPhase
{
public:
    //! Construct and initialize an EdgePhase directly from an input file
    /*!
     * @param infile name of the input file. If blank, an empty phase will be created.
     * @param id     name of the phase id in the file.
     *               If this is blank, the first phase in the file is used.
     */
    explicit EdgePhase(const string& infile="", const string& id="");

    string type() const override {
        return "edge";
    }
};
}

#endif
