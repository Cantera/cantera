/**
 *  @file UnityLewisTransport.h
 *    Headers for the UnityLewisTransport object, which extends the MixTransport
 *    object, which models transport properties in ideal gas solutions using a 
 *    mixture averaged approximation. This object assumes Lewis numbers to be 
 *    unity for all species.
 *    (see \ref tranprops and \link Cantera::MixTransport MixTransport \endlink) .
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_UNITYLEWISTRAN_H
#define CT_UNITYLEWISTRAN_H

#include "GasTransport.h"
#include "MixTransport.h"
#include "cantera/numerics/DenseMatrix.h"

namespace Cantera
{
//! Class UnityLewisTransport implements mixture-averaged transport properties for
//! ideal gas mixtures.
class UnityLewisTransport : public MixTransport
{
public:
    //! Default constructor.
    UnityLewisTransport();

    virtual void init(thermo_t* thermo, int mode=0, int log_level=0);

    virtual std::string transportType() const {
        return "UnityLewisNumber";
    }

    virtual void getMixDiffCoeffs(doublereal* const d);

    virtual void getMixDiffCoeffsMole(doublereal* const d);

    virtual void getMixDiffCoeffsMass(doublereal* const d);


private:

};

}
#endif
