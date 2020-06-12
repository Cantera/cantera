/**
 *  @file TransportParams.h
 *  Class that holds the data that is read in from the input file, and which is used for
 *  processing of the transport object
 *  (see \ref tranprops and \link Cantera::TransportParams TransportParams \endlink).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_TRANSPORTPARAMS_H
#define CT_TRANSPORTPARAMS_H

#include "cantera/numerics/DenseMatrix.h"
#include "TransportBase.h"

namespace Cantera
{

//! Base structure to hold transport model parameters.
/*!
 * This structure is used by TransportFactory.
 *
 * @deprecated Unused. To be removed after Cantera 2.5.
 */
class TransportParams
{
public:
    TransportParams();
    virtual ~TransportParams() {}

    //! Local storage of the number of species
    size_t nsp_;

    //!  Pointer to the ThermoPhase object: shallow pointer
    thermo_t* thermo;

    //! Local storage of the molecular weights of the species
    /*!
     *  Length is nsp_ and units are kg kmol-1.
     */
    vector_fp mw;

    //! A basis for the average velocity can be specified.
    /*!
     *  Valid bases include "mole", "mass", and "species" names.
     */
    VelocityBasis velocityBasis_;

    //! Maximum temperatures for parameter fits
    doublereal tmax;

    //! Minimum temperatures for parameter fits
    doublereal tmin;

    //!  Mode parameter
    int mode_;

    //! Log level
    int log_level;
};

} // End of namespace Cantera

#endif //CT_TRANSPORTPARAMS_H
