/**
 * @file TortuosityBase.h
 * Virtual base class to compute the increase in diffusive path length associated with
 * tortuous path diffusion through, for example, porous media.
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef CT_TORTUOSITYPERCOLATION_H
#define CT_TORTUOSITYPERCOLATION_H

#include "TortuosityBase.h"

namespace Cantera
{

//!  This class implements transport coefficient corrections
//! appropriate for porous media where percolation theory applies.
class TortuosityPercolation : public TortuosityBase
{

public:
    //! Default constructor uses Percolation exponent of 1.5
    /*!
     *  @param setPower       Exponent in the Percolation factor. The default is 1.5
     */
    TortuosityPercolation(double percolationThreshold = 0.4, double conductivityExponent = 2.0);

    //! Copy Constructor
    /*!
     * @param right  Object to be copied
     */
    TortuosityPercolation(const TortuosityPercolation& right);

    //! Assignment operator
    /*!
     * @param right Object to be copied
     */
    TortuosityPercolation& operator=(const TortuosityPercolation& right);

    //! Duplication operator
    /*!
     *  @return  Returns a pointer to a duplicate of the current object given a
     *           base class pointer
     */
    virtual TortuosityBase* duplMyselfAsTortuosityBase() const;

    //! The tortuosity factor models the effective increase in the
    //! diffusive transport length.
    /*!
     * This method returns \f$ 1/\tau^2 \f$ in the description of the  flux
     *
     *    \f$  C_T D_i \nabla X_i / \tau^2 \f$.
     *
     *
     */
    virtual doublereal tortuosityFactor(doublereal porosity);

    //! The McMillan number is the ratio of the flux-like
    //! variable to the value it would have without porous flow.
    /*!
     * The McMillan number combines the effect of tortuosity
     * and volume fraction of the transported phase.  The net flux
     * observed is then the product of the McMillan number and the
     * non-porous transport rate.  For a conductivity in a non-porous
     * media, \f$ \kappa_0 \f$, the conductivity in the porous media
     * would be \f$ \kappa = (\rm McMillan) \kappa_0 \f$.
     */
    virtual doublereal McMillanFactor(doublereal porosity);


protected:

    //! Critical volume fraction / site density for percolation
    double percolationThreshold_;

    //! Conductivity exponent
    /*!
     *   The McMillan number (ratio of effective conductivity  to non-porous conductivity) is
     *           \f[ \kappa/\kappa_0 = ( \phi - \phi_c )^\mu \f]
     *   where \f$ \mu \f$ is the conductivity exponent (typical values range from 1.6 to 2.0) and \f$ \phi_c \f$
     *   is the percolation threshold.
     */
    double conductivityExponent_;


};



}

#endif

