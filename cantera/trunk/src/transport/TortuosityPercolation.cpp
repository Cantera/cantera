/**
 * @file TortuosityPercolation.cpp
 *   Base class to compute the increase in diffusive path length associated with
 *   tortuous path diffusion through, for example, porous media.
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "TortuosityPercolation.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

//====================================================================================================================
// Default constructor
TortuosityPercolation::TortuosityPercolation(double percolationThreshold, double conductivityExponent) :
    TortuosityBase(),
    percolationThreshold_(percolationThreshold),
    conductivityExponent_(conductivityExponent)
{

}
//====================================================================================================================
// Copy Constructor
/*
 * @param right  Object to be copied
 */
TortuosityPercolation::TortuosityPercolation(const TortuosityPercolation& right) :
    TortuosityBase(),
    percolationThreshold_(right.percolationThreshold_),
    conductivityExponent_(right.conductivityExponent_)
{
    *this = right;
}
//====================================================================================================================
// Assignment operator
/*
 * @param right Object to be copied
 */
TortuosityPercolation& TortuosityPercolation::operator=(const TortuosityPercolation& right)
{
    if (&right == this) {
        return *this;
    }
    TortuosityBase::operator=(right);

    percolationThreshold_   = right.percolationThreshold_;
    conductivityExponent_   = right.conductivityExponent_;

    return *this;
}
//====================================================================================================================
// Duplication operator
/*
 *  @return  Returns a pointer to a duplicate of the current object given a
 *           base class pointer
 */
TortuosityBase* TortuosityPercolation::duplMyselfAsTortuosityBase() const
{
    return new TortuosityPercolation(*this);
}
//====================================================================================================================
// The tortuosity factor models the effective increase in the diffusive transport length.
/*
 * This method returns \f$ 1/\tau^2 \f$ in the description of the  flux
 *
 *    \f$  C_T D_i \nabla X_i / \tau^2 \f$.
 *
 */
doublereal TortuosityPercolation::tortuosityFactor(doublereal porosity)
{
    return McMillanFactor(porosity) / porosity;
}
//====================================================================================================================
// The McMillan number is the ratio of the flux-like variable to the value it would have without porous flow.
/*
 * The McMillan number combines the effect of tortuosity
 * and volume fraction of the transported phase.  The net flux
 * observed is then the product of the McMillan number and the
 * non-porous transport rate.  For a conductivity in a non-porous
 * media, \f$ \kappa_0 \f$, the conductivity in the porous media
 * would be \f$ \kappa = (\rm McMillan) \kappa_0 \f$.
 */
doublereal TortuosityPercolation::McMillanFactor(doublereal porosity)
{
    doublereal tmp = pow(((porosity - percolationThreshold_)
                          / (1.0 - percolationThreshold_)) ,
                         conductivityExponent_);
    return tmp;
}
//====================================================================================================================
}
