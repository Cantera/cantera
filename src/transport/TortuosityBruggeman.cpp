/**
 * @file TortuosityBase.cpp
 *   Base class to compute the increase in diffusive path length associated with
 *   tortuous path diffusion through, for example, porous media.
 */

/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "TortuosityBruggeman.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

//====================================================================================================================
// Default constructor
TortuosityBruggeman::TortuosityBruggeman(doublereal setPower) :
    TortuosityBase(),
    expBrug_(setPower)
{
}
//====================================================================================================================
// Copy Constructor
/*
 * @param right  Object to be copied
 */
TortuosityBruggeman::TortuosityBruggeman(const TortuosityBruggeman& right) :
    TortuosityBase(),
    expBrug_(right.expBrug_)
{
    *this = right;
}
//====================================================================================================================
// Assignment operator
/*
 * @param right Object to be copied
 */
TortuosityBruggeman&   TortuosityBruggeman::operator=(const TortuosityBruggeman& right)
{
    if (&right == this) {
        return *this;
    }
    TortuosityBase::operator=(right);

    expBrug_                = right.expBrug_;

    return *this;
}
//====================================================================================================================
// Duplication operator
/*
 *  @return  Returns a pointer to a duplicate of the current object given a
 *           base class pointer
 */
TortuosityBase* TortuosityBruggeman::duplMyselfAsTortuosityBase() const
{
    return new TortuosityBruggeman(*this);
}
//====================================================================================================================
// The tortuosity factor models the effective increase in the diffusive transport length.
/*
 * This method returns \f$ 1/\tau^2 \f$ in the description of the  flux
 *
 *    \f$  C_T D_i \nabla X_i / \tau^2 \f$.
 *
 *
 */
doublereal TortuosityBruggeman::tortuosityFactor(doublereal porosity)
{
    return pow(porosity, expBrug_ - 1.0);
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
doublereal TortuosityBruggeman::McMillanFactor(doublereal porosity)
{
    return pow(porosity, expBrug_);
}
//====================================================================================================================
}
