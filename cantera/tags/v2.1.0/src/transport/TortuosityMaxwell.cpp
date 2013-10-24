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

#include "TortuosityMaxwell.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

//====================================================================================================================
// Default constructor
TortuosityMaxwell::TortuosityMaxwell(doublereal relativeConductivities) :
    TortuosityBase(),
    relativeConductivities_(relativeConductivities)
{
}
//====================================================================================================================
// Copy Constructor
/*
 * @param right  Object to be copied
 */
TortuosityMaxwell::TortuosityMaxwell(const TortuosityMaxwell& right) :
    TortuosityBase(),
    relativeConductivities_(right.relativeConductivities_)
{
    *this = right;
}
//====================================================================================================================
// Assignment operator
/*
 * @param right Object to be copied
 */
TortuosityMaxwell&   TortuosityMaxwell::operator=(const TortuosityMaxwell& right)
{
    if (&right == this) {
        return *this;
    }
    TortuosityBase::operator=(right);

    relativeConductivities_ = right.relativeConductivities_;

    return *this;
}
//====================================================================================================================
// Duplication operator
/*
 *  @return  Returns a pointer to a duplicate of the current object given a
 *           base class pointer
 */
TortuosityBase* TortuosityMaxwell::duplMyselfAsTortuosityBase() const
{
    return new TortuosityMaxwell(*this);
}
//====================================================================================================================
// The tortuosity factor models the effective increase in the diffusive transport length.
/*
 * This method returns \f$ 1/\tau^2 \f$ in the description of the  flux
 *
 *    \f$  C_T D_i \nabla X_i / \tau^2 \f$.
 *
 */
doublereal TortuosityMaxwell::tortuosityFactor(doublereal porosity)
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
doublereal TortuosityMaxwell::McMillanFactor(doublereal porosity)
{
    return 1 + 3 * (1.0 - porosity) * (relativeConductivities_ - 1.0) / (relativeConductivities_ + 2);
}
//====================================================================================================================
}
