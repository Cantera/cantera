/**
 * @file Tortuosity.h
 *  Class to compute the increase in diffusive path length in porous media
 *  assuming the Bruggeman exponent relation
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#ifndef CT_TORTUOSITY_H
#define CT_TORTUOSITY_H

namespace Cantera
{

//! Specific Class to handle tortuosity corrections for diffusive transport
//! in porous media using the Bruggeman exponent
/*!
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 *
 * @deprecated To be removed after Cantera 2.4
 *
 * Class to compute the increase in diffusive path length associated with
 * tortuous path diffusion through, for example, porous media. This base class
 * implementation relates tortuosity to volume fraction through a power-law
 * relationship that goes back to Bruggeman.  The exponent is referred to as the
 * Bruggeman exponent.
 *
 * Note that the total diffusional flux is generally written as
 *
 * \f[
 *   \frac{ \phi C_T D_i \nabla X_i }{ \tau^2 }
 * \f]
 *
 * where \f$ \phi \f$ is the volume fraction of the transported phase,
 * \f$ \tau \f$ is referred to as the tortuosity.  (Other variables are
 * \f$ C_T \f$, the total concentration, \f$ D_i \f$, the diffusion coefficient,
 * and \f$ X_i \f$, the mole fraction with Fickian transport assumed.)
 */
class Tortuosity
{
public:
    //! Default constructor uses Bruggeman exponent of 1.5
    Tortuosity(double setPower = 1.5) : expBrug_(setPower) {
    }

    //! The tortuosity factor models the effective increase in the
    //! diffusive transport length.
    /**
     * This method returns \f$ 1/\tau^2 \f$ in the description of the
     * flux \f$ \phi C_T D_i \nabla X_i / \tau^2 \f$.
     */
    virtual double tortuosityFactor(double porosity) {
        return pow(porosity, expBrug_ - 1.0);
    }

    //! The McMillan number is the ratio of the flux-like
    //! variable to the value it would have without porous flow.
    /**
     * The McMillan number combines the effect of tortuosity and volume fraction
     * of the transported phase.  The net flux observed is then the product of
     * the McMillan number and the non-porous transport rate.  For a
     * conductivity in a non-porous media, \f$ \kappa_0 \f$, the conductivity in
     * the porous media would be \f$ \kappa = (\rm McMillan) \kappa_0 \f$.
     */
    virtual double McMillan(double porosity) {
        return pow(porosity, expBrug_);
    }

protected:
    //! Bruggeman exponent: power to which the tortuosity depends on the volume
    //! fraction
    double expBrug_;
};


/**
 * This class implements transport coefficient corrections appropriate for
 * porous media where percolation theory applies.
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 */
class TortuosityPercolation : public Tortuosity
{
public:
    //! Default constructor uses Bruggeman exponent of 1.5
    TortuosityPercolation(double percolationThreshold = 0.4, double conductivityExponent = 2.0) : percolationThreshold_(percolationThreshold), conductivityExponent_(conductivityExponent) {
    }

    double tortuosityFactor(double porosity) {
        return McMillan(porosity) / porosity;
    }

    double McMillan(double porosity) {
        return pow((porosity - percolationThreshold_)
                   / (1.0 - percolationThreshold_),
                   conductivityExponent_);
    }

protected:
    //! Critical volume fraction / site density for percolation
    double percolationThreshold_;
    //! Conductivity exponent
    /**
     * The McMillan number (ratio of effective conductivity
     * to non-porous conductivity) is
     * \f[
     *     \kappa/\kappa_0 = ( \phi - \phi_c )^\mu
     * \f]
     * where \f$ \mu \f$ is the conductivity exponent (typical values range from
     * 1.6 to 2.0) and \f$ \phi_c \f$ is the percolation threshold.
     */
    double conductivityExponent_;
};


/**
 * This class implements transport coefficient corrections appropriate for
 * porous media with a dispersed phase. This model goes back to Maxwell.  The
 * formula for the conductivity is expressed in terms of the volume fraction of
 * the continuous phase, \f$ \phi \f$, and the relative conductivities of the
 * dispersed and continuous phases, \f$ r = \kappa_d / \kappa_0 \f$.  For dilute
 * particle suspensions the effective conductivity is
 * \f[
 *    \kappa / \kappa_0 = 1 + 3 ( 1 - \phi ) ( r - 1 ) / ( r + 2 ) + O(\phi^2)
 * \f]
 *
 * @attention This class currently does not have any test cases or examples. Its
 *     implementation may be incomplete, and future changes to Cantera may
 *     unexpectedly cause this class to stop working. If you use this class,
 *     please consider contributing examples or test cases. In the absence of
 *     new tests or examples, this class may be deprecated and removed in a
 *     future version of Cantera. See
 *     https://github.com/Cantera/cantera/issues/267 for additional information.
 */
class TortuosityMaxwell : public Tortuosity
{
public:
    //! Default constructor uses Bruggeman exponent of 1.5
    TortuosityMaxwell(double relativeConductivites = 0.0) : relativeConductivites_(relativeConductivites) {
    }

    double tortuosityFactor(double porosity) {
        return McMillan(porosity) / porosity;
    }

    double McMillan(double porosity) {
        return 1 + 3 * (1.0 - porosity) * (relativeConductivites_ - 1.0) / (relativeConductivites_ + 2);
    }

protected:
    //! Relative conductivities of the dispersed and continuous phases,
    //! `relativeConductivites_` \f$ = \kappa_d / \kappa_0 \f$.
    double relativeConductivites_;
};

}
#endif
