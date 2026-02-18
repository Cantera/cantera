/**
 * @file MMCollisionInt.h
 *  Monchick and Mason collision integrals
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_MMCOLLISIONINT_H
#define CT_MMCOLLISIONINT_H

#include "cantera/base/ct_defs.h"

namespace Cantera
{

//! Calculation of Collision integrals
/**
 * This class provides functions that interpolate the tabulated collision integrals in
 * Monchick and Mason @cite monchick1961.
 *
 * The collision integrals computed by Monchick and Mason use the Stockmayer potential,
 * which models a polar molecule as a spherical potential with a point dipole at the
 * center). Equation 16 of Monchick and Mason @cite monchick1961 gives the potential
 * as:
 *
 * @f[
 *  \phi(r) = 4 \epsilon_0 \left[ \left(\frac{\sigma_0}{r}\right)^{12} - \left(\frac{\sigma_0}{r}\right)^6 + \delta \left(\frac{\sigma_0}{r}\right)^3 \right]
 * @f]
 *
 * Where @f$ \epsilon_0 @f$ is the depth of the potential well, @f$ \sigma_0 @f$ is the
 * distance at which the potential between the two molecules is zero, @f$ \delta @f$ is
 * defined as:
 *
 * @f[
 *  \delta = \frac{1}{4} (\mu^*)^2 \zeta \left( \theta_1, \theta_2, \phi \right)
 * @f]
 *
 * @f$ \mu^* @f$ is the reduced dipole moment. @f$ \theta_1 @f$ , @f$ \theta_2 @f$ ,
 * and @f$ \phi @f$ are angles related to trajectories of colliding molecules. In the
 * work of Monchick and Mason, these details are not what is presented in the tables.
 * Instead, the tables are presented as being functions of the reduced temperature,
 * @f$ T^* @f$, and the @f$ \delta @f$ parameter. The reduced dipole moment,
 * @f$ \mu^* @f$ is defined as:
 *
 * @f[
 *  \mu^* = \frac{\mu}{\sqrt{\epsilon_0 \sigma_0^3}}
 * @f]
 *
 * Where @f$ \mu @f$ is the dipole moment of the molecule and the other parameters
 * have been defined earlier. This work considers only the collisions of like
 * molecules, so only a single value is needed.
 *
 * The tabulated data comes from the averaged collision integrals in tables
 * V through VIII of Monchick and Mason @cite monchick1961.
 *
 * @ingroup tranprops
 */
class MMCollisionInt
{
public:
    MMCollisionInt() {}
    virtual ~MMCollisionInt() {}

    //! Initialize the object for calculation
    /*!
     *  @param tsmin       Minimum value of Tstar to carry out the fitting
     *  @param tsmax       Maximum value of Tstar to carry out the fitting
     */
    void init(double tsmin, double tsmax);

    double omega22(double ts, double deltastar);
    double astar(double ts, double deltastar);
    double bstar(double ts, double deltastar);
    double cstar(double ts, double deltastar);
    void fit(int degree, double deltastar, span<double> astar, span<double> bstar,
             span<double> cstar);
    void fit_omega22(int degree, double deltastar, span<double> om22);
    double omega11(double ts, double deltastar) {
        return omega22(ts, deltastar)/astar(ts, deltastar);
    }

private:
    double fitDelta(int table, int ntstar, int degree, span<double> c);
    double quadInterp(double x0, span<const double> x, span<const double> y);

    vector<vector<double>> m_o22poly;
    vector<vector<double>> m_apoly;
    vector<vector<double>> m_bpoly;
    vector<vector<double>> m_cpoly;

    static double delta[8];
    static double tstar22[37];

    //! Table of omega22 values
    static double omega22_table[37*8];

    //! T* values (reduced temperature)
    static double tstar[39];

    //! astar table
    static double astar_table[39*8];

    //! bstar table
    static double bstar_table[39*8];

    //! cstar table
    static double cstar_table[39*8];

    //! Log temp
    vector<double> m_logTemp;

    //! Index of the tstar array that encompasses the minimum temperature
    //! fitting range value of tsmin.
    int m_nmin;

    //! Index of the tstar array that encompasses the maximum temperature
    //! fitting range value of tsmax.
    int m_nmax;
};

}
#endif
