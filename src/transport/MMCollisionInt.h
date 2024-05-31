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
/*!
 * This class provides functions that interpolate the tabulated collision
 * integrals in Monchick and Mason @cite monchick1961.
 *
 * The collision integrals computed by Monchick and Mason use the Stockmayer
 * potential, which models a polar molecule as a spherical potential with a
 * point dipole at the center).
 *
 * @f[
 *    \phi(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 + \delta \left(\frac{\sigma}{r}\right)^3 \right]
 * @f]
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
    void fit(int degree, double deltastar, double* astar, double* bstar, double* cstar);
    void fit_omega22(int degree, double deltastar, double* om22);
    double omega11(double ts, double deltastar) {
        return omega22(ts, deltastar)/astar(ts, deltastar);
    }

private:
    double fitDelta(int table, int ntstar, int degree, double* c);
    double quadInterp(double x0, double* x, double* y);

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
