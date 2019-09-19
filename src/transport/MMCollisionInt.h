/**
 * @file MMCollisionInt.h
 *  Monk and Monchick collision integrals
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
 * integrals in Monchick and Mason, "Transport Properties of Polar Gases," J.
 * Chem. Phys. (1961)
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
     *  @param loglevel    Set the loglevel for the object. The default
     *                     loglevel is zero, indicating no output.
     */
    void init(doublereal tsmin, doublereal tsmax, int loglevel = 0);

    doublereal omega22(double ts, double deltastar);
    doublereal astar(double ts, double deltastar);
    doublereal bstar(double ts, double deltastar);
    doublereal cstar(double ts, double deltastar);
    void fit(int degree, doublereal deltastar,
             doublereal* astar, doublereal* bstar, doublereal* cstar);
    void fit_omega22(int degree, doublereal deltastar, doublereal* om22);
    doublereal omega11(double ts, double deltastar) {
        return omega22(ts, deltastar)/astar(ts, deltastar);
    }

private:
    doublereal fitDelta(int table, int ntstar, int degree, doublereal* c);

    std::vector<vector_fp> m_o22poly;
    std::vector<vector_fp> m_apoly;
    std::vector<vector_fp> m_bpoly;
    std::vector<vector_fp> m_cpoly;

    static doublereal delta[8];
    static doublereal tstar22[37];

    //! Table of omega22 values from MM
    static doublereal omega22_table[37*8];

    //! table of tstar values
    static doublereal tstar[39];

    //! astar table from MM
    static doublereal astar_table[39*8];

    //! bstar table from MM
    static doublereal bstar_table[39*8];

    //! cstar table from MM
    static doublereal cstar_table[39*8];

    //! Log temp
    vector_fp m_logTemp;

    int m_nmin;
    int m_nmax;

    //! loglevel
    int m_loglevel;
};
}
#endif
