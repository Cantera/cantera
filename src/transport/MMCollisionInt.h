/**
 * @file MMCollisionInt.h
 *  Monk and Monchick collision integrals
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_MMCOLLISIONINT_H
#define CT_MMCOLLISIONINT_H

#include "cantera/base/ct_defs.h"

#include <vector>
#include <iostream>

namespace Cantera
{

class XML_Writer;

//! Calculation of Collision integrals
/*!
 * This class provides functions that
 * interpolate the tabulated collision integrals in Monchick and
 * Mason, "Transport Properties of Polar Gases," J. Chem. Phys. (1961)
 *
 * @ingroup transportgroup
 */
class MMCollisionInt
{

public:

    //! Default Constructor
    MMCollisionInt();

    //! Destructor
    virtual ~MMCollisionInt();


    //! Initialize the object for calculation
    /*!
     *
     *  @param xml         Pointer to the log file that will receive the debug output
     *                     messages
     *  @param tsmin       Minimum value of Tstar to carry out the fitting
     *  @param tsmax       Maximum value of Tstar to carry out the fitting
     *  @param loglevel    Set the loglevel for the object. The default
     *                     loglevel is zero, indicating no output.
     */
    void init(XML_Writer* xml, doublereal tsmin,  doublereal tsmax, int loglevel = 0);

    //! omega22
    /*!
     * @param ts
     *  @param deltastar
     */
    doublereal omega22(double ts, double deltastar);

    //! astar
    /*!
     * @param ts
     *  @param deltastar
     */
    doublereal astar(double ts, double deltastar);

    //! bstar
    /*!
     * @param ts
     *  @param deltastar
     */
    doublereal bstar(double ts, double deltastar);

    //! cstar
    /*!
     * @param ts
     *  @param deltastar
     */
    doublereal cstar(double ts, double deltastar);

    //! fit
    /*!
     *  @param logfile
     *  @param degree
     *  @param deltastar
     *  @param astar
     *  @param bstar
     *  @param cstar
     */
    void fit(std::ostream& logfile, int degree, doublereal deltastar,
             doublereal* astar, doublereal* bstar, doublereal* cstar);

    //! fit_omega22
    /*!
     *    @param logfile
     *    @param degree
     *    @param deltastar
     *    @param om22
     */
    void fit_omega22(std::ostream& logfile, int degree, doublereal deltastar, doublereal* om22);

    //! omega11
    /*!
     *  @param ts
     *  @param deltastar
     */
    doublereal omega11(double ts, double deltastar) {
        return omega22(ts, deltastar)/astar(ts, deltastar);
    }

private:

    //! Fit delta
    /*!
     *  @param table
     *  @param ntstar
     *  @param degree
     *  @param c         C is probable the output vector
     *
     *  @return
     */
    doublereal fitDelta(int table, int ntstar, int degree, doublereal* c);

    //! m_o22poly
    std::vector<vector_fp>  m_o22poly;

    //! m_apoly
    std::vector<vector_fp>  m_apoly;
    //! m_bpoly
    std::vector<vector_fp>  m_bpoly;

    //! m_cpoly
    std::vector<vector_fp>  m_cpoly;

    //! delta
    static doublereal delta[8];

    //! tstar22
    static doublereal tstar22[37];

    //! Table of omega22 values from MM
    static doublereal omega22_table[37*8];

    //! tstar
    /*!
     *   table of tstar values
     */
    static doublereal tstar[39];

    //! astar table from MM
    static doublereal astar_table[39*8];

    //! bstar table from MM
    static doublereal bstar_table[39*8];

    //! cstar table from MM
    static doublereal cstar_table[39*8];

    //! Log temp
    vector_fp  m_logTemp;

    //! nmin
    int m_nmin;

    //! nmax
    int m_nmax;

    //! XML_Writer pointer
    XML_Writer* m_xml;

    //! loglevel
    int m_loglevel;
};
}
#endif


