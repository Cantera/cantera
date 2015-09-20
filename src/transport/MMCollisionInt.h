/**
 * @file MMCollisionInt.h
 *  Monk and Monchick collision integrals
 */
// Copyright 2001  California Institute of Technology

#ifndef CT_MMCOLLISIONINT_H
#define CT_MMCOLLISIONINT_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <map>
#include <iomanip>

#include "cantera/base/ct_defs.h"
#include "StringFunct.h"
#include "TransportReadExt.h"
#include "TransportCharged.h"


using namespace std;

namespace Cantera
{

//! Calculation of Collision integrals
/*!
 * This class provides functions that
 * interpolate the tabulated collision integrals in Monchick and
 * Mason, "Transport Properties of Polar Gases," J. Chem. Phys. (1961)
 *
 * @ingroup tranprops
 *
 * NOTE: new functions have been added for a more accurate description at high temperature
 * and for ionized mixtures
 */
class MMCollisionInt
{
public:
    MMCollisionInt();
    virtual ~MMCollisionInt();

    //! Initialize the object for calculation
    /*!
     *  @param tsmin       Minimum value of Tstar to carry out the fitting
     *  @param tsmax       Maximum value of Tstar to carry out the fitting
     *  @param loglevel    Set the loglevel for the object. The default
     *                     loglevel is zero, indicating no output.
     */
    void init(doublereal tsmin,  doublereal tsmax, int loglevel = 0);

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

    // Collision integrals with (m,6)-potential and Born-Meyer involving neutral species --> _hT
    // Collision integrals involving charged species --> _charged
    doublereal omega11_hT(const string* species, double T, int nS);
    doublereal omega22_hT(const string* species, double T, int nS);
    doublereal omega11_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega22_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega12_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega13_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega14_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega15_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega23_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal omega24_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);
    doublereal bstar_hT(const string* species, double T, int nS);
    doublereal bstar_charged(const string species1, const string species2, int charge1, int charge2, double T, double Te, double Xe, double P);

    void getCoefficients_omega11(string species1, string species2, vector<double> &coefficients, double T);
    void getCoefficients_omega22(string species1, string species2, vector<double> &coefficients, double T);

    // in the external file the minus sign (-) is used to separate collisional pairs
    // this function locally replace the minus sign with m
    void replaceMinus(string &species1, string &species2);
    void getCoefficients_bstar(string species1, string species2, vector<double> &coefficients, double T);

    vector<double> C_MMomega11;
    vector<double> C_MMomega22;
    vector<string> PairsMM;
    vector<double> C_MMbstar;


private:
    doublereal fitDelta(int table, int ntstar, int degree, doublereal* c);

    std::vector<vector_fp>  m_o22poly;

    std::vector<vector_fp>  m_apoly;
    std::vector<vector_fp>  m_bpoly;

    std::vector<vector_fp>  m_cpoly;

    static doublereal delta[8];

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

    int m_nmin;

    int m_nmax;

    //! loglevel
    int m_loglevel;
};
}
#endif
