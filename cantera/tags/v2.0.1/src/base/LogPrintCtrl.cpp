/**
 * @file LogPrintCtrl.cpp
 *    Declarations for a simple class that augments the logfile printing capabilities
 *   (see \ref Cantera::LogPrintCtrl).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include <cmath>

#include <iostream>
#include <fstream>

#include "LogPrintCtrl.h"
#include "cantera/base/global.h"

using namespace std;

namespace Cantera
{

LogPrintCtrl::LogPrintCtrl(int Ndec) :
    m_ffss(0),
    m_pc(0)
{
    m_ffss = new std::ostream(m_os.rdbuf());
    m_pc = new PrintCtrl(*m_ffss, Ndec);
}

LogPrintCtrl::~LogPrintCtrl()
{
    delete m_pc;
    delete m_ffss;
}

// Print a double using scientific notation
/*
 * Prints a double using scientific notation in a
 * fixed number of spaces
 *
 *
 *  @param d  double to be printed
 *  @param w  Number of spaces to use
 *  @param p  Precision
 *
 *
 */
void LogPrintCtrl::pr_de_c10(const double din, int p, const int wMin,
                             const int wMax)
{
    m_pc->pr_de_c10(din, p, wMin, wMax);
    writelog(m_os.str());
    m_os.str("");
}

// Print a double using scientific notation
/*
 * Prints a double using scientific notation in a
 * fixed number of spaces. Rounding of the last digit is carried out
 * by the standard c++ printing utilities.
 *
 *  @param d  double to be printed
 *  @param w  Number of spaces to use
 *  @param p  Precision
 */
void LogPrintCtrl::pr_de(const double d, int sigDigIn, const int wMinIn,
                         const int wMaxIn)
{
    m_pc->pr_de(d, sigDigIn, wMinIn, wMaxIn);
    writelog(m_os.str());
    m_os.str("");
}

// Croup a double at a certain decade level
/*
 *    This routine will crop a floating point number at a certain
 *  decade lvl. In other words everything below a power of 10^Ndec
 *  will be deleted.
 *  Note, it currently does not do rounding of the last digit.
 *
 *   @param d Double to be cropped
 *   @param nSig Number of significant digits
 *   example:
 *    d = 1.1305E-15;
 *    Ndec = -16;
 *   This routine will return 1.1E-15
 *
 *    d = 8.0E-17
 *    Ndec = -16
 *   This routine will return 0.0
 */
double LogPrintCtrl::cropAbs10(const double d, int Ndec) const
{
    return m_pc->cropAbs10(d, Ndec);
}

// Crop a double at a certain number of significant digits
/*
 *  This routine will crop a floating point number at a certain
 *  number of significant digits. Note, it currently does
 *  rounding up of the last digit.
 *
 *  example:
 *    d = 1.0305E-15;
 *    nsig = 3;
 *   This routine will return 1.03E-15
 */
double LogPrintCtrl::cropSigDigits(const double d, int nSig) const
{
    return m_pc->cropSigDigits(d, nSig);
}

// Set the default value of N decade
/*
 * @param Ndec new value of Ndec
 *
 * @return returns the old value of Ndec
 */
int LogPrintCtrl::setNdec(int Ndec)
{
    return m_pc->setNdec(Ndec);
}

// Set the default significant digits to output
/*
 * @param nSigDigits new value of the sig digits
 *
 * @return returns the old value of Ndec
 */
int LogPrintCtrl::setSigDigits(int nSigDigits)
{
    return m_pc->setSigDigits(nSigDigits);
}

// Set the default minimum width
/*
 * @param wmin Default minimum width
 *
 * @return returns the old default
 */
int LogPrintCtrl::setWmin(int wmin)
{
    return m_pc->setWmin(wmin);
}


// Set the default maximum width
/*
 * @param wmin Default maximum width
 *
 * @return returns the old default
 */
int LogPrintCtrl::setWmax(int wmax)
{
    return m_pc->setWmax(wmax);
}


}
