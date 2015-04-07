/**
 * @file LogPrintCtrl.h
 *    Declarations for a simple class that augments the logfile printing capabilities
 *   (see \ref Cantera::LogPrintCtrl).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_LOGPRINTCTRL_H
#define CT_LOGPRINTCTRL_H

#include <sstream>

#include "cantera/base/PrintCtrl.h"

namespace Cantera
{

//! This class provides some printing and cropping utilities
//! for writing to the logfile.
/*!
 *  This class writes its output to Cantera's logfile
 *  utility. It's a wrapper around PrintCtrl object.
 *  First, we direct PrintCtrl to write to a string
 *  and then we redirect the string to the logfile utility.
 *  This is a first cut, and it's pretty much a kluge.
 *  The logfile utility, however, demands a string, and this
 *  is what I came up with.
 *
 * @ingroup globalUtilFuncs
 * @deprecated To be removed in Cantera 2.2.
 */
class LogPrintCtrl
{
public:
    //! Constructor
    /*!
     * This also serves to initialize the ticks within the object
     *
     * @param Ndec value of Ndec. Defaults to -1000, i.e., no decade cropping
     */
    LogPrintCtrl(int Ndec = -1000);

    //! Destructor
    ~LogPrintCtrl();

    //! Print a double using scientific notation
    /*!
     * Prints a double using scientific notation in a
     * fixed number of spaces.
     *
     * The precision of the number will be adjusted to
     * fit into the maximum space.
     *
     *  @param d  double to be printed
     *  @param sigDigits Number of significant digits (-1 = default, means to
     *          use the default number for the object, which is initially set
     *          to 13.
     *  @param wMin Minimum number of spaces to print out
     *  @param wMax Maximum number of spaces to print out
     */
    void pr_de(const double d, int sigDigits = -1,
               const int wMin = -1, const int wMax = -1);

    //! Print a double using scientific notation cropping
    //! decade values
    /*!
     * Prints a double using scientific notation in a
     * fixed number of spaces. This routine also crops
     * number below the default decade level.
     *
     * The precision of the number will be adjusted to
     * fit into the maximum space.
     *
     *  @param d  double to be printed
     *  @param sigDigits Number of significant digits (-1 = default, means to
     *          use the default number for the object, which is initially set
     *          to 13.
     *  @param wMin Minimum number of spaces to print out
     *  @param wMax Maximum number of spaces to print out
     */
    void pr_de_c10(const double d, int sigDigits = -1,
                   const int wMin = -1, const int wMax = -1);

    //! Crop a double at a certain number of significant digits
    /*!
     *  This routine will crop a floating point number at a certain
     *  number of significant digits. Note, it does
     *  rounding up of the last digit.
     *
     * @param d         Double to be cropped
     * @param sigDigits Number of significant digits
     *  example:
     *    d = 1.0305E-15;
     *    nsig = 3;
     *   This routine will return 1.03E-15
     */
    double cropSigDigits(const double d, int sigDigits) const;

    //! Crop a double at a certain decade level
    /*!
     *  This routine will crop a floating point number at a certain decade
     *  lvl. In other words everything below a power of 10^Ndec will be
     *  deleted. Note, it does rounding up of the last digit.
     *
     *   @param d Double to be cropped
     *   @param nDecades Number of significant digits
     *   example:
     *    d = 1.1305E-15;
     *    nDecades = -16;
     *   This routine will return 1.1E-15
     *
     *    d = 8.0E-17
     *    nDecades = -16
     *   This routine will return 0.0
     */
    double cropAbs10(const double d, const int nDecades) const;

    //! Set the default value of N decade
    /*!
     * @param nDecades new value of Ndec
     * @return returns the old value of Ndec
     */
    int setNdec(int nDecades);


    //! Set the default significant digits to output
    /*!
     * @param sigDigits new value of the sig digits
     * @return returns the old value of Ndec
     */
    int setSigDigits(int sigDigits);

    //! Set the default minimum width
    /*!
     * @param wMin Default minimum width
     * @return returns the old default
     */
    int setWmin(int wMin);

    //! Set the default maximum width
    /*!
     * @param wMax Default maximum width
     * @return returns the old default
     */
    int setWmax(int wMax);

private:
    //! local stringstream class for temp output
    std::ostringstream m_os;

    //! Pointer to the ostream where this class actually
    //! prints its information
    std::ostream* m_ffss;

    //! Pointer to the PrintCtrl class
    PrintCtrl* m_pc;
};
}

#endif
