/**
 * @file PrintCtrl.h
 *    Declarations for a simple class that augments the streams printing capabilities
 *   (see \ref Cantera::PrintCtrl).
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef CT_PRINTCTRL_H
#define CT_PRINTCTRL_H

#include <iostream>

namespace Cantera
{

//! This class provides some printing and cropping utilities
/*!
 *  The class is used to provide some formatting options for
 *  printing out real numbers to files and to standard output.
 *  Specifically, it can make sure that a max and min field
 *  width is honored when conducting IO of numbers and strings.
 *  Basically, it's the spot to house all wrappers around
 *  commonly used printing facilities.
 *
 *  It can also handle cropping of numbers below a certain
 *  decade level. This is useful for IO for testing purposes.
 *  For example, if you don't care about anything below
 *  1.0E-20, you can set up the IO so that it won't print out
 *  any digits below 1.0E-20, even digits that are in numbers
 *  greater than 1.0E-20. In other words the number
 *
 *  1.12345E-19
 *
 *  would be cropped to the value
 *
 *  1.1000E-19
 *
 *  The class wraps around a single std::ostream class. Its
 *  cropping functions are also available as a "double"
 *  conversion utility.
 *
 * @ingroup globalUtilFuncs
 * @deprecated To be removed in Cantera 2.2.
 */
class PrintCtrl
{
public:
    //! enum for cropping control
    enum CROP_TYPE {
        //! Turn off cropping always
        CT_OFF=0,
        //! Turn off cropping, unless the global toggle is turned on
        CT_OFF_GLOBALOBEY,
        //! Turn on cropping unless the global toggle is turned off
        CT_ON_GLOBALOBEY,
        //! Turn on cropping always
        CT_ON
    };

    //! enum for global cropping control
    enum CROP_TYPE_GLOBAL {
        //! no preference for global cropping
        GCT_NOPREF = 0,
        //! global toggle for turning on cropping
        GCT_CROP,
        //! global toggle for turning off cropping
        GCT_NOCROP
    };

    //! static enum for turning on and off cropping
    /*!
     *  The default is to not have a preference for cropping
     */
    static CROP_TYPE_GLOBAL GlobalCrop;

    //! Constructor
    /*!
     * This also serves to initialize the ticks within the object
     *
     * @param coutProxy  This is a reference to the ostream
     *                   to use for all IO from ths object.
     * @param Ndec value of Ndec. Defaults to -1000, i.e., no decade cropping
     * @param ctlocal    The default is to turn on cropping all the time.
     */
    PrintCtrl(std::ostream& coutProxy = std::cout, int Ndec = -1000,
              CROP_TYPE ctlocal = CT_ON);

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

    //! Set the cropping control flag
    /*!
     * @param ctlocal Local enum value for the cropping type
     */
    void setCropCntrl(CROP_TYPE  ctlocal);

private:
    //! private function to figure out cropping logic
    /*!
     *  @return Returns the decision as to whether to crop or not
     */
    bool doCrop() const;

    //! This is the ostream to send all output from the object
    /*!
     * It defaults to cout
     */
    std::ostream& m_cout;

    //! Default decade level to use for decade cropping
    /*!
     * This is initially set to -1000, which means that
     * no cropping will be carried out
     */
    int m_Ndec;

    //! default precision level to use in printing
    /*!
     * This actually is one less than the number of significant digits.
     *
     * Initially set to 12
     */
    int m_precision;

    //! default minimimum field width
    /*!
     *  Initially, this is set to 9
     */
    int m_wMin;

    //! Default maximum field width
    /*!
     *  Initially this is set to 19
     */
    int m_wMax;

    //! Local Cropping Control
    CROP_TYPE m_cropCntrl;
};
}

#endif
