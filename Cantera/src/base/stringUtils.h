/**
 *  @file std::stringUtils.h
 *       Contains declarations for string manipulation functions
 *       within Cantera.
 */

/*
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_STRINGUTILS_H
#define CT_STRINGUTILS_H

#include "ct_defs.h"

#include <string>

namespace Cantera {

  class Phase;
  class ThermoPhase;

  //! Convert a double into a c++ string
  /*!
   *  This routine doesn't assume a formatting. You
   *  must supply the formatting
   *
   * @param x double to be converted
   * @param fmt   Format to be used (printf style)
   */
  std::string fp2str(const double x, const std::string &fmt);

  //! Convert a double into a c++ string
  /*!
   * The default format to use is equivalent to the default
   * format used by printf's %g formatting.
   *
   * @param x double to be converted
   */
  std::string fp2str(const double x);

  //!  Convert an int to a string using a format converter
  /*!
   *  @param n          int to be converted
   *  @param fmt        format converter for an int int the printf command
   */
  std::string int2str(const int n, const std::string &fmt);

  //!  Convert an int to a string 
  /*!
   *  @param n          int to be converted
   */
  std::string int2str(const int n);

  //! Strip the leading and trailing white space
  //! from a string
  /*!
   *  The command isprint() is used to determine printable
   *  characters.
   *
   *    @param   s       Input string
   *    @return  Returns a copy of the string, stripped
   *             of leading and trailing white space
   */
  std::string stripws(const std::string &s);

  //! Strip non-printing characters wherever they are
  /*!
   *   @param s        Input string
   *   @return         Returns a copy of the string,
   *                   stripped of all non-printing characters.
   */
  std::string stripnonprint(const std::string &s);

  //! Cast a copy of a string to lower case
  /*!
   *   @param s        Input string
   *   @return         Returns a copy of the string,
   *                   with all characters lowercase.
   */
  std::string lowercase(const std::string &s);

  //! Parse a composition string into a map consisting of individual key:composition
  //! pairs.
  /*!
   *  The composition is a double.
   * Example
   *
   *  Input is
   *
   *    "fire:0   ice:1   snow:2"
   *
   *  Output is
   *             x["fire"] = 0
   *             x["ice"]  = 1
   *             x["snow"] = 2
   *
   *     @param ss   original string consisting of multiple key:composition
   *                 pairs on multiple lines
   *     @param x    Output map consisting of a composition
   *                 map, which is a string to double map
   */
  void parseCompString(const std::string &ss, Cantera::compositionMap& x);


  //! Parse a composition string into individual key:composition
  //! pairs
  /*!
   *
   *     @param ss   original string consisting of multiple key:composition
   *                 pairs on multiple lines
   *     @param w    Output vector consisting of single key:composition
   *                 items in each index.
   */
  void split(const std::string &ss, std::vector<std::string>& w);

  int fillArrayFromString(const std::string& str, doublereal* a, char delim = ' ');
  std::string report(const ThermoPhase& th, bool show_thermo = true);
  std::string formatCompList(const Phase& mix, int xyc);
  std::string logfileName(const std::string& infile);    
  std::string getFileName(const std::string& path);

  //! Translate a string into one integer value
  /*!
   *  No error checking is done on the conversion. The c stdlib function
   *  atoi() is used.
   *
   *  @param val   String value of the integer
   *
   *  @return      Returns an integer
   */
  int intValue(std::string val);

  //! Translate a string into one doublereal value
  /*!
   *  No error checking is done on the conversion. The c stdlib function
   *  atof() is used.
   *
   *  @param val   String value of the double
   *
   *  @return      Returns a doublereal value
   */
  doublereal fpValue(std::string val);

  //! Translate a string into one doublereal value
  /*!
   *  Error checking is carried on the conversion. 
   *
   *  @param val   String value of the double
   *
   *  @return      Returns a doublereal value
   */
  doublereal fpValueCheck(std::string val);

  std::string wrapString(const std::string& s, int len=70);

  int stripLTWScstring(char str[]);
  double atofCheck(const char *dptr);
  doublereal strSItoDbl(const std::string& strSI); 


  //! This function separates a string up into tokens
  //! according to the location of white space.
  /*!
   *  White space includes the new line character. tokens
   *  are stripped of leading and trailing white space.
   *
   *  The separate tokens are returned in a string vector, v.
   *
   *  @param oval   String to be broken up
   *  @param v     Output vector of tokens.
   */
  void tokenizeString(const std::string& oval,
                      std::vector<std::string>& v);

}

#endif
