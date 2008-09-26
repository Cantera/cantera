/**
 *  @file std::stringUtils.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CT_STRINGUTILS_H
#define CT_STRINGUTILS_H

#include "ct_defs.h"
#include <string>

namespace Cantera {

  class Phase;
  class ThermoPhase;

  std::string fp2str(double x, std::string fmt);

  //! Convert a double into a c++ string
  /*!
   * The default format to use is equivalent to the default
   * format used by printf's %g formatting.
   *
   * @param x double to be converted
   */
  std::string fp2str(double x);
  std::string int2str(int n, std::string fmt);
  std::string int2str(int n);
  std::string stripws(std::string s);
  std::string stripnonprint(std::string s);
  std::string lowercase(const std::string &s);
  void parseCompString(const std::string ss, compositionMap& x);
  void split(const std::string ss, std::vector<std::string>& w);
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

  std::string wrapString(const std::string& s, int len=70);

  int stripLTWScstring(char str[]);
  double atofCheck(const char *dptr);
  doublereal strSItoDbl(const std::string& strSI); 

}

#endif
