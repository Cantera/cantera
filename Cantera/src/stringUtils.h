/**
 *  @file stringUtils.h
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

    string fp2str(double x, string fmt = "%g");
    string int2str(int n, string fmt = "%d");
    string stripws(string s);
    string stripnonprint(string s);
    string lowercase(string s);
    void parseCompString(const string ss, compositionMap& x);
    void split(const string ss, vector<string>& w);
    int fillArrayFromString(const string& str, doublereal* a, char delim = ' ');
    string report(const ThermoPhase& th, bool show_thermo = true);
    string formatCompList(const Phase& mix, int xyc);
    string logfileName(const string& infile);    
    string getFileName(const string& path);

    inline int intValue(string val) {
        return atoi(stripws(val).c_str());
    }

    inline doublereal fpValue(string val) {
        return atof(stripws(val).c_str());
    }
    string wrapString(const string& s, int len=70);

    int stripLTWScstring(char str[]);
    double atofCheck(const char *dptr);

}

#endif
