/**
 *  @file stringUtils.cpp
 *
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "ctexceptions.h"
#include <stdio.h>

namespace Cantera {


    /**
     * Convert a floating point number to a string using sprintf.
     */
    string fp2str(double x, string fmt) {
        char buf[30];
        sprintf(buf, fmt.c_str(), x);
        return string(buf);
    }

    /**
     * Convert an integer number to a string using sprintf.
     */
    string int2str(int n, string fmt) {
        char buf[30];
        sprintf(buf, fmt.c_str(), n);
        return string(buf);
    }

    string lowercase(string s) {
        int n = static_cast<int>(s.size());
        string lc(s);
        for (int i = 0; i < n; i++) lc[i] = tolower(s[i]);
        return lc;
    }

    /** 
     * Strip white space.
     */
    string stripws(string s) {
        int i;
        bool sempty = true;
        int n = static_cast<int>(s.size());
        string ss = "";
        for (i = 0; i < n; i++) {
            if (s[i] != ' ' && isprint(s[i])) {
                ss += s[i];
                sempty = false;
            }
            else if (!sempty) break;
        }
        return ss;
    }

    /** 
     * Strip non-printing characters.
     */
    string stripnonprint(string s) {
        int i;
        int n = static_cast<int>(s.size());
        string ss = "";
        for (i = 0; i < n; i++) {
            if (isprint(s[i])) {
                ss += s[i];
            }
        }
        return ss;
    }

            
    /**
     * Parse a composition string.
     */
    void parseCompString(const string ss, compositionMap& x) {
        string s = ss;
		string::size_type icolon, ibegin, iend;
        string name, num, nm;
        do {
            ibegin = s.find_first_not_of(", ;\n\t");
			if (ibegin != string::npos) {
                s = s.substr(ibegin,s.size());
                icolon = s.find(':');
                iend = s.find_first_of(", ;\n\t");
                //icomma = s.find(',');
                if (icolon > 0) {
                    name = s.substr(0, icolon);
                    if (iend > 0) {
                        num = s.substr(icolon+1, iend-icolon);
                        s = s.substr(iend+1, s.size());
                    }
                    else {
                        num = s.substr(icolon+1, s.size());
                        s = "";
                    }
                    nm = stripws(name);
                    if (x[nm] == 0.0) {
		       throw CanteraError("parseCompString",
                                          "unknown species " + nm);
		    }
                    x[nm] = atof(num.c_str());
                }
                else s = "";
            }
        }
        while (s != "");
    }

    int fillArrayFromString(const string& str, doublereal* a, char delim) {
		string::size_type iloc;
        int count = 0;
        string num, s;
        s = str;
        while (s.size() > 0) {
            iloc = s.find(delim);
            if (iloc > 0) {
                num = s.substr(0, iloc);
                s = s.substr(iloc+1,s.size());
            }
            else {
                num = s;
                s = "";
            }
            a[count] = atof(num.c_str());
            count++;
        }
        return count;
    }

    /**
     *  Get the file name without the path or extension
     */
    string getFileName(const string& path) {
        string file;
        size_t idot = path.find_last_of('.');
        size_t islash = path.find_last_of('/');
        if (idot > 0 && idot < path.size()) {
            if (islash > 0 && islash < idot) {
                file = path.substr(islash+1, idot-islash-1);
            }
            else {
                file = path.substr(0,idot);
            }
        }
        else {
            file = path;
        }       
        return file;
    }


    /**
     *  Generate a logfile name based on an input file name
     */
    string logfileName(const string& infile) {
        string logfile;
        size_t idot = infile.find_last_of('.');
        size_t islash = infile.find_last_of('/');
        if (idot > 0 && idot < infile.size()) {
            if (islash > 0 && islash < idot) {
                logfile = infile.substr(islash+1, idot-islash-1) + ".log";
            }
            else {
                logfile = infile.substr(0,idot) + ".log";
            }
        }
        else {
            logfile = infile + ".log";
        }       
        return logfile;
    }
}
