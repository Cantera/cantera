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
#include <string.h>
#include <stdlib.h>

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

    static int firstChar(string s) {
        int i;
        int n = static_cast<int>(s.size());
        for (i = 0; i < n; i++) 
            if (s[i] != ' ' && isprint(s[i])) break;
        return i;
    }

    static int lastChar(string s) {
        int i;
        int n = static_cast<int>(s.size());
        for (i = n-1; i >= 0; i--) 
            if (s[i] != ' ' && isprint(s[i])) break;
        return i;
    }

    /** 
     * Strip leading and trailing white space.
     */
    string stripws(string s) {
        int ifirst = firstChar(s);
        int ilast = lastChar(s);
        return s.substr(ifirst, ilast - ifirst + 1); 
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
                if (icolon != string::npos) {
                    name = s.substr(0, icolon);
                    if (iend != string::npos) {
                        num = s.substr(icolon+1, iend-icolon);
                        s = s.substr(iend+1, s.size());
                    }
                    else {
                        num = s.substr(icolon+1, s.size());
                        s = "";
                    }
                    nm = stripws(name);
                    if (x.find(nm) == x.end()) {
                        //if (x[nm] == 0.0) {
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


            
    /**
     * Parse a composition string.
     */
    void split(const string ss, vector<string>& w) {
        string s = ss;
        string::size_type ibegin, iend;
        string name, num, nm;
        do {
            ibegin = s.find_first_not_of(", ;\n\t");
            if (ibegin != string::npos) {
                s = s.substr(ibegin,s.size());
                iend = s.find_first_of(", ;\n\t");
                if (iend != string::npos) {
                    w.push_back(s.substr(0, iend));
                    s = s.substr(iend+1, s.size());
                }
                else {
                    w.push_back(s.substr(0, s.size()));
                    return;
                }
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

    string wrapString(const string& s, int len) {
        int nc = s.size();
        int n, count=0;
        string r;
        for (n = 0; n < nc; n++) {
            if (s[n] == '\n') count = 0;
            else count++;
            if (count > len && s[n] == ' ') {
                r += "\n     ";
                count = 0;
            }
            r += s[n];
        }
        return r;
    }


    /*
     *     This routine strips off blanks and tabs (only leading and trailing
     *     characters) in 'str'.  On return, it returns the number of
     *     characters still included in the string (excluding the null character).
     *
     *      Comments are excluded -> All instances of the comment character, '!',
     *                               are replaced by '\0' thereby terminating
     *                               the string
     *
     *     Parameter list:
     *
     *     str == On output 'str' contains the same characters as on
     *               input except the leading and trailing white space and
     *               comments have been removed.
     */
    int stripLTWScstring(char str[]) {
	int  i = 0, j = 0;
	char ch;
	const char COM_CHAR='\0';
	/*
	 *    Quick Returns
	 */
	if ((str == 0) || (str[0] == '\0')) return (0);

	/* Find first non-space character character */
	while(((ch = str[i]) != '\0') && isspace(ch)) i++;
	
	/*
	 * Move real part of str to the front by copying the string
	 *   - Comments are handled here, by terminating the copy at the
	 *     first comment indicator, and inserting the null character at
	 *     that point.
	 */
 
	while ( (ch = str[j+i]) != '\0' &&
		(ch != COM_CHAR)) {
	  str[j] = ch;
	  j++;
	}
	str[j] = '\0';
	j--;
	/* Remove trailing white space by inserting a null character */    
	while( (j != -1 ) && isspace(str[j])) j--;
	j++;
	str[j] = '\0';
	return (j);
    }


    /*
     * atofCheck is a wrapper around the C stdlib routine atof().
     * It does quite a bit more error checking than atof() or
     * strtod(), and is quite a bit more restrictive.
     *
     *   First it interprets both E, e, d, and D as exponents.
     *   atof() only interprets e or E as an exponent character.
     *
     *   It only accepts a string as well formed if it consists as a 
     *   single token. Multiple words will produce an error message
     *
     *   It will produce an error for NAN and inf entries as well,
     *   in contrast to atof() or strtod().
     *   The user needs to know that a serious numerical issue
     *   has occurred.
     *
     *   It does not accept hexadecimal numbers.
     *
     * On any error, it will throw a CanteraError signal.
     */
    double atofCheck(const char *dptr) {
	if (!dptr) {
	  throw CanteraError("atofCheck", "null pointer to string");
	}
	char *eptr = (char *) malloc(strlen(dptr)+1);
	strcpy(eptr, dptr);
	int ll = stripLTWScstring(eptr);
	if (ll == 0) {
	  throw CanteraError("atofCheck", "string has zero length");
	}
	int numDot = 0;
	int numExp = 0;
	char ch;
	int istart = 0;
	ch = eptr[0];
	if (ch == '+' || ch == '-') {
	  istart = 1;
	}
	for (int i = istart; i < ll; i++) {
	  ch = eptr[i];
	  if (isdigit(ch)) {
	  } else if (ch == '.') {
	    numDot++;
	    if (numDot > 1) {
	      free(eptr);
	      throw CanteraError("atofCheck",
				 "string has more than one .");
	    }
	  } else if (ch == 'e' || ch == 'E' || ch == 'd' || ch == 'D') {
	    numExp++;
	    eptr[i] = 'E';
	    if (numExp > 1) {
	      free(eptr);
	      throw CanteraError("atofCheck", 
				 "string has more than one exp char");
	    }
	    ch = eptr[i+1];
	    if (ch == '+' || ch == '-') {
	      i++;
	    }
	  } else {
	    string hh(dptr);
	    free(eptr);
	    throw CanteraError("aotofCheck",
			       "Trouble processing string, " + hh);
	  }
	}
	double rval = atof(eptr);
	free(eptr);
	return rval;
    }

}
