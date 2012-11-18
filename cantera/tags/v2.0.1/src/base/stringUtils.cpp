/**
 *  @file stringUtils.cpp
 * Contains definitions for string manipulation functions
 *       within Cantera.
 */
// Copyright 2001  California Institute of Technology

//@{
#include "cantera/base/ct_defs.h"

#ifdef _MSC_VER
#define SNPRINTF _snprintf
#else
#define SNPRINTF snprintf
#endif
//@}

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"
#include "cantera/base/ctml.h"

#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cctype>

namespace Cantera
{

//================================================================================================
// Convert a double into a c++ string
/*
 *  This routine doesn't assume a formatting. You
 *  must supply the formatting
 *
 * @param x double to be converted
 * @param fmt   Format to be used (printf style)
 */
std::string fp2str(const double x, const std::string& fmt)
{
    char buf[64];
    int n = SNPRINTF(buf, 63, fmt.c_str(), x);
    if (n > 0) {
        buf[63] = '\0';
        return std::string(buf);
    }
    return std::string(" ");
}
std::string fp2str(const double x)
{
    char buf[64];
    int n = SNPRINTF(buf, 64, "%g" , x);
    if (n > 0) {
        buf[29] = '\0';
        return std::string(buf);
    }
    return std::string(" ");
}
//================================================================================================
/*
 * Convert an integer number to a std::string using sprintf.
 */
std::string int2str(const int n, const std::string& fmt)
{
    char buf[30];
    int m = SNPRINTF(buf, 30, fmt.c_str(), n);
    if (m > 0) {
        buf[29] = '\0';
        return std::string(buf);
    }
    return std::string(" ");
}
//================================================================================================
//  Convert an int to a string
/*
 *  @param n          int to be converted
 */
std::string int2str(const int n)
{
    char buf[30];
    int m = SNPRINTF(buf, 30, "%d", n);
    if (m > 0) {
        buf[29] = '\0';
        return std::string(buf);
    }
    return std::string(" ");
}
//================================================================================================
//  Convert an int to a string
/*
 *  @param n          int to be converted
 */
std::string int2str(const size_t n)
{
    std::stringstream ss;
    ss << n;
    return ss.str();
}
//================================================================================================
std::string lowercase(const std::string& s)
{
    int n = static_cast<int>(s.size());
    std::string lc(s);
    for (int i = 0; i < n; i++) {
        lc[i] = (char) tolower(s[i]);
    }
    return lc;
}
//================================================================================================
//! Return the position of the first printable
//! character in the string
/*!
 *    @param  s    input string
 *    @return      Returns an int representing the first
 *                 printable string. If none returns the
 *                 size of the string.
 */
static int firstChar(const std::string& s)
{
    int i;
    int n = static_cast<int>(s.size());
    for (i = 0; i < n; i++) {
        if (s[i] != ' ' && isprint(s[i])) {
            break;
        }
    }
    return i;
}
//================================================================================================
//! Return the position of the last printable
//! character in the string
/*!
 *    @param  s    input string
 *    @return      Returns an int representing the first
 *                 printable string. If none returns
 *                 -1.
 */
static int lastChar(const std::string& s)
{
    int i;
    int n = static_cast<int>(s.size());
    for (i = n-1; i >= 0; i--)
        if (s[i] != ' ' && isprint(s[i])) {
            break;
        }
    return i;
}
//================================================================================================
// Strip the leading and trailing white space
// from a string
/*
 *  The command isprint() is used to determine printable
 *  characters.
 *
 *    @param   s       Input string
 *    @return  Returns a copy of the string, stripped
 *             of leading and trailing white space
 */
std::string stripws(const std::string& s)
{
    int ifirst = firstChar(s);
    int ilast = lastChar(s);
    return s.substr(ifirst, ilast - ifirst + 1);
}
//================================================================================================
// Strip non-printing characters wherever they are
/*
 *   @param s        Input string
 *   @return         Returns a copy of the string,
 *                   stripped of all non-printing characters.
 */
std::string stripnonprint(const std::string& s)
{
    int i;
    int n = static_cast<int>(s.size());
    std::string ss = "";
    for (i = 0; i < n; i++) {
        if (isprint(s[i])) {
            ss += s[i];
        }
    }
    return ss;
}
//================================================================================================
// Parse a composition string into a map consisting of individual key:composition
// pairs.
/*
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
void parseCompString(const std::string& ss, Cantera::compositionMap& x)
{
    std::string s = ss;
    std::string::size_type icolon, ibegin, iend;
    std::string name, num, nm;
    do {
        ibegin = s.find_first_not_of(", ;\n\t");
        if (ibegin != std::string::npos) {
            s = s.substr(ibegin,s.size());
            icolon = s.find(':');
            iend = s.find_first_of(", ;\n\t");
            //icomma = s.find(',');
            if (icolon != std::string::npos) {
                name = s.substr(0, icolon);
                if (iend != std::string::npos) {
                    num = s.substr(icolon+1, iend-icolon);
                    s = s.substr(iend+1, s.size());
                } else {
                    num = s.substr(icolon+1, s.size());
                    s = "";
                }
                nm = stripws(name);
                if (x.find(nm) == x.end()) {
                    throw CanteraError("parseCompString",
                                       "unknown species " + nm);
                }
                x[nm] = atof(num.c_str());
            } else {
                s = "";
            }
        }
    } while (s != "");
}
//================================================================================================
//   Parse a composition string into individual key:composition
//   pairs
/*
 *
 *     @param ss   original string consisting of multiple key:composition
 *                 pairs on multiple lines
 *     @param w    Output vector consisting of single key:composition
 *                 items in each index.
 */
void split(const std::string& ss, std::vector<std::string>& w)
{
    std::string s = ss;
    std::string::size_type ibegin, iend;
    std::string name, num, nm;
    do {
        ibegin = s.find_first_not_of(", ;\n\t");
        if (ibegin != std::string::npos) {
            s = s.substr(ibegin,s.size());
            iend = s.find_first_of(", ;\n\t");
            if (iend != std::string::npos) {
                w.push_back(s.substr(0, iend));
                s = s.substr(iend+1, s.size());
            } else {
                w.push_back(s.substr(0, s.size()));
                return;
            }
        }
    } while (s != "");
}
//================================================================================================
int fillArrayFromString(const std::string& str,
                        doublereal* const a, const char delim)
{
    std::string::size_type iloc;
    int count = 0;
    std::string num;
    std::string s = str;
    while (s.size() > 0) {
        iloc = s.find(delim);
        if (iloc > 0) {
            num = s.substr(0, iloc);
            s = s.substr(iloc+1,s.size());
        } else {
            num = s;
            s = "";
        }
        a[count] = atofCheck(num.c_str());
        count++;
    }
    return count;
}
//================================================================================================
// Get the file name without the path or extension
/*
 *   @param fullPath   Input file name consisting
 *                     of the full file name
 *
 *  @return Returns the basename
 */
std::string getBaseName(const std::string& path)
{
    std::string file;
    size_t idot = path.find_last_of('.');
    size_t islash = path.find_last_of('/');
    if (idot > 0 && idot < path.size()) {
        if (islash > 0 && islash < idot) {
            file = path.substr(islash+1, idot-islash-1);
        } else {
            file = path.substr(0,idot);
        }
    } else {
        file = path;
    }
    return file;
}
//================================================================================================
int intValue(std::string val)
{
    return std::atoi(stripws(val).c_str());
}
//================================================================================================
doublereal fpValue(std::string val)
{
    return std::atof(stripws(val).c_str());
}
//================================================================================================
doublereal fpValueCheck(std::string val)
{
    return atofCheck(stripws(val).c_str());
}
//================================================================================================
//  Generate a logfile name based on an input file name
/*
 *   It tries to find the basename. Then, it appends a .log
 *   to it.
 *
 *   @param infile      Input file name
 *
 *  @return Returns a logfile name
 */
std::string logfileName(const std::string& infile)
{
    std::string logfile = getBaseName(infile);
    logfile += ".log";
    return logfile;
}
//================================================================================================
//    Line wrap a string via a copy operation
/*
 *   @param s   Input string to be line wrapped
 *   @paramlen  Length at which to wrap. The
 *              default is 70.
 */
std::string wrapString(const std::string& s, const int len)
{
    int count=0;
    std::string r;
    for (size_t n = 0; n < s.size(); n++) {
        if (s[n] == '\n') {
            count = 0;
        } else {
            count++;
        }
        if (count > len && s[n] == ' ') {
            r += "\n     ";
            count = 0;
        }
        r += s[n];
    }
    return r;
}
//================================================================================================
// Parse a name string, separating out the phase name from the species name
/*
 *   Name strings must not contain these internal characters "; \n \t "
 *   Only one colon is allowed, the one separating the phase name from the
 *   species name. Therefore, names may not include a colon.
 *
 *   @param nameStr   (input) Name string containing the phase name and the species
 *                            name separated by a colon. The phase name is optional.
 *                             example:   "silane:SiH4"
 *   @param phaseName (output) Name of the phase, if specified. If not specified,
 *                             a blank string is returned.
 *   @return          (output) Species name is returned. If nameStr is blank
 *                             an empty string is returned.
 */
std::string parseSpeciesName(const std::string& nameStr, std::string& phaseName)
{
    std::string s = stripws(nameStr);
    std::string::size_type ibegin, iend, icolon;
    phaseName = "";
    ibegin = s.find_first_not_of(" ;\n\t");
    if (ibegin != std::string::npos) {
        s = s.substr(ibegin,s.size());
        icolon = s.find(':');
        iend = s.find_first_of(" ;\n\t");
        if (icolon != std::string::npos) {
            phaseName = s.substr(0, icolon);
            s = s.substr(icolon+1, s.size());
            icolon =  s.find(':');
            if (icolon != std::string::npos) {
                throw CanteraError("parseSpeciesName()", "two colons in name: " + nameStr);
            }
        }
        if (iend != std::string::npos) {
            throw CanteraError("parseSpeciesName()",
                               "Species name has \" ;/\n/\t\" in the middle of it: " + nameStr);
        }
    }
    return s;
}
//================================================================================================
// Routine strips off white space from a c character string
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
 * @param  str   On output 'str' contains the same characters as on
 *               input except the leading and trailing white space and
 *               comments have been removed.
 */
int stripLTWScstring(char str[])
{
    int  i = 0, j = 0;
    char ch;
    const char COM_CHAR='\0';
    /*
     *    Quick Returns
     */
    if ((str == 0) || (str[0] == '\0')) {
        return (0);
    }

    /* Find first non-space character character */
    while (((ch = str[i]) != '\0') && isspace(ch)) {
        i++;
    }

    /*
     * Move real part of str to the front by copying the string
     *   - Comments are handled here, by terminating the copy at the
     *     first comment indicator, and inserting the null character at
     *     that point.
     */

    while ((ch = str[j+i]) != '\0' &&
            (ch != COM_CHAR)) {
        str[j] = ch;
        j++;
    }
    str[j] = '\0';
    j--;
    /* Remove trailing white space by inserting a null character */
    while ((j != -1) && isspace(str[j])) {
        j--;
    }
    j++;
    str[j] = '\0';
    return (j);
}
//================================================================================================
// Translate a char string into a single double
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
 *  @param dptr  pointer to the input c string
 *  @return      Returns the double
 *
 * On any error, it will throw a CanteraError signal.
 */
doublereal atofCheck(const char* const dptr)
{
    if (!dptr) {
        throw CanteraError("atofCheck", "null pointer to string");
    }
    char* eptr = (char*) malloc(strlen(dptr)+1);
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
            std::string hh(dptr);
            free(eptr);
            throw CanteraError("atofCheck",
                               "Trouble processing string, " + hh);
        }
    }
    doublereal rval = atof(eptr);
    free(eptr);
    return rval;
}
//================================================================================================
// Interpret one or two token string as a single double
/*
 *   This is similar to atof(). However, the second token
 *   is interpreted as an MKS units string and a conversion
 *   factor to MKS is applied.
 *
 *   Example
 *  " 1.0 atm"
 *
 *   results in the number 1.01325e5
 *
 *   @param strSI string to be converted. One or two tokens
 *
 *   @return returns a converted double
 */
doublereal strSItoDbl(const std::string& strSI)
{
    std::vector<std::string> v;
    tokenizeString(strSI, v);
    doublereal fp = 1.0;
    size_t n = v.size();
    if (n > 2 || n < 1) {
        throw CanteraError("strSItoDbl",
                           "number of tokens is too high");
    } else if (n == 2) {
        fp = toSI(v[1]);
    }
    doublereal val = atofCheck(v[0].c_str());
    return (val * fp);
}
//================================================================================================
//!  Find the first white space in a string
/*!
 *   Returns the location of the first white space character in a string
 *
 *   @param   val    Input string to be parsed
 *   @return  In a size_type variable, return the location of the first white space character.
 *             Return npos if none is found
 */
static std::string::size_type findFirstWS(const std::string& val)
{
    std::string::size_type ibegin = std::string::npos;
    int j = 0;
    std::string::const_iterator i = val.begin();
    for (; i != val.end(); i++) {
        char ch = *i;
        int ll = (int) ch;
        if (isspace(ll)) {
            ibegin = (std::string::size_type) j;
            break;
        }
        j++;
    }
    return ibegin;
}
//================================================================================================
//!  Find the first non-white space in a string
/*!
 *   Returns the location of the first non-white space character in a string
 *
 *   @param   val    Input string to be parsed
 *   @return  In a size_type variable, return the location of the first nonwhite space character.
 *             Return npos if none is found
 */
static std::string::size_type findFirstNotOfWS(const std::string& val)
{
    std::string::size_type ibegin = std::string::npos;
    int j = 0;
    std::string::const_iterator i = val.begin();
    for (; i != val.end(); i++) {
        char ch = *i;
        int ll = (int) ch;
        if (!isspace(ll)) {
            ibegin = (std::string::size_type) j;
            break;
        }
        j++;
    }
    return ibegin;
}
//================================================================================================
// This function  separates a string up into tokens
// according to the location of white space.
/*
 *    The separate tokens are returned in a string vector, v.
 *
 *  @param oval   String to be broken up
 *  @param v     Output vector of tokens.
 */
void tokenizeString(const std::string& oval,
                    std::vector<std::string>& v)
{
    std::string val(oval);
    std::string::size_type ibegin, iend;
    v.clear();
    while (1 > 0) {
        ibegin = findFirstNotOfWS(val);
        if (ibegin != std::string::npos) {
            val = val.substr(ibegin,val.size());
            iend = findFirstWS(val);
            if (iend == std::string::npos) {
                v.push_back(val);
                break;
            } else {
                v.push_back(val.substr(0,iend));
                val = val.substr(iend+1,val.size());
            }
        } else {
            break;
        }
    }
}
//================================================================================================

}
