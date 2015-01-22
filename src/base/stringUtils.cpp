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
#include "cantera/base/ctml.h"

#include <sstream>
#include <cstdio>

namespace Cantera
{

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

std::string int2str(const size_t n)
{
    std::stringstream ss;
    ss << n;
    return ss.str();
}

std::string vec2str(const vector_fp& v, const std::string& fmt,
                    const std::string& sep)
{
    char buf[64];
    std::stringstream o;
    for (size_t i = 0; i < v.size(); i++) {
        SNPRINTF(buf, 63, fmt.c_str(), v[i]);
        o << v[i];
        if (i != v.size() - 1) {
            o << sep;
        }
    }
    return o.str();
}


std::string lowercase(const std::string& s)
{
    std::string lc(s);
    for (size_t i = 0; i < s.size(); i++) {
        lc[i] = (char) tolower(s[i]);
    }
    return lc;
}

//! Return the position of the first printable character in the string
/*!
 *    @param  s    input string
 *    @return      Returns an int representing the first printable string. If
 *                 none returns the size of the string.
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

//! Return the position of the last printable character in the string
/*!
 *    @param  s    input string
 *    @return      Returns an int representing the first printable string. If
 *                 none returns -1.
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

std::string stripws(const std::string& s)
{
    int ifirst = firstChar(s);
    int ilast = lastChar(s);
    return s.substr(ifirst, ilast - ifirst + 1);
}

std::string stripnonprint(const std::string& s)
{
    std::string ss = "";
    for (size_t i = 0; i < s.size(); i++) {
        if (isprint(s[i])) {
            ss += s[i];
        }
    }
    return ss;
}

compositionMap parseCompString(const std::string& ss,
                               const std::vector<std::string>& names)
{
    compositionMap x;
    for (size_t k = 0; k < names.size(); k++) {
        x[names[k]] = 0.0;
    }

    size_t start = 0;
    size_t stop = 0;
    while (stop < ss.size()) {
        size_t colon = ss.find(':', start);
        if (colon == npos) {
            break;
        }
        size_t valstart = ss.find_first_not_of(" \t\n", colon+1);
        stop = ss.find_first_of(", ;\n\t", valstart);
        std::string name = stripws(ss.substr(start, colon-start));
        if (!names.empty() && x.find(name) == x.end()) {
            throw CanteraError("parseCompString",
                "unknown species '" + name + "'");
        }
        x[name] = fpValueCheck(ss.substr(valstart, stop-colon-1));
        start = ss.find_first_not_of(", ;\n\t", stop+1);
    }
    if (stop != npos && !stripws(ss.substr(stop)).empty()) {
        throw CanteraError("parseCompString", "Found non-key:value data "
            "in composition string: '" + ss.substr(stop) + "'");
    }
    return x;
}

void split(const std::string& ss, std::vector<std::string>& w)
{
    warn_deprecated("split", "To be removed after Cantera 2.2.");
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

int fillArrayFromString(const std::string& str,
                        doublereal* const a, const char delim)
{
    warn_deprecated("fillArrayFromString", "To be removed after Cantera 2.2.");
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
        a[count] = fpValueCheck(num);
        count++;
    }
    return count;
}

std::string getBaseName(const std::string& path)
{
    warn_deprecated("getBaseName", "To be removed after Cantera 2.2.");
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

int intValue(const std::string& val)
{
    return std::atoi(stripws(val).c_str());
}

doublereal fpValue(const std::string& val)
{
    doublereal rval;
    std::stringstream ss(val);
    ss.imbue(std::locale("C"));
    ss >> rval;
    return rval;
}

doublereal fpValueCheck(const std::string& val)
{
    std::string str = stripws(val);
    if (str.empty()) {
        throw CanteraError("fpValueCheck", "string has zero length");
    }
    int numDot = 0;
    int numExp = 0;
    char ch;
    int istart = 0;
    ch = str[0];
    if (ch == '+' || ch == '-') {
        istart = 1;
    }
    for (size_t i = istart; i < str.size(); i++) {
        ch = str[i];
        if (isdigit(ch)) {
        } else if (ch == '.') {
            numDot++;
            if (numDot > 1) {
                throw CanteraError("fpValueCheck",
                                   "string has more than one .");
            }
            if (numExp > 0) {
                throw CanteraError("fpValueCheck",
                                   "string has decimal point in exponent");
            }
        } else if (ch == 'e' || ch == 'E' || ch == 'd' || ch == 'D') {
            numExp++;
            str[i] = 'E';
            if (numExp > 1) {
                throw CanteraError("fpValueCheck",
                                   "string has more than one exp char");
            }
            ch = str[i+1];
            if (ch == '+' || ch == '-') {
                i++;
            }
        } else {
            throw CanteraError("fpValueCheck",
                               "Trouble processing string, " + str);
        }
    }
    return fpValue(str);
}

std::string logfileName(const std::string& infile)
{
    warn_deprecated("logfileName", "To be removed after Cantera 2.2.");
    std::string logfile = getBaseName(infile);
    logfile += ".log";
    return logfile;
}

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
    doublereal val = fpValueCheck(v[0]);
    return val * fp;
}

//!  Find the first white space in a string
/*!
 *   Returns the location of the first white space character in a string
 *
 *   @param   val    Input string to be parsed
 *   @return  In a size_type variable, return the location of the first white
 *             space character. Return npos if none is found
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

//!  Find the first non-white space in a string
/*!
 *   Returns the location of the first non-white space character in a string
 *
 *   @param   val    Input string to be parsed
 *   @return  In a size_type variable, return the location of the first
 *             nonwhite space character. Return npos if none is found
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

void copyString(const std::string& source, char* dest, size_t length)
{
    const char* c_src = source.c_str();
    size_t N = std::min(length, source.length()+1);
    std::copy(c_src, c_src + N, dest);
    if (length != 0) {
        dest[length-1] = '\0';
    }
}

}
