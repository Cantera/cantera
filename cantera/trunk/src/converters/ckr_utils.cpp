/**
 *  @file ckr_utils.cpp
 *
 */

// Copyright 2001  California Institute of Technology

#include <ctype.h>

#include <math.h>
#include "ckr_utils.h"
#include <string.h>

using namespace std;

namespace ckr
{


bool match(const std::string& s1, const std::string& s2)
{
    size_t n = s2.size();
    if (s1.size() < n) {
        return false;
    }
    for (size_t i = 0; i < n; i++)
        if (s2[i] != '*' && (toupper(s1[i]) != toupper(s2[i]))) {
            return false;
        }
    return true;
}


void removeWhiteSpace(std::string& s)
{
    string r;
    int ssize = static_cast<int>(s.size());
    for (int n = 0; n < ssize; n++) if (s[n] != ' '
                                            && s[n] != '\t' && s[n] != '\n') {
            r += s[n];
        }
    s = r;
}

/**
 * Get tokens from char array of length n beginning at 'begin'. Tokens are
 * delimited by character 'delim' (space by default), with the exception of
 * text contained within forward slashes ('/'), which is treated literally.
 *
 * code
 * vector<string> tokens;
 * char line[] = "a b/3.0 txt/  c /4.0 6/ d";
 * int n = strlen(line);
 * getTokens(line, n, tokens);
 * for (int i = 0; i < tokens.size(); i++) cout << tokens[i] << endl;
 * endcode
 *
 */

void getTokens(std::string& s, int n, std::vector<std::string>& toks, char delim)
{
    string::iterator q, p = s.begin(), end = p + n;
    vector<string> tokk;
    int inslash = -1;

    //p = begin;
    while (1 > 0) {
        for (; p < end; p++) {
            if (*p != delim) {
                break;
            }
        }
        q = p;
        for (; q < end; q++) {
            if (*q == '/') {
                inslash *= -1;
            }
            if (inslash < 0 && *q == delim) {
                break;
            }
        }
        if (p != q) {
            tokk.push_back(s.substr(p - s.begin(), q - p));
        }
        p = q;
        if (p == end) {
            break;
        }
    }

    toks.clear();
    int nt = static_cast<int>(tokk.size());
    string t = "";
    for (int i = 0; i < nt; i++) {
        if (tokk[i][0] == '/') {
            t += tokk[i];
        } else {
            if (t != "") {
                toks.push_back(t);
            }
            t = tokk[i];
        }
    }
    if (t != "") {
        toks.push_back(t);
    }
}


/**
 * Look for a slash-delimited number in string s,
 * and if found set v to the numerical value, and
 * set s to the portion of the string before the
 * first slash. Return true if slash data found,
 * false otherwise.
 */
bool extractSlashData(std::string& s, std::string& name, std::string& data)
{
    int slen = static_cast<int>(s.size());
    string::size_type n = s.find_first_of("/");
    if (n != string::npos && (static_cast<int>(n) < slen)) {
        int m;
        for (m = static_cast<int>(n)+1; m < slen; m++) if (s[m] == '/') {
                break;
            }
        if (m < slen) {
            data = s.substr(n+1,m-n-1);
            name = s.substr(0,n);
            removeWhiteSpace(name);
            s = s.substr(m+1,1000);
            return true;
        } else {
            name = s;
            removeWhiteSpace(name);
            data = "";
            s = "";
            return false;
        }
    } else {
        name = s;
        removeWhiteSpace(name);
        data = "";
        s = "";
        return false;
    }
}

/**
 * Return a modified version of string word, in which
 * the first letter is upper case, and the rest are
 * lower case.
 */
string capitalize(const std::string& word)
{
    string cap = word;
    int n = static_cast<int>(word.size());
    if (n > 0) {
        cap[0] = (char) toupper(word[0]);
        for (int m = 1; m < n; m++) {
            cap[m] = (char) tolower(word[m]);
        }
    }
    return cap;
}

}

