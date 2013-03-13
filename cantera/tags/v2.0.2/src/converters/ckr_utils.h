/**
 *  @file ckr_utils.h
 *
 */

// Copyright 2001  California Institute of Technology

#ifndef CKR_UTILS_H
#define CKR_UTILS_H

#include <math.h>
#include <string>
#include <map>
#include <vector>

namespace ckr
{

/**
 *  Fill vector 'keys' with the keys of map 'mp'
 */
template<class K, class V>
void getMapKeys(const std::map<K,V>& mp, std::vector<K>& keys)
{
    keys.clear();
    typename std::map<K,V>::const_iterator i = mp.begin();
    for (; i != mp.end(); ++i) {
        keys.push_back(i->first);
    }
}


/**
 *  Fill vector 'values' with the values of map 'mp'
 */
template<class K, class V>
void getMapValues(const std::map<K,V>& mp, std::vector<V>& values)
{
    values.clear();
    typename std::map<K,V>::const_iterator i = mp.begin();
    for (; i != mp.end(); ++i) {
        values.push_back(i->second);
    }
}


/**
 *
 *  Template to compare two objects a and b, possibly of different
 *  types, and return the greater of the two, converted to the type of
 *  a. The '<' operator is used for the comparison, and must be
 *  defined for the two types in question.
 *
 */

template<class T, class S>
inline T max(T a, S b)
{
    return (a < b ? b : a);
}



/**
 *
 *  Template to compare two objects a and b, possibly of different
 *  types, and * return the lesser of the two, converted to the type
 *  of a. The '<' operator is used for the comparison, and must be
 *  defined for the two types in question.
 *
 */

template<class T, class S>
inline T min(T a, S b)
{
    return (a < b ? a : b);
}


/**
 *  Template to return a string equal to s, but padded with spaces on
 *  the right as necessary to make the length n.
 */
template<class S>
inline S pad(const S& s, size_t n)
{
    S output;
    output.resize(max(n,s.size()),' ');
    copy(s.begin(), s.end(), output.begin());
    return output;
}

/// Absolute value.
template<class T>
inline T absval(T x)
{
    if (x < 0) {
        return -x;
    }
    return x;
}

/**
 *
 * Iterate through a list of objects that have a numeric member named
 * 'valid', and return false if for any object this attribute is not
 * greater than 0.  Otherwise return true.
 *
 */

template<class L>
inline bool valid(L& list)
{
    size_t i;
    for (i=0; i < list.size(); i++) if (list[i].valid <= 0) {
            return false;
        }
    return true;
}


/// Remove all white space from string s.
void removeWhiteSpace(std::string& s);

void getTokens(std::string& begin,
               int n, std::vector<std::string>& toks, char delim=' ');


/**
 *  Perform a case-insensitive comparison of the first n2 characters
 *  of strings s1 and s2, where n2 is the length of s2. Typically, s1
 *  is an unknown string and s2 is the significant portion of a
 *  keyword. Returns true if a match is found, false otherwise. An asterisk
 *  in string s2 matches any character at that position.
 *
 *  Example: if s1 = "elements", then match(s1, "ELEM") would return true.
 */
bool match(const std::string& s1, const std::string& s2);

/**
 * Check whether string 'word' begins with a Chemkin keyword.
 */
inline bool isKeyword(std::string word)
{
    return (match(word, "ELEM") ||
            match(word, "SPEC") ||
            match(word, "THERM") ||
            match(word, "REAC") ||
            match(word, "END"));
}


bool extractSlashData(std::string& s, std::string& name, std::string& data);
std::string capitalize(const std::string& word);

}

#endif





