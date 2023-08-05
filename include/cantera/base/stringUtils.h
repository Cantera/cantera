/**
 *  @file stringUtils.h Contains declarations for string manipulation
 *       functions within Cantera.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STRINGUTILS_H
#define CT_STRINGUTILS_H

#include "ct_defs.h"

namespace Cantera
{

//! @defgroup StringConversion String Conversion
//! Utility functions for conversion of formatted strings.
//! @ingroup globalUtilFuncs
//! @{

//! Convert a vector to a string (separated by commas)
/*!
 * @param v     vector to be converted
 * @param fmt   Format to be used (printf style) for each element
 * @param sep   Separator
 */
string vec2str(const vector<double>& v, const string& fmt="%g", const string& sep=", ");

//! Strip non-printing characters wherever they are
/*!
 * @param s        Input string
 * @returns a copy of the string, stripped of all non- printing characters.
 */
string stripnonprint(const string& s);

//! Parse a composition string into a map consisting of individual
//! key:composition pairs.
/*!
 * Elements present in *names* but not in the composition string will have
 * a value of 0. Elements present in the composition string but not in *names*
 * will generate an exception. The composition is a double. Example:
 *
 * Input is
 *
 *    "ice:1   snow:2"
 *    names = ["fire", "ice", "snow"]
 *
 * Output is
 *             x["fire"] = 0
 *             x["ice"]  = 1
 *             x["snow"] = 2
 *
 * @param ss    original string consisting of multiple key:composition
 *               pairs on multiple lines
 * @param names (optional) valid names for elements in the composition map. If
 *      empty or unspecified, all values are allowed.
 * @return     map of names to values
 */
Composition parseCompString(const string& ss,
                            const vector<string>& names=vector<string>());

//! Translate a string into one double value
/*!
 * No error checking is done on the conversion.
 *
 * @param val   String value of the double
 * @returns     a double
 */
double fpValue(const string& val);

//! Translate a string into one double value, with error checking
/*!
 * fpValueCheck is a wrapper around the C++ stringstream double parser. It
 * does quite a bit more error checking than atof() or strtod(), and is quite
 * a bit more restrictive.
 *
 * First it interprets both E, e, d, and D as exponents. stringstreams only
 * interpret e or E as an exponent character.
 *
 * It only accepts a string as well formed if it consists as a single token.
 * Multiple words will raise an exception. It will raise a CanteraError for
 * NAN and inf entries as well, in contrast to atof() or strtod(). The user
 * needs to know that a serious numerical issue has occurred.
 *
 * It does not accept hexadecimal numbers.
 *
 * It always use the C locale, regardless of any locale settings.
 *
 * @param val   String representation of the number
 * @returns     a double
 */
double fpValueCheck(const string& val);

//! This function separates a string up into tokens according to the location of
//! white space.
/*!
 * White space includes the new line character. tokens are stripped of leading
 * and trailing white space.
 *
 * The separate tokens are returned in a string vector, v.
 *
 * @param oval   String to be broken up
 * @param v     Output vector of tokens.
 */
void tokenizeString(const string& oval, vector<string>& v);

//! This function separates a string up into tokens according to the location of
//! path separators.
/*!
 * The separate tokens are returned in a string vector, v.
 *
 * @param oval   String to be broken up
 * @param v     Output vector of tokens.
 *
 * @since New in %Cantera 3.0.
 */
void tokenizePath(const string& oval, vector<string>& v);

//! Copy the contents of a string into a char array of a given length
/*!
 *  If *length* is less than the size of *source*, the string will be truncated
 *  and the function will return the length of the buffer required to hold
 *  *source*. Otherwise, returns 0.
 */
size_t copyString(const string& source, char* dest, size_t length);

//! Trim.
/*!
 *  Remove all leading and trailing spaces (with default locale).
 */
string trimCopy(const string &input);

//! Convert to lower case.
/*!
 *  Convert the given string to lower case (with default locale).
 */
string toLowerCopy(const string& input);

//! Case insensitive equality predicate.
/*!
 *  Returns true if and only if all elements in both strings are the same
 *  when compared case insensitively (with default locale).
 */
bool caseInsensitiveEquals(const string &input, const string &test);

//! @}

}

#endif
