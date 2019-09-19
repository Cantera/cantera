/**
 *  @file stringUtils.h Contains declarations for string manipulation
 *       functions within Cantera.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_STRINGUTILS_H
#define CT_STRINGUTILS_H

#include "ct_defs.h"
#include "cantera/base/fmt.h"

#include <string>

namespace Cantera
{

//! Convert a vector to a string (separated by commas)
/*!
 * @param v     vector to be converted
 * @param fmt   Format to be used (printf style) for each element
 * @param sep   Separator
 */
std::string vec2str(const vector_fp& v, const std::string& fmt="%g",
                    const std::string& sep=", ");

//! Strip non-printing characters wherever they are
/*!
 * @param s        Input string
 * @returns a copy of the string, stripped of all non- printing characters.
 */
std::string stripnonprint(const std::string& s);

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
compositionMap parseCompString(const std::string& ss,
        const std::vector<std::string>& names=std::vector<std::string>());

//! Translate a string into one integer value
/*!
 * No error checking is done on the conversion. The c stdlib function atoi() is
 * used.
 *
 * @param val   String value of the integer
 * @returns     an integer
 */
int intValue(const std::string& val);

//! Translate a string into one doublereal value
/*!
 * No error checking is done on the conversion.
 *
 * @param val   String value of the double
 * @returns     a double
 */
doublereal fpValue(const std::string& val);

//! Translate a string into one doublereal value, with error checking
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
doublereal fpValueCheck(const std::string& val);

//! Parse a name string, separating out the phase name from the species name
/*!
 * Name strings must not contain these internal characters "; \n \t ," Only one
 * colon is allowed, the one separating the phase name from the species name.
 * Therefore, names may not include a colon.
 *
 * @param[in] nameStr    Name string containing the phase name and the species
 *                       name separated by a colon. The phase name is optional.
 *                       example: "silane:SiH4"
 * @param[out] phaseName Name of the phase, if specified. If not specified, a
 *                       blank string is returned.
 * @returns species name. If nameStr is blank an empty string is returned.
 */
std::string parseSpeciesName(const std::string& nameStr, std::string& phaseName);

//! Interpret one or two token string as a single double
/*!
 * This is similar to atof(). However, the second token is interpreted as an
 * MKS units string and a conversion factor to MKS is applied.
 *
 * Example: "1.0 atm" results in the number 1.01325e5.
 *
 * @param strSI string to be converted. One or two tokens
 * @returns a converted double
 */
doublereal strSItoDbl(const std::string& strSI);

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
void tokenizeString(const std::string& oval,
                    std::vector<std::string>& v);

//! Copy the contents of a std::string into a char array of a given length
/*!
 *  If *length* is less than the size of *source*, the string will be truncated
 *  and the function will return the length of the buffer required to hold
 *  *source*. Otherwise, returns 0.
 */
size_t copyString(const std::string& source, char* dest, size_t length);

//! Trim.
/*!
 *  Remove all leading and trailing spaces (with default locale).
 */
std::string trimCopy(const std::string &input);

//! Convert to lower case.
/*!
 *  Convert the given string to lower case (with default locale).
 */
std::string toLowerCopy(const std::string& input);

//! Case insensitive equality predicate.
/*!
 *  Returns true if and only if all elements in both strings are the same
 *  when compared case insensitively (with default locale).
 */
bool caseInsensitiveEquals(const std::string &input, const std::string &test);

}

#endif
