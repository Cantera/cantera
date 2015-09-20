#ifndef UTILITIES_STRING_UTILS_H
#define UTILITIES_STRING_UTILS_H

#include <string>
#include <vector>

namespace Cantera {
namespace String {
/**
 * Tokenizes a string into a vector of sub strings that are separated by the 
 * characters in teh delim string.
 */
void tokenize(
    const std::string &str, std::vector<std::string> &tokens, 
    const std::string &delim = " ", const bool multi_delim = true);


/**
 * Removes all characters from the string belonging to to_erase.
 */
std::string& eraseAll(
    std::string &str, const std::string &to_erase = " \n\t");

} // namespace String
} // namespace Cantera

#endif // UTILITIES_STRING_UTILS_H
