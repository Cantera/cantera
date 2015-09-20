#include <algorithm>
#include "StringFunct.h"

using std::string;
using std::vector;

namespace Cantera {
namespace String {
//==============================================================================

void tokenize(
    const string &str, vector<string> &tokens, const string &delim,
    const bool multi_delim)
{
    if (str.empty()) return;
    
    string::size_type last_pos  = 0;
    string::size_type delim_pos;
    if (multi_delim)
        delim_pos = str.find_first_of(delim);
    else
        delim_pos = str.find(delim);
    
    tokens.clear();
    
    while (last_pos != string::npos) {
        if (last_pos == delim_pos) {
            // Skip over delimeter(s)
            if (multi_delim) {
                last_pos++;
                delim_pos = str.find_first_of(delim, last_pos);
            } else {
                last_pos += delim.length();
                delim_pos = str.find(delim, last_pos);
            }
            
            // Special end case to prevent adding empty string to tokens
            if (last_pos >= str.length())
                last_pos = string::npos;
        } else {
            // Copy token
            if (delim_pos == string::npos)
                tokens.push_back(str.substr(last_pos));
            else
                tokens.push_back(str.substr(last_pos, delim_pos-last_pos));
            
            last_pos = delim_pos;
        }
    }        
}

//==============================================================================

string& eraseAll(string& str, const string& to_erase)
{
    string::size_type pos = str.find(to_erase);
    
    while (pos != string::npos) {
        str.erase(pos, to_erase.length());
        pos = str.find(to_erase, pos);
    }
    
    return str;
}

//==============================================================================

} //namespace String
} // namespace Cantera
