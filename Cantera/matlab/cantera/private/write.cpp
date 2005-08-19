
#include "mex.h"
#include <string>

static std::string ss = "disp('";

namespace Cantera {

    void writelog(const std::string& s) {
      // debug
      mexEvalString(s.c_str());
      // end debug
        char ch = s[0];
        int n = 0;
        while (ch != '\0') {
            if (ch =='\n') {
                ss += " ');";
                mexEvalString(ss.c_str());
                ss = "disp('doda: ";
            }
            else 
                ss += ch;
            if (ch == '\'') ss += ch;
            n++;
            ch = s[n];
        }
    }

    void error(const std::string& msg) {
        std::string err = "error("+msg+");";
        mexEvalString(err.c_str());
    }

    int userInterface() {
        return 1;
    }
    
}
