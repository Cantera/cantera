#include "Python.h"
#include <string>
using namespace std;

static std::string ss = "print \"\"\" ";

namespace Cantera {

    void writelog(const string& s) {
        char ch = s[0];
        int n = 0;
        while (ch != '\0') {
            if (ch =='\n') {
                ss += "\"\"\"";
                PyRun_SimpleString((char *)ss.c_str());            
                ss = "print \"\"\"";
            }
            else 
                ss += ch;
            n++;
            ch = s[n];
        }
    }

    void error(const std::string& msg) {
        string err = "raise \""+msg+"\"";
        PyRun_SimpleString((char *)err.c_str());
    }

    int userInterface() {
        return 2;
    }
}
