#include "Python.h"
#include <string>
#include "../../src/logger.h"

using namespace std;

static std::string ss = "print \"\"\" ";

namespace Cantera {

    class Py_Logger : public Logger {
    public:
        Py_Logger() { cout << "created Py_Logger" << endl;}
        virtual ~Py_Logger() {}

        virtual void write(const string& s) {
            cout << "write..." << endl;
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

        virtual void error(const std::string& msg) {
            string err = "raise \""+msg+"\"";
            PyRun_SimpleString((char *)err.c_str());
        }

        virtual int env() {
            return 2;
        }
    };
}

