#ifndef MLLOGGER_H
#define MLLOGGER_H

#include "mex.h"
#include <string>
#include "cantera/kernel/logger.h"

#include <iostream>
using namespace std;

static string ss = "disp(' ";

namespace Cantera {

    class ML_Logger : public Logger {
    public:
        ML_Logger() {}
        virtual ~ML_Logger() {}


        virtual void write(const string& s) {
            char ch = s[0];
            int n = 0;
            while (ch != '\0') {
                if (ch =='\n') {
                    ss += "');";

                    mexEvalString(ss.c_str());
                    ss = "disp(' ";
                }
                else 
                    ss += ch;
                if (ch == '\'') ss += ch;
                n++;
                ch = s[n];
            }
        }


        virtual void error(const string& msg) {
            string err = "error("+msg+");";
            mexEvalString(err.c_str());
        }

        virtual int env() {
            return 1;
        }
    };
    
}

#endif
