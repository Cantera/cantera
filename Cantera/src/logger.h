#ifndef CT_LOGGER_H
#define CT_LOGGER_H

#include <iostream>
using namespace std;

namespace Cantera {

    class Logger {
    public:

        Logger() {}
        virtual ~Logger() {}

        virtual void write(const string& msg) {
            cout << msg;
        }

        virtual void error(const string& msg) {
            cerr << msg << endl;
            exit(-1);
        }

        virtual int env() { return 0; }
    };

}
#endif
