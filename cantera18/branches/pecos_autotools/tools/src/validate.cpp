/**
 *  @file validate.cpp
 *  validate a CK-format reaction mechanism file
 */

// Copyright 2001  California Institute of Technology

#include <iostream>
#include "Cantera.h"
#include "IdealGasMix.h"
#include <math.h>
using namespace Cantera;

int main(int argc, char** args) {

    try {
        IdealGasMix gas;
        char* infile = 0; 
        char* dbfile = 0;
        bool ok = true;
        if (argc == 1 || argc > 3) { 
            cout << "usage: validate <input_file> <database>" << endl;
            return 0;
        }
        else if (argc == 2) {
            infile = args[1];
            ok = gas.import(string(infile),"",true);
        }
        else if (argc == 3) {
            infile = args[1];
            dbfile = args[2];
            ok = gas.import(string(infile),string(dbfile),true);
        }
        showErrors(cerr);
        if (!ok) 
            cerr << infile << " contains errors" << endl;
        else
            cout << infile << " is valid." << endl;
        return int(!ok);
    }
    catch (...) {
        cerr << "exception raised." << endl;
        showErrors(cerr);
    }
}
