/**
 * @file ck2cti.cpp
 *
 * Program to convert CK-format reaction mechanism files to Cantera input format.
 *
 */
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <string>
using namespace std;

#include "ct_defs.h"
#include "global.h"
#include "converters/ck2ct.h"

using namespace Cantera;

int showHelp() {
    cout << "\nck2cti: convert a CK-format reaction mechanism file to Cantera input format.\n"
         << "\n   D. G. Goodwin, Caltech \n"
         << "   Version 1.0, August 2003.\n\n"
         << endl;
    cout << "options:" << endl;
    cout << "    -i  <input file> \n"
         << "    -t  <thermo database> \n"
         << "    -tr <transport database> \n"
         << "    -id  <identifier> \n\n"
         << "The results are written to the standard output.\n";
    return 0;
}

int main(int argc, char** argv) {
    string infile="chem.inp", dbfile="", trfile="", logfile;
    string idtag = "gas";
    int i=1;
    if (argc == 1) return showHelp();
 
    while (i < argc) {
        string arg = string(argv[i]);
        if (i < argc-1) {
            if (arg == "-i") {
                infile = argv[i+1];
                ++i;
            }
            else if (arg == "-t") {
                dbfile = argv[i+1];
                ++i;
            }
            else if (arg == "-tr") {
                trfile = argv[i+1];
                ++i;
            }
            else if (arg == "-id") {
                idtag = argv[i+1];
            }
        }
        else if (arg == "-h" || argc < 3) {
            return showHelp(); 
        }
        ++i;
    }

    int ierr = pip::convert_ck(infile.c_str(), dbfile.c_str(), trfile.c_str(), 
        idtag.c_str());

    if (ierr < 0) {
        showErrors(cerr);
    }
    return ierr;
}
