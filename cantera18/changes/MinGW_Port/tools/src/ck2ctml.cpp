/**
 * @file ck2ctml.cpp
 *
 * Program to convert CK-format reaction mechanism files to CTML format.
 *
 */
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <string>
using namespace std;

#include "converters/ck2ctml.h"
#include "converters/ck2ct.h"

using namespace ctml;

int showHelp() {
    cout << "\nck2ckml: convert a CK-format reaction mechanism file to CTML.\n"
         << "\n   D. G. Goodwin, Caltech \n"
         << "   Version 1.0, August 2002.\n\n"
         << endl;
    cout << "options:" << endl;
    cout << "    -i  <input file> \n"
         << "    -o  <output file> \n"
         << "    -t  <thermo database> \n"
         << "    -tr <transport database> \n"
         << "    -id  <identifier> \n";
    return 0;
}

int main(int argc, char** argv) {
    string infile="chem.inp", dbfile="", trfile="", logfile, outfile="";
    string idtag = "gas";
    //    ckr::CKReader r;
    //r.validate = true;
    int i=1;
    if (argc == 1) return showHelp();
 
    while (i < argc) {
        string arg = string(argv[i]);
        if (i < argc-1) {
            if (arg == "-i") {
                infile = argv[i+1];
                ++i;
            }
            else if (arg == "-o") {
                outfile = argv[i+1];
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

//#define MAKE_CT_INPUT
#ifdef MAKE_CT_INPUT
    int ierr = pip::convert_ck(infile.c_str(), dbfile.c_str(), trfile.c_str(), 
        idtag.c_str());
#else
    int ierr = convert_ck(infile.c_str(), dbfile.c_str(), trfile.c_str(), 
        outfile.c_str(), idtag.c_str());
#endif
    if (ierr < 0) {
        showErrors(cerr);
    }
    return ierr;
}
