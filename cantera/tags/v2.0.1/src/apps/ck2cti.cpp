/**
 * @file ck2cti.cpp
 *
 * Program to convert Chemkin-II-format reaction mechanism files to
 * Cantera input format. The resulting Cantera input file contains a
 * definition of one ideal_gas entry that represents an ideal gas
 * mixture corresponding to the Chemkin-II reaction mechanism.  The
 * file also contains Cantera-format definitions for each species and
 * each reaction in the input reaction mechanism file.
 *
 * Usage: ck2cti -i input -t thermo -tr transport -id idtag
 *
 * The Cantera-format text is written to the standard output.
 *
 * @param input Chemkin-II reaction mechanism file to be converted. Required.
 *
 * @param thermo Thermodynamic property database. If the THERMO section of the
 * input file is missing or does not have entries for one or more species,
 * this file will be searched for the required thermo data. This file may
 * be another reaction mechanism file containing a THERMO section, or
 * a Chemkin-II-compatible thermodynamic database file.
 *
 * @param transport Transport property database. If this file name is supplied,
 * transport property parameters will be taken from this file and
 * included in the output Cantera-format file. If this parameter is omitted,
 * no transport property parameters will be included in the output.
 *
 * @param id idtag. The ideal_gas entry in the Cantera-format output
 * has name \i idtag. If this parameter is omitted, it will be set to the
 * input file name without the extension. Since only one phase definition
 * is present in the ck2cti output, this parameter is not required.
 */
#include <iostream>
#include <string>
using namespace std;

#include "cantera/base/ct_defs.h"
#include "cantera/base/global.h"
#include "converters/ck2ct.h"

using namespace Cantera;

int showHelp()
{
    cout << "\nck2cti: convert a CK-format reaction mechanism file to Cantera input format.\n"
         << "\n   D. G. Goodwin, Caltech \n"
         << "   Version 1.0, August 2003.\n\n"
         << endl;
    cout << "options:" << endl;
    cout << "    -i  <input file> \n"
         << "    -t  <thermo database> \n"
         << "    -tr <transport database> \n"
         << "    -id  <identifier> \n"
         << "    -d  print debugging output \n\n"
         << "    -v  validate the input file \n\n"
         << "The results are written to the standard output.\n";
    return 0;
}

string getp(int& i, int argc, char** args)
{
    string a="--";
    if (i < argc-1) {
        a = string(args[i+1]);
    }
    if (a[0] == '-') {
        a = "<missing>";
    } else {
        i += 1;
    }
    return a;
}


int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    string infile="chem.inp", dbfile="", trfile="", logfile;
    string idtag = "gas";
    bool debug = false;
    bool validate = false;
    int i=1;
    if (argc == 1) {
        return showHelp();
    }

    while (i < argc) {
        string arg = string(argv[i]);
        if (arg == "-i") {
            infile = getp(i,argc,argv);
        } else if (arg == "-t") {
            dbfile = getp(i,argc,argv);
        } else if (arg == "-tr") {
            trfile = getp(i,argc,argv);
        } else if (arg == "-id") {
            idtag = getp(i,argc,argv);
        } else if (arg == "-d") {
            debug = true;
            cout << "### DEBUG MODE ###" << endl;
        } else if (arg == "-v") {
            validate = true;
            cout << "### VALIDATION ENABLED ###" << endl;
        } else if (arg == "-h" || argc < 3) {
            return showHelp();
        } else {
            cout << "unknown option:" << arg << endl;
            exit(-1);
        }
        ++i;
    }

    int ierr = pip::convert_ck(infile.c_str(), dbfile.c_str(), trfile.c_str(),
                               idtag.c_str(), debug, validate);

    if (ierr < 0) {
        showErrors(std::cerr);
    }
    return ierr;
}
