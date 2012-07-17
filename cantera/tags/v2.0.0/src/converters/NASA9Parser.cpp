/**
 *  @file CKParser.cpp
 *
 */

// Copyright 2001  California Institute of Technology

#include <numeric>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <math.h>

#include "CKParser.h"
#include "ckr_utils.h"
#include "writelog.h"
#include <stdio.h>
//#include "cantera/base/stringUtils.h"

#include <ctype.h>
#include <string>
#include <cstdlib>

using namespace std;

namespace ckr
{



/**
 *  Add an element to a species.
 *  @param symbol  element symbol
 *  @param atoms   number of atoms of this element in the
 *                 species (may be non-integral)
 *  @param sp      Species object to add element to
 *  @param log     log file output stream
 */
static void addElement(std::string symbol, double atoms,
                       Species& sp, std::ostream& log)
{

    if (atoms != 0.0) {
        Constituent e;
        e.name = symbol;
        e.number = atoms;
        sp.elements.push_back(e);
        sp.comp[symbol] = atoms;
    }
}


/**
 *  Throw an exception if number string is bad
 */
static void illegalNumber(std::ostream& f,
                          std::string s, int linenum = -1)
{
    string msg = "illegal number: "+s;
    throw CK_SyntaxError(f, msg, linenum);
};


void CKParser::checkSpeciesName(std::string spname)
{
    if (spname.size() <= 0) {
        string sss =  "Empty for string name";
        throw CK_SyntaxError(*m_log, sss, m_line);
    }
    char first = spname[0];
    if (isdigit(first)) {
        string sss = "First char of string name is number";
        throw CK_SyntaxError(*m_log, sss, m_line);
    }
    if (isspace(first)) {
        string sss = "First char of  string name is white space";
        throw CK_SyntaxError(*m_log, sss, m_line);
    }
}

static std::string d2e(std::string s)
{
    size_t n;
    size_t sz = s.size();
    string r = s;
    char ch;
    for (n = 0; n < sz; n++) {
        ch = s[n];
        if (ch == 'D') {
            r[n] = 'E';
        } else if (ch == 'd') {
            r[n] = 'e';
        }
    }
    return r;
}

static double de_atof(std::string s)
{
    string r = d2e(s);
    //double rval = Cantera::atofCheck(r.c_str());
    double rval = atof(r.c_str());
    return rval;
}

static double getNumberFromString(std::string s)
{
    bool inexp = false;
    removeWhiteSpace(s);
    int sz = static_cast<int>(s.size());
    char ch;
    for (int n = 0; n < sz; n++) {
        ch = s[n];
        if (!inexp && (ch == 'E' || ch == 'e' || ch == 'D' || ch == 'd')) {
            inexp = true;
        } else if (ch == '+' || ch == '-') {
            if (n > 0 && (s[n-1] != 'E' && s[n-1]
                          != 'e' && s[n-1] != 'd' && s[n-1] != 'D')) {
                return UNDEF;
            }
        } else if (ch != '.' && (ch < '0' || ch > '9')) {
            return UNDEF;
        }
    }
    return de_atof(s);
}

static int de_atoi(std::ostream& log, std::string s, int line = -1)
{
    double val = getNumberFromString(s);
    int ival = (int) val;
    double val2 = (double) ival;
    if (fabs(val - val2) >= 0.00001 * (val + val2)) {
        string sss = "de_atoi: Conversion of int failed: " + s;
        throw CK_SyntaxError(log, sss, line);
    }
    return ival;
}

/**
 *
 * Read species data from THERMO section records.
 *
 * @param names        List of species names (input).
 * @param species      Table of species objects holding data from records
 *                     in THERMO section (output).
 * @param temp         Default vector of temperature region boundaries
 *                     There are one more temperatures than there are
 *                     temperature regions.
 *
 * @return            True, if the THERMO section exists and the species
 *                    have all been successfully processed. False, if
 *                    the THERMO section doesn't exist or there were
 *                    additional problems.
 */
bool CKParser::readNASA9ThermoSection(std::vector<string>& names,
                                      speciesTable& species, vector_fp& temp,
                                      int& optionFlag, std::ostream& log)
{
    // String buffer for lines
    string s;
    vector<string> toks;
    string defaultDate="";
    int nsp = static_cast<int>(names.size());

    // Comment string
    string comment;

    // Check to see that we expect to be reading a NASA9 formatted file
    if (!m_nasa9fmt) {
        throw CK_SyntaxError(log,
                             "In NASA9 parser. However, we expect a different file format",
                             -1);
    }

    // now read in all species records that have names in list 'names'

    bool getAllSpecies = (nsp > 0 && match(names[0], "<ALL>"));
    if (getAllSpecies) {
        names.clear();
    }

    // Map between the number of times a species name appears in the database
    map<string, int> dup; // used to check for duplicate THERMO records
    bool already_read;

    while (1 > 0) {
        // If we don't have any more species to read, break
        if (nsp == 0) {
            break;
        }
        already_read = false;

        // Read a new species record from the section
        Species spec;
        readNASA9ThermoRecord(spec);

        // we signal the end of the section by putting <END> as a
        // species name. Break if you find the end of the section.
        if (spec.name == "<END>") {
            break;
        }

        // check for duplicate thermo data
        if (dup[spec.name] == 2) {
            log << "Warning: more than one THERMO record for "
                << "species " << spec.name << endl;
            log << "Record at line " << m_line
                << " of " << m_ckfilename << " ignored." << endl;
            already_read = true;
        }
        // Set the record in the map to 2 to create a signal for the
        // next time.
        dup[spec.name] = 2;

        // Check to see whether we need this particular species name
        if (!already_read && (getAllSpecies
                              || (find(names.begin(), names.end(), spec.name)
                                  < names.end()))) {

            // Add the species object to the map. Note we are
            // doing a copy constructor here, so we create a
            // lasting entry.
            species[spec.name] = spec;

            if (verbose) {
                log << endl << "found species " << spec.name;
                log << " at line " << m_line
                    << " of " << m_ckfilename;
                writeSpeciesData(log, spec);
            }
            if (getAllSpecies) {
                names.push_back(spec.name);
                nsp = static_cast<int>(names.size());
            } else {
                nsp--;
            }
        }
    }
    return true;
}


/**
 *
 * Read one species definition in a NASA9 string.
 *
 */
void CKParser::readNASA9ThermoRecord(Species& sp)
{
    string s;
    string numstr;
    double cf;
    // Set to the NASA9 polynomial format
    sp.thermoFormatType = 1;

    // look for line 1, but if a keyword is found first or the end of
    // the file is reached, return "<END>" as the species name
    string comment;

    // Name of the species
    string nameid;
    vector<string> toks;
    size_t nToks = 0;

    // Loop forward until we get to the next nonempty line.
    do {
        getCKLine(s, comment);
        if (isKeyword(s) || match(s, "<EOF>")) {
            sp.name = "<END>";
            putCKLine(s, comment);
            return;
        }

        // The first 18 spaces are devoted to the name of the species
        string nameid = s.substr(0,18);
        getTokens(nameid, static_cast<int>(nameid.size()), toks);
        nToks = toks.size();
    }  while (nToks == 0);

    //------------- line 1 ---------------------------
    //  Everything after the first 18 spaces is a comment.
    size_t nt = s.size();
    sp.m_commentsRef = s.substr(18, nt-18);

    // Parse the species name
    sp.name = toks[0];
    sp.id = "";
    if (nToks > 1) {
        throw CK_SyntaxError(*m_log,
                             "Illegal number of tokens for string name", m_line);
    }
    checkSpeciesName(sp.name);


    //------------- line 2 ---------------------------

    getCKLine(s, comment);
    if (s.size() < 79) {
        throw CK_SyntaxError(*m_log,
                             "Size of second line is too small", m_line);
    }
    // Read the number of temperature regions.
    string sN = s.substr(0,2);
    sp.nTempRegions = de_atoi(*m_log, sN);

    string refDataCode = s.substr(3,6);

    // elemental composition (first 5 elements)
    for (int i = 0; i < 5; i++) {
        string elementSym = "";
        int iloc = 10 + 8*i;
        if (s[iloc] != ' ') {
            if (s[iloc+1] != ' ') {
                elementSym = s.substr(iloc,2);
            } else {
                elementSym = s.substr(iloc,1);
            }
        } else if (s[iloc+1] != ' ') {
            elementSym = s.substr(iloc+1,1);
        }
        double atoms = de_atof(s.substr(iloc+2,6));
        addElement(elementSym, atoms, sp, *m_log);
    }


    // single-character phase descriptor
    sp.phase = s.substr(50,2);

    // Molecular weight in gm per gmol
    string molecWeight = s.substr(52, 13);

    // Heat of formation at 298.15 K in J / gmol
    string Hf298_Jgmol = s.substr(65, 15);

    vector_fp* coeffs_ptr;
    for (int i = 0; i < sp.nTempRegions; i++) {

        coeffs_ptr = new vector_fp(9);
        vector_fp& coeffs = *coeffs_ptr;

        //------------- line 3 ---------------------------
        getCKLine(s, comment);
        if (s.size() < 79) {
            throw CK_SyntaxError(*m_log,
                                 "Size of third line is too small", m_line);
        }

        string sTlow = s.substr(0, 11);
        double tLow = de_atof(sTlow);

        string sTHigh = s.substr(11, 11);
        double tHigh = de_atof(sTHigh);

        string sNCoeff = s.substr(22, 1);
        int nCoeff = de_atoi(*m_log, sNCoeff);
        if (nCoeff != 7) {
            throw CK_SyntaxError(*m_log, "ncoeff ne 7", m_line);
        }

        string sTCoeff1 = s.substr(24, 5);
        double TCoeff1 = de_atof(sTCoeff1);
        if (TCoeff1 != -2.0) {
            throw CK_SyntaxError(*m_log, "TCoeff1 ne -2.0", m_line);
        }

        string sTCoeff2 = s.substr(29, 5);
        double TCoeff2 = de_atof(sTCoeff2);
        if (TCoeff2 != -1.0) {
            throw CK_SyntaxError(*m_log, "TCoeff2 ne -1.0", m_line);
        }

        string sTCoeff3 = s.substr(34, 5);
        double TCoeff3 = de_atof(sTCoeff3);
        if (TCoeff3 != 0.0) {
            throw CK_SyntaxError(*m_log, "TCoeff3 ne 0.0", m_line);
        }

        string sTCoeff4 = s.substr(39, 5);
        double TCoeff4 = de_atof(sTCoeff4);
        if (TCoeff4 != 1.0) {
            throw CK_SyntaxError(*m_log, "TCoeff4 ne 1.0", m_line);
        }

        string sTCoeff5 = s.substr(44, 5);
        double TCoeff5 = de_atof(sTCoeff5);
        if (TCoeff5 != 2.0) {
            throw CK_SyntaxError(*m_log, "TCoeff5 ne 2.0", m_line);
        }

        string sTCoeff6 = s.substr(49, 5);
        double TCoeff6 = de_atof(sTCoeff6);
        if (TCoeff6 != 3.0) {
            throw CK_SyntaxError(*m_log, "TCoeff6 ne 3.0", m_line);
        }

        string sTCoeff7 = s.substr(54, 5);
        double TCoeff7 = de_atof(sTCoeff7);
        if (TCoeff7 != 4.0) {
            throw CK_SyntaxError(*m_log, "TCoeff7 ne 4.0", m_line);
        }

        string sHf298mHF0 = s.substr(65, 15);

        //------------- line 4 ---------------------------
        getCKLine(s, comment);
        if (s.size() < 79) {
            throw CK_SyntaxError(*m_log,
                                 "Size of third line is too small", m_line);
        }

        numstr = s.substr(0, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[0] = cf;

        numstr = s.substr(16, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[1] = cf;

        numstr = s.substr(32, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[2] = cf;

        numstr = s.substr(48, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[3] = cf;

        numstr = s.substr(64, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[4] = cf;

        //------------- line 5 ---------------------------
        getCKLine(s, comment);
        if (s.size() < 79) {
            throw CK_SyntaxError(*m_log,
                                 "Size of fourth line is too small", m_line);
        }

        numstr = s.substr(0, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[5] = cf;

        numstr = s.substr(16, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[6] = cf;

        numstr = s.substr(48, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[7] = cf;

        numstr = s.substr(64, 16);
        cf = getNumberFromString(numstr);
        if (cf == UNDEF) {
            illegalNumber(*m_log, numstr, m_line);
        }
        coeffs[8] = cf;

        // Store the coefficients.
        sp.minTemps.push_back(tLow);
        sp.maxTemps.push_back(tHigh);

        sp.region_coeffs.push_back(coeffs_ptr);

    }

    sp.valid = 1;
}



}  // ckr namespace

