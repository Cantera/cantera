/**
 *  @file CKParser.h
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_CKPARSER_H
#define CKR_CKPARSER_H

#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include <fstream>
#include <string>
#include <iostream>
using namespace std;

#include "ckr_defs.h"
#include "Element.h"
#include "Species.h"
#include "Reaction.h"

namespace ckr {


    // typedefs

    /// readability constants
    //@{ 
    const int NoThermoDatabase = 10;
    const int HasTempRange = 11;
    //@}


    /// Exception class for syntax errors.
    class CK_SyntaxError : public CK_Exception {
    public:
        CK_SyntaxError(ostream& f, const string& s, int linenum = -1);
        ostream& m_out;
    };

    /**
     * Chemkin mechanism file parser. For internal use by class CKReader.
     */
    class CKParser {
        
    public:
        
        CKParser(string ckfile, ostream* log);
        CKParser(istream* infile, const string& fname, ostream* log);
    
        /// Destructor.
        ~CKParser(){}

        bool readElementSection(elementList& elements);
        bool readSpeciesSection(speciesList& species);
        bool readThermoSection(vector<string>& names, 
            speciesTable& speciesData, vector_fp& temp, 
            int& optionFlag, ostream& log);
        bool readReactionSection(const vector<string>& speciesNames, 
            vector<string>& elementNames, 
            reactionList& reactions, ReactionUnits& units);
        bool advanceToKeyword(const string& kw, const string& stop);
        bool verbose;
        bool debug;
    
    private:

        int m_line;
        string m_buf;
        string m_comment;
        istream* m_ckfile;
        string m_ckfilename;
        ostream* m_log;
        bool m_nasafmt;
        char m_last_eol;
        void readThermoRecord(Species& sp);
        void getCKLine(string& s, string& comment);    
        void putCKLine(string& s, string& comment);
        void missingAuxData(const string& kw);
    };
}


#endif







