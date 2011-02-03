/**
 *  @file filter.cpp
 *
 */

// Copyright 2001  California Institute of Technology


#ifndef CKR_FILTER_H
#define CKR_FILTER_H

#ifdef WIN32
#pragma warning(disable:4786)
#endif

//
//   STL includes
//
#include <map>
#include <string>
#include <vector>
#include <fstream>

using namespace std;
#include <time.h>

//
//   CKReader includes. 
//
#include "CKReader.h"
namespace ckr {

    /**
     *  Edit a mechanism.
     */

    bool filter(const string& infile, const string& database,
		const string& outfile, const vector<int>& species,
		const vector<int>& reactions) {

        bool ok = true;
        
        //
        //  read the input file
        //
        ckr::CKReader ck;
        ck.validate = false;
        string logfile = "filter";
        try {
            if (!ck.read(infile, database, logfile)) {
                cerr << "error encountered while parsing " << infile << endl;
                cerr << "see log file " << logfile << " for details." 
                     << endl << endl;
                return false;
            }
        }
        catch (ckr::CK_SyntaxError) {
            cerr << "syntax error encountered while parsing " 
		 << infile << endl;
            cerr << "see log file " << logfile << " for details." 
                 << endl << endl;
            return false;
        }
        
        ofstream fout(outfile.c_str());
        if (!fout) {
            cerr << "could not open " << outfile << endl;
            return false;
        }
        int nel = static_cast<int>(ck.elements.size());
        int n;

        // write header information
        struct tm *newtime;
        time_t aclock;
        time( &aclock );                  /* Get time in seconds */
        newtime = localtime( &aclock );   /* Convert time to struct tm form */

        fout << "!\n! reduced mechanism generated from " << infile << endl 
             << "! " << asctime(newtime) << endl;

        fout << "ELEMENTS " << endl;
        for (n = 0; n < nel; n++) fout << ck.elements[n];
        fout << "END" << endl;

        int nsp = static_cast<int>(species.size());
        fout << "SPECIES" << endl;
        for (n = 0; n < nsp; n++) {
            fout << ck.species[ species[n] ].name << " ";
            if (5*((n+1)/5) == n+1) fout << endl;
        }
        fout << "END" << endl;
        fout << "REACTIONS" << endl;
        int nrxns = static_cast<int>(reactions.size());
        for (n = 0; n < nrxns; n++) {
            vector<string>& lines = ck.reactions[reactions[n]].lines;
            int nl = static_cast<int>(lines.size());
            for (int j = 0; j < nl; j++) fout << lines[j] << endl;
        }
        fout.close();
        return ok;
    }
}


#endif

    
