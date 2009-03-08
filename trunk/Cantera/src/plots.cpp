/**
 * @file plots.cpp
 */

/*
 * $Author: hkmoffa $
 * $Revision: 1.2 $
 * $Date: 2004/07/01 23:46:52 $
 */

// Copyright 2002  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "plots.h"

namespace Cantera {

    void writePlotFile(string fname, string fmt, 
        string plotTitle, const vector<string> names, const Array2D& data) {
        ofstream f(fname.c_str());
        if (!f) {
            throw CanteraError("writePlotFile","could not open file "+fname+
                " for writing.");
        }
        if (fmt == "TEC") {
            outputTEC(f, plotTitle, names, data);
            f.close();
        }
        else if (fmt == "XL") {
            outputExcel(f, plotTitle, names, data);
            f.close();
        }
        else {
            throw CanteraError("writePlotFile",
                "unsupported plot type:" + fmt);
        }
    }
            

    /**
     * Write a Tecplot data file.  
     * @param s output stream 
     * @param title plot title 
     * @param names vector of variable names 
     * @param data N x M data array. 
     * data(n,m) is the m^th value of the n^th variable.
     */
    void outputTEC(ostream &s, string title, const vector<string>& names,  
        const Array2D& data) {
        int i,j;
        int npts = static_cast<int>(data.nColumns());
        int nv = static_cast<int>(data.nRows());
        s << "TITLE     = \"" + title + "\"" << endl;
        s << "VARIABLES = " << endl;
        for (i = 0; i < nv; i++) {
            s << "\"" << names[i] << "\"" << endl;
        }
        s << "ZONE T=\"zone1\"" << endl;
        s << " I=" << npts << ",J=1,K=1,F=POINT" << endl;
        s << "DT=( ";
        for (i = 0; i < nv; i++) s << " SINGLE";
        s << " )" << endl;
        for (i = 0; i < npts; i++) {
            for (j = 0; j < nv; j++) {
                s << data(j,i) << " ";
            }
            s << endl;
        }
    }


    /**
     * Write an Excel spreadsheet in 'csv' form.  
     * @param s output stream 
     * @param title plot title 
     * @param names vector of variable names @param data N
     * x M data array. data(n,m) is the m^th value of the n^th variable.
     */
    void outputExcel(ostream &s, string title, const vector<string>& names,  
        const Array2D& data) {
        int i,j;
        int npts = static_cast<int>(data.nColumns());
        int nv = static_cast<int>(data.nRows());
        s << title + "," << endl;
        for (i = 0; i < nv; i++) {
            s << names[i] << ",";
        }
        s << endl;
        for (i = 0; i < npts; i++) {
            for (j = 0; j < nv; j++) {
                s << data(j,i) << ",";
            }
            s << endl;
        }
    }


}
