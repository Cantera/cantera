/**
 * @file plots.cpp
 */
// Copyright 2002  California Institute of Technology

#include "cantera/base/plots.h"

using namespace std;

namespace Cantera
{

void writePlotFile(const std::string& fname, const std::string& fmt,
                   const std::string& plotTitle,
                   const std::vector<std::string> &names,
                   const Array2D& data)
{
    ofstream f(fname.c_str());
    if (!f) {
        throw CanteraError("writePlotFile","could not open file "+fname+
                           " for writing.");
    }
    if (fmt == "TEC") {
        outputTEC(f, plotTitle, names, data);
        f.close();
    } else if (fmt == "XL" || fmt == "CSV") {
        outputExcel(f, plotTitle, names, data);
        f.close();
    } else {
        throw CanteraError("writePlotFile",
                           "unsupported plot type:" + fmt);
    }
}

void outputTEC(std::ostream& s, const std::string& title,
               const std::vector<std::string>& names,
               const Array2D& data)
{
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
    for (i = 0; i < nv; i++) {
        s << " SINGLE";
    }
    s << " )" << endl;
    for (i = 0; i < npts; i++) {
        for (j = 0; j < nv; j++) {
            s << data(j,i) << " ";
        }
        s << endl;
    }
}

void outputExcel(std::ostream& s, const std::string& title,
                 const std::vector<std::string>& names,
                 const Array2D& data)
{
    int i,j;
    int npts = static_cast<int>(data.nColumns());
    int nv = static_cast<int>(data.nRows());
    s << title + "," << endl;
    for (i = 0; i < nv; i++) {
        s << names[i];
        if (i != nv-1) {
            s << ",";
        }
    }
    s << endl;
    for (i = 0; i < npts; i++) {
        for (j = 0; j < nv; j++) {
            s << data(j,i);
            if (j != nv-1) {
                s << ",";
            }
        }
        s << endl;
    }
}

}
