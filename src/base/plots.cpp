/**
 * @file plots.cpp
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/plots.h"

using namespace std;

namespace Cantera
{

void writePlotFile(const std::string& fname, const std::string& fmt,
                   const std::string& plotTitle,
                   const std::vector<std::string> &names,
                   const Array2D& data)
{
    ofstream f(fname);
    if (!f) {
        throw CanteraError("writePlotFile","could not open file "+fname+
                           " for writing.");
    }
    if (fmt == "TEC") {
        outputTEC(f, plotTitle, names, data);
    } else if (fmt == "XL" || fmt == "CSV") {
        outputExcel(f, plotTitle, names, data);
    } else {
        throw CanteraError("writePlotFile",
                           "unsupported plot type:" + fmt);
    }
}

void outputTEC(std::ostream& s, const std::string& title,
               const std::vector<std::string>& names,
               const Array2D& data)
{
    s << "TITLE     = \"" + title + "\"" << endl;
    s << "VARIABLES = " << endl;
    for (size_t i = 0; i < data.nRows(); i++) {
        s << "\"" << names[i] << "\"" << endl;
    }
    s << "ZONE T=\"zone1\"" << endl;
    s << " I=" << data.nColumns() << ",J=1,K=1,F=POINT" << endl;
    s << "DT=( ";
    for (size_t i = 0; i < data.nRows(); i++) {
        s << " SINGLE";
    }
    s << " )" << endl;
    for (size_t i = 0; i < data.nColumns(); i++) {
        for (size_t j = 0; j < data.nRows(); j++) {
            s << data(j,i) << " ";
        }
        s << endl;
    }
}

void outputExcel(std::ostream& s, const std::string& title,
                 const std::vector<std::string>& names,
                 const Array2D& data)
{
    s << title + "," << endl;
    for (size_t i = 0; i < data.nRows(); i++) {
        s << names[i];
        if (i != data.nRows()-1) {
            s << ",";
        }
    }
    s << endl;
    for (size_t i = 0; i < data.nColumns(); i++) {
        for (size_t j = 0; j < data.nRows(); j++) {
            s << data(j,i);
            if (j != data.nRows()-1) {
                s << ",";
            }
        }
        s << endl;
    }
}

}
