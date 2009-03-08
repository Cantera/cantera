
#ifndef CT_PLOTS_H
#define CT_PLOTS_H

#include <vector>
#include <string>
#include <fstream>
using namespace std;

#include "Array.h"
#include "ctexceptions.h"

namespace Cantera {

    void writePlotFile(string fname, string fmt, 
        string plotTitle, const vector<string> names, const Array2D& data);

    void outputTEC(ostream &s, string title, const vector<string>& names,  
        const Array2D& data);

    void outputExcel(ostream &s, string title, const vector<string>& names,  
        const Array2D& data);
}

#endif
