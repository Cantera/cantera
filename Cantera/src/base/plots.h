
#ifndef CT_PLOTS_H
#define CT_PLOTS_H

#include <vector>
#include <string>
#include <fstream>
//using namespace std;

#include "Array.h"
#include "ctexceptions.h"

namespace Cantera {

    void writePlotFile(std::string fname, std::string fmt, 
        std::string plotTitle, const std::vector<std::string> names, const Array2D& data);

    void outputTEC(std::ostream &s, std::string title, const std::vector<std::string>& names,  
        const Array2D& data);

    void outputExcel(std::ostream &s, std::string title, const std::vector<std::string>& names,  
        const Array2D& data);
}

#endif
