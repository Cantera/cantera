#include "ct_defs.h"
#include <ofstream>
#include <string>
#include <stdlib.h>

namespace Cantera {
void ct2ctml(char* file) {
    ofstream f(".cttmp.py");
    f << "from Cantera import *\n"
      << "from Cantera.ctml_writer import *\n"
      << "import sys, os, os.path\n"
      << "file = " << file << endl
      << "base = os.path.basename(file)\n"
      << "root, ext = os.path.splitext(base)\n"
      << "dataset(root)\n"
      << "execfile(file)\n"
      << "write()\n";
    f.close();
    string PY_CMD = getenv("PYTHON_CMD");
    if (!PY_CMD) PY_CMD = "python";
    int ierr = system(PY_CMD+" .cttmp.py");
    if (ierr != 0) 
        throw CanteraError("ct2ctml", "could not convert input file to CTML.");
}

}
