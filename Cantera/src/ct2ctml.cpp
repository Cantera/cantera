/**
 * @file ct2ctml.cpp
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#include "ct_defs.h"
#include "ctexceptions.h"
#include <fstream>
#include <string>
#include <stdlib.h>
#include "ctml.h"

using namespace Cantera;

namespace ctml {

    static string pypath() {
        string s = "python";
#ifdef PYTHON_EXE
        s = string(PYTHON_EXE);
        return s;
#else
        const char* py = getenv("PYTHON_CMD");
        if (py) s = string(py);
        return s;
#endif
    }

    static bool checkPython() {
        string path = tmpDir() + "/check.py";
        ofstream f(path.c_str());
        if (!f) {
            throw CanteraError("checkPython","cannot write to "+tmpDir());
        }
        f << "from Cantera import *\n";
        f.close();
        int ierr = 0;
        string cmd = pypath() + " " + path + " &> " + tmpDir() + "/log";
        try {
            ierr = system(cmd.c_str());
            if (ierr != 0) {
                string msg;
                msg = cmd + "\n\n########################################################################\n\n"
                      "The Cantera Python interface is required in order to process\n"
                      "Cantera input files, but it does not seem to be correctly installed.\n\n"
                      "Check that you can invoke the Python interpreter with \n"
                      "the command \"python\", and that typing \"from Cantera import *\"\n"
                      "at the Python prompt does not produce an error. If Python on your system\n"
                      "is invoked with some other command, set environment variable PYTHON_CMD\n"
                      "to the full path to the Python interpreter. \n\n"
                      "#########################################################################\n\n";
                writelog(msg);
                return false;
            }
        }
        catch (...) {
            return false;
        }
        return true;
    }


    void ct2ctml(const char* file) {
        string path = tmpDir()+"/.cttmp.py";
        ofstream f(path.c_str());
        if (!f) {
            throw CanteraError("ct2ctml","cannot write to "+tmpDir());
        }
        f << "from Cantera import *\n"
          << "from Cantera.ctml_writer import *\n"
          << "import sys, os, os.path\n"
          << "file = \"" << file << "\"\n"
          << "base = os.path.basename(file)\n"
          << "root, ext = os.path.splitext(base)\n"
          << "dataset(root)\n"
          << "execfile(file)\n"
          << "write()\n";
        f.close();
        string cmd = pypath() + " " + path + " &> ct2ctml.log";
        int ierr;
        try {
            ierr = system(cmd.c_str());
        }
        catch (...) {
            ierr = -10;
        }
        //char line[90];
        //if (ierr != 0) {
        try {
            char ch;
            string s = "";
            ifstream ferr("ct2ctml.log");
            if (ferr) {
                while (!ferr.eof()) {
                    //msg += "\n";
                    ferr.get(ch);
                    s += ch;
                    if (ch == '\n') {
                        writelog(s);
                        s = "";
                    }
                }
                ferr.close();
            }
        }
        catch (...) {
            ; //writelog("could not print error message.\n");
        }
        if (ierr != 0) {
            string msg = cmd;
            bool pyok = checkPython();
            if (!pyok) 
                msg += "\nError in Python installation.";
            else
                msg += "\nCheck error messages above for syntax errors.";
            throw CanteraError("ct2ctml", 
			       "could not convert input file to CTML\n "
			       "command line was: " + msg);
        }
    }


    /**
     * Get an CTML tree from a file, possibly preprocessing the file
     * first. 
     */
    void get_CTML_Tree(XML_Node* rootPtr, string file) {
        string ff, ext = "";

        // find the input file on the Cantera search path
        string inname = findInputFile(file);
        if (inname == "") 
            throw CanteraError("get_CTML_tree", "file "+file+" not found");

        /** 
         * Check whether or not the file is XML. If not, it will be first
         * processed with the preprocessor.
         */
        string::size_type idot = inname.rfind('.');
	if (idot != string::npos) {
	  ext = inname.substr(idot, inname.size());
	}
        if (ext != ".xml" && ext != ".ctml") {
	  ct2ctml(inname.c_str());
	  ff = inname.substr(0,idot) + ".xml";
        }
        else {
	  ff = inname;
	}
        //        rootPtr = &get_XML_File(ff)->child("ctml");
        // rootPtr->write(cout);
        ifstream fin(ff.c_str());
        rootPtr->build(fin);
        fin.close();
    }
}
