/**
 * @file ct2ctml.cpp
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "ctexceptions.h"
#include <fstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "ctml.h"

//#define DEBUG_PATHS

using namespace Cantera;

namespace ctml {

    // return the full path to the Python interpreter.  Use
    // environment variable PYTHON_CMD if it is set, otherwise use the
    // definition PYTHON_EXE if it is defined. Failing that, return
    // the string 'python'.
    static string pypath() {
        string s = "python";
        const char* py = getenv("PYTHON_CMD");
        if (py) {
            string sp = stripws(string(py));
            if (sp.size() > 0) s = sp;
        }
        //#ifdef PYTHON_EXE
            //else {
        //string se = stripws(string(PYTHON_EXE));
        //  if (se.size() > 0) s = se;
        //}
        //#endif
        return s;
    }

#ifdef INCL_CHECKPYTHON
    static bool checkPython() {
        time_t aclock;
        time( &aclock );
        int ia = static_cast<int>(aclock);
        string path = tmpDir() + "/.check"+int2str(ia)+".pyw";
        ofstream f(path.c_str());
        if (!f) {
            throw CanteraError("checkPython","cannot open "+path+" for writing");
        }
        f << "import ctml_writer\n";
        f.close();
        int ierr = 0;
#ifdef WIN32
        string cmd = "cmd /C "+pypath() + " " + path + ">> log 2>&1";
#else
        string cmd = pypath() + " " + path + " &> " + tmpDir() + "/log";
#endif
        try {
            ierr = system(cmd.c_str());
            if (ierr != 0) {
                string msg;
                msg = cmd +"\n\n########################################################################\n\n"
                      "The Cantera Python interface is required in order to process\n"
                      "Cantera input files, but it does not seem to be correctly installed.\n\n"
                      "Check that you can invoke the Python interpreter with \n"
                      "the command \"python\", and that typing \"from Cantera import *\"  \n"
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
#ifdef WIN32
        cmd = "cmd /C rm " + path;
#else
        cmd = "rm -f " + path;
        try { 
            system(cmd.c_str());
        }
        catch (...) { ; }
#endif
        return true;
    }
#endif

    void ct2ctml(const char* file) {

        time_t aclock;
        time( &aclock );
        int ia = static_cast<int>(aclock);
        string path = tmpDir()+"/.cttmp"+int2str(ia)+".pyw";
        ofstream f(path.c_str());
        if (!f) {
            throw CanteraError("ct2ctml","cannot open "+path+" for writing.");
        }
        f << "from ctml_writer import *\n"
          << "import sys, os, os.path\n"
          << "file = \"" << file << "\"\n"
          << "base = os.path.basename(file)\n"
          << "root, ext = os.path.splitext(base)\n"
          << "dataset(root)\n"
          << "execfile(file)\n"
          << "write()\n";
        f.close();
#ifdef WIN32
        string cmd = pypath() + " " + path + "> ct2ctml.log 2>&1";
#else
        string cmd = "sleep " + sleep() + "; " + pypath() + 
                     " " + path + " &> ct2ctml.log";
#endif
#ifdef DEBUG_PATHS
        writelog("ct2ctml: executing the command " + cmd + "\n");
#endif
        int ierr = 0;
        try {
            ierr = system(cmd.c_str());
        }
        catch (...) {
            ierr = -10;
        }

        /*
         * HKM -> This next section may seem a bit weird. However,
         *        it is in response to an issue that arises when
         *        running cantera with cygwin, using cygwin's
         *        python intepreter. Basically, the xml file is
         *        written to the local directory by the last
         *        system command. Then, the xml file is read
         *        immediately after by an ifstream() c++
         *        command. Unfortunately, it seems that the
         *        directory info is not being synched fast enough
         *        so that the ifstream() read fails, even
         *        though the file is actually there. Putting in a
         *        sleep system call here fixes this problem. Also,
         *        having the xml file pre-existing fixes the
         *        problem as well. There may be more direct ways
         *        to fix this bug; however, I am not aware of them.
         */
#ifdef CYGWIN 
#ifdef DEBUG_PATHS
        writelog("sleeping for 3 secs+\n");
#endif
        cmd = "sleep 3";
        try {
            ierr = system(cmd.c_str());
        }
        catch (...) {
            ierr = -10;
        }
#endif

        try {
            char ch=0;
            string s = "";
            ifstream ferr("ct2ctml.log");
            if (ferr) {
                while (!ferr.eof()) {
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
            ; 
        }
        if (ierr != 0) {
            string msg = cmd;
            //bool pyok = checkPython();
            //if (!pyok) 
            //    msg += "\nError in Python installation.";
            //else
            //    msg += "\nCheck error messages above for syntax errors.";
            throw CanteraError("ct2ctml", 
			       "could not convert input file to CTML\n "
			       "command line was: " + msg);
        }
#ifndef DEBUG_PATHS
#ifdef WIN32
        cmd = "cmd /C rm " + path;
#else
        cmd = "rm -f " + path;
        try {
            if (ierr == 0) 
                system(cmd.c_str());
        }
        catch (...) { ; }
#endif
#endif
    }


    /**
     * Get an CTML tree from a file, possibly preprocessing the file
     * first. 
     */
    void get_CTML_Tree(XML_Node* rootPtr, string file) {
        string ff, ext = "";

        // find the input file on the Cantera search path
        string inname = findInputFile(file);
#ifdef DEBUG_PATHS
        writelog("Found file: "+inname+"\n");
#endif
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
#ifdef DEBUG_PATHS
        writelog("Attempting to parse xml file " + ff + "\n");
#endif
        ifstream fin(ff.c_str());
        if (!fin) {
          throw 
            CanteraError("get_CTML_tree",
                         "XML file " + ff + " not found");
        }
        rootPtr->build(fin);
        fin.close();
    }
}
