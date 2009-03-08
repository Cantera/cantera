/**
 * @file ct2ctml.cpp
 *
 * $Author: dggoodwin $
 * $Revision: 1.38 $
 * $Date: 2006/05/03 20:49:22 $
 */

// Copyright 2001-2005  California Institute of Technology

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
    // environment variable PYTHON_CMD if it is set. If not, return
    // the string 'python'.
    static string pypath() {
        string s = "python";
        const char* py = getenv("PYTHON_CMD");
        if (py) {
            string sp = stripws(string(py));
            if (sp.size() > 0) s = sp;
        }
        return s;
    }


    void ct2ctml(const char* file) {

#ifdef HAS_NO_PYTHON
       /*
        *  Section to bomb out if python is not
        *  present in the computation environment.
        */
        string ppath = file;
        throw CanteraError("ct2ctml", 
                           "python cti to ctml conversion requested for file, " + ppath +
                           ", but not available in this computational environment");
#endif
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
         * This next section may seem a bit weird. However, it is in
         * response to an issue that arises when running cantera with
         * cygwin, using cygwin's python intepreter. Basically, the
         * xml file is written to the local directory by the last
         * system command. Then, the xml file is read immediately
         * after by an ifstream() c++ command. Unfortunately, it seems
         * that the directory info is not being synched fast enough so
         * that the ifstream() read fails, even though the file is
         * actually there. Putting in a sleep system call here fixes
         * this problem. Also, having the xml file pre-existing fixes
         * the problem as well. There may be more direct ways to fix
         * this bug; however, I am not aware of them.
         * HKM -> During the solaris port, I found the same thing.
         *        It probably has to do with NFS syncing problems.
         *        3/3/06
         */
        string sss = sleep();
#ifdef DEBUG_PATHS
        writelog("sleeping for " + sss + " secs+\n");
#endif
#ifndef WIN32
        cmd = "sleep " + sss;
        try {
            ierr = system(cmd.c_str());
        }
        catch (...) {
            ierr = -10;
        }
#endif
        // show the contents of the log file on the screen
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
            throw CanteraError("ct2ctml", 
			       "could not convert input file to CTML.\n "
			       "Command line was: \n" + msg);
        }

        // if the conversion succeeded and DEBUG_PATHS is not defined,
        // then clean up by deleting the temporary Python file.
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
