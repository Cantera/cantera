/**
 * @file ct2ctml.cpp
 * Driver for the system call to the python executable that converts
 * cti files to ctml files (see \ref inputfiles).
 */
/*
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001-2005  California Institute of Technology

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#pragma warning(disable:4996)
#endif

#include "ct_defs.h"
#include "ctexceptions.h"
#include <fstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include "ctml.h"

using namespace std;

#undef DEBUG_PATHS

using namespace Cantera;

namespace ctml {

    //! return the full path to the Python interpreter. 
    /*!
     * Use the environment variable PYTHON_CMD if it is set. If not, return
     * the string 'python'.
     * 
     * Note, there are hidden problems here that really direct us to use
     * a full pathname for the location of python. Basically the system
     * call will use the shell /bin/sh, in order to launch python.
     * This default shell may not be the shell that the user is employing.
     * Therefore, the default path to python may be different during 
     * a system call than during the default user shell environment.
     * This is quite a headache. The answer is to always set the 
     * PYTHON_CMD environmental variable in the user environment to 
     * an absolute path to locate the python executable. Then this 
     * issue goes away. 
     */
    static string pypath() {
        string s = "python";
        const char* py = getenv("PYTHON_CMD");
        if (!py) {
            const char* hm = getenv("HOME");
            string home = stripws(string(hm));
            string cmd = string("source ")+home
                +string("/setup_cantera &> /dev/null");
            system(cmd.c_str());
            py = getenv("PYTHON_CMD");
        }        
        if (py) {
            string sp = stripws(string(py));
            if (sp.size() > 0) {
              s = sp;
            }
        }
        //else {
        //    throw CanteraError("ct2ctml", 
        //        "set environment variable PYTHON_CMD");
        //}
        return s;
    }


    void ct2ctml(const char* file, int debug) {

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
        string path =  tmpDir()+"/.cttmp"+int2str(ia)+".pyw";
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
        string logfile = tmpDir()+"/ct2ctml.log";
#ifdef WIN32
        string cmd = pypath() + " " + "\"" + path + "\"" + "> " + logfile + " 2>&1";
#else
        string cmd = "sleep " + sleep() + "; " + "\"" + pypath() + "\"" + 
                     " " + "\"" + path + "\"" + " &> " + logfile;
#endif
#ifdef DEBUG_PATHS
        writelog("ct2ctml: executing the command " + cmd + "\n");
#endif
        if (debug > 0) {
            writelog("ct2ctml: executing the command " + cmd + "\n");
            writelog("ct2ctml: the Python command is: " + pypath() + "\n");
        }

        int ierr = 0;
        try {
            ierr = system(cmd.c_str());
        }
        catch (...) {
            ierr = -10;
            if (debug > 0)
                writelog("ct2ctml: command execution failed.\n");
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
#ifndef WIN32
        string sss = sleep();
#ifdef DEBUG_PATHS
        writelog("sleeping for " + sss + " secs+\n");
#endif
        if (debug > 0)
            writelog("sleeping for " + sss + " secs+\n");            
        cmd = "sleep " + sss;
        try {
            ierr = system(cmd.c_str());
        }
        catch (...) {
            ierr = -10;
            writelog("ct2ctml: command execution failed.\n");
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
            else {
                if (debug > 0)
                    writelog("cannot open ct2ctml.log for reading.\n");
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
        //#ifdef WIN32
        //cmd = "cmd /C rm " + path;
        if (debug == 0)
            remove(path.c_str());
        else {
            writelog("ct2ctml: retaining temporary file "+path+"\n");
        }
#else
        writelog("ct2ctml: retaining temporary file "+path+"\n");
        //#else
            //cmd = "rm -f \"" + path + "\"";
        //try {
        //    if (ierr == 0) 
        //        system(cmd.c_str());
        //}
        //catch (...) { ; }
        //#endif
#endif
    }


    /**
     * Get an CTML tree from a file, possibly preprocessing the file
     * first. 
     */
    void get_CTML_Tree(XML_Node* rootPtr, string file, int debug) {
        string ff, ext = "";

        // find the input file on the Cantera search path
        string inname = findInputFile(file);
#ifdef DEBUG_PATHS
        writelog("Found file: "+inname+"\n");
#endif
        if (debug > 0)
            writelog("Found file: "+inname+"\n");

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
            ct2ctml(inname.c_str(), debug);
            ff = inname.substr(0,idot) + ".xml";
        }
        else {
            ff = inname;
	}
#ifdef DEBUG_PATHS
        writelog("Attempting to parse xml file " + ff + "\n");
#endif
        if (debug > 0)
            writelog("Attempting to parse xml file " + ff + "\n");
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
