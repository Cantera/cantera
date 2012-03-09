/**
 * @file ct2ctml.cpp
 * Driver for the system call to the python executable that converts
 * cti files to ctml files (see \ref inputfiles).
 */
// Copyright 2001-2005  California Institute of Technology

#include "cantera/base/ct_defs.h"
#include "cantera/base/ctexceptions.h"
#include "cantera/base/ctml.h"
#include "cantera/base/global.h"
#include "cantera/base/stringUtils.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>

// These defines are needed for the windows Sleep() function
// - comment them out if you don't want the Sleep function.
//#ifdef _WIN32
//#include "Windows.h"
//#include "Winbase.h"
//#endif

using namespace Cantera;
using namespace std;

namespace ctml
{

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
static string pypath()
{
    string s = "python";
    const char* py = getenv("PYTHON_CMD");

    if (py) {
        string sp = stripws(string(py));
        if (sp.size() > 0) {
            s = sp;
        }
    }
    return s;
}

// Convert a cti file into a ctml file
/*
 *
 *  @param   file    Pointer to the file
 *  @param   debug   Turn on debug printing
 *
 *  @ingroup inputfiles
 */
void ct2ctml(const char* file, const int debug)
{

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
    time(&aclock);
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
#ifdef _WIN32
    string cmd = pypath() + " " + "\"" + path + "\"" + "> " + logfile + " 2>&1";
#else
    string cmd = "sleep " + sleep() + "; " + "\"" + pypath() + "\"" +
                 " " + "\"" + path + "\"" + " &> " + logfile;
#endif
    if (debug > 0) {
        writelog("ct2ctml: executing the command " + cmd + "\n");
        writelog("ct2ctml: the Python command is: " + pypath() + "\n");
    }

    int ierr = 0;
    try {
        ierr = system(cmd.c_str());
    } catch (...) {
        ierr = -10;
        if (debug > 0) {
            writelog("ct2ctml: command execution failed.\n");
        }
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
#ifndef _WIN32
    string sss = sleep();
    if (debug > 0) {
        writelog("sleeping for " + sss + " secs+\n");
    }
    cmd = "sleep " + sss;
    try {
        ierr = system(cmd.c_str());
    } catch (...) {
        ierr = -10;
        writelog("ct2ctml: command execution failed.\n");
    }
#else
    // This command works on windows machines if Windows.h and Winbase.h are included
    // Sleep(5000);
#endif

    if (ierr != 0) {
        // Generate an error message that includes the contents of the
        // ct2ctml log file.
        stringstream message;
        message << "Error converting input file \"" << file << "\" to CTML.\n";
        message << "Command was:\n\n";
        message << cmd << "\n\n";
        ifstream ferr(logfile.c_str());
        if (ferr) {
            message << "-------------- start of ct2ctml.log --------------\n";
            message << ferr.rdbuf();
            message << "--------------- end of ct2ctml.log ---------------";
        } else {
            message << "Additionally, the contents of ct2ctml.log"
                "could not be read";
        }
        throw CanteraError("ct2ctml", message.str());
    }

    // If the conversion succeeded and no debugging information is needed,
    // clean up by deleting the temporary Python file and the log file.
    if (debug == 0) {
        remove(path.c_str());
        remove(logfile.c_str());
    } else {
        writelog("ct2ctml: retaining temporary file "+path+"\n");
        writelog("ct2ctml: retaining temporary file "+logfile+"\n");
    }
}


// Read an ctml file from a file and fill up an XML tree
/*
 *  This is the main routine that reads a ctml file and puts it into
 *  an XML_Node tree
 *
 *  @param node    Root of the tree
 *  @param file    Name of the file
 *  @param debug   Turn on debugging printing
 */
void get_CTML_Tree(Cantera::XML_Node* rootPtr, const std::string file, const int debug)
{
    std::string ff, ext = "";

    // find the input file on the Cantera search path
    std::string inname = findInputFile(file);
    if (debug > 0) {
        writelog("Found file: "+inname+"\n");
    }

    if (inname == "") {
        throw CanteraError("get_CTML_Tree", "file "+file+" not found");
    }

    /*
     * Check whether or not the file is XML. If not, it will be first
     * processed with the preprocessor.
     */
    std::string::size_type idot = inname.rfind('.');
    if (idot != string::npos) {
        ext = inname.substr(idot, inname.size());
    }
    if (ext != ".xml" && ext != ".ctml") {
        try {
            ctml::ct2ctml(inname.c_str(), debug);
        } catch (std::exception& err) {
            writelog("get_CTML_Tree: caught an exception:\n");
            writelog(err.what());
        }
        string ffull = inname.substr(0,idot) + ".xml";
        ff = "./" + getBaseName(ffull) + ".xml";
        if (debug > 0) {
            writelogf("ffull name = %s\n", ffull.c_str());
            writelogf("ff name = %s\n", ff.c_str());
        }
    } else {
        ff = inname;
    }
    if (debug > 0) {
        writelog("Attempting to parse xml file " + ff + "\n");
    }
    ifstream fin(ff.c_str());
    if (!fin) {
        throw
        CanteraError("get_CTML_Tree",
                     "XML file " + ff + " not found");
    }
    rootPtr->build(fin);
    fin.close();
}
}
