/**
 *  @file misc.cpp
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "global.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "units.h"
#include "xml.h"
#include "ctml.h"

#ifndef WIN32
#include "ctdir.h"
#endif

#include <fstream>
using namespace std;

namespace Cantera {

    /**
     * Class to hold global data. Class Application is the top-level
     * class that stores data that should persist for the duration of
     * the process. The class should not be instantiated directly;
     * instead, it is instantiated as needed by the functions declared
     * here. At most one instance is created, and it is not destroyed
     * until the process terminates.
     */
    class Application {
    public:
        Application() : linelen(0), stop_on_error(false), 
                        tmp_dir("/tmp") {}
        virtual ~Application(){}
        vector<string> inputDirs;
        vector<string> errorMessage;
        vector<string> warning;
        vector<string> errorRoutine;
        string msglog;
        size_t linelen;
        bool stop_on_error;
        //bool write_log_to_cout;
        //bool matlab;
        map<string, string>     options;
        string tmp_dir;
        map<string, XML_Node*> xmlfiles;
    };


    /// Returns a pointer to the one and only instance of Application
    Application* app();

    void setDefaultDirectories();

    static Application* __app = 0;
    Unit* Unit::__u = 0;

    static void appinit() {
        if (__app == 0) __app = new Application;
    }

    Application* app() {
        if (__app == 0) {
            __app = new Application;
            setDefaultDirectories();
        }
        return __app;
    }

    XML_Node* get_XML_File(string file) {
        if (app()->xmlfiles.find(file) 
            == app()->xmlfiles.end()) {
            string path = findInputFile(file);

            /*
             * Check whether or not the file is XML. If not, it will
             * be first processed with the preprocessor.
             */
            string::size_type idot = path.rfind('.');
            string ext, ff;
            if (idot != string::npos) {
                ext = path.substr(idot, path.size());
            }
            if (ext != ".xml" && ext != ".ctml") {
                ctml::ct2ctml(path.c_str());
                string::size_type islash = path.rfind('/');
                if (islash != string::npos) 
                    ff = string("./")+path.substr(islash+1,idot-islash - 1) + ".xml";
                else
                    ff = string("./")+path.substr(0,idot) + ".xml";
                writelog("ff = "+ff+"\n");
            }
            else {
                ff = path;
            }
            ifstream s(ff.c_str());
            XML_Node* x = new XML_Node("doc");
            if (s) {
                x->build(s);
                __app->xmlfiles[file] = x;
                __app->xmlfiles[ff] = x;
            }
            else
                throw CanteraError("get_XML_File","cannot open "+ff+" for reading.");

	    /*
	     * Add the built XML Tree to the map, xmlfiles.
	     * It stores the pointer to the tree, with the
	     * key being the name of the file.
	     *
	     * HKM Note: shouldn't the key be the full pathname of
	     *           the file, i.e., ff?
	     */
        }
        return __app->xmlfiles[file];
    }

    void close_XML_File(string file) {
        if (file == "all") {
            map<string, XML_Node*>::iterator 
                b = app()->xmlfiles.begin(), e = app()->xmlfiles.end();
            for(; b != e; ++b) {
                delete b->second;
                __app->xmlfiles.erase(b->first);
            }
        }
        else if (app()->xmlfiles.find(file) 
            != app()->xmlfiles.end()) {
            delete __app->xmlfiles[file];
            __app->xmlfiles.erase(file);
        }
    }

    void setTmpDir(string tmp) { app()->tmp_dir = tmp; }
    string tmpDir() { appinit(); return app()->tmp_dir; }

    int nErrors() {return app()->errorMessage.size();}

    void popError() {
        appinit();
        if (nErrors() > 0) {
            __app->errorMessage.pop_back();
            __app->errorRoutine.pop_back();
        }
    }

    string lastErrorMessage() {
        appinit();
        if (nErrors() > 0)
            return "\nProcedure: "+__app->errorRoutine.back()+"\nError:   "+__app->errorMessage.back();
        else
            return "<no Cantera error>";
    }

    void showErrors(ostream& f) {
        appinit(); 
        int i = __app->errorMessage.size();
        if (i == 0) return;
        f << endl << endl;
        f << "************************************************" << endl;
        f << "                Cantera Error!                  " << endl;
        f << "************************************************" << endl << endl;
        int j;
        for (j = 0; j < i; j++) {
            f << endl;
            f << "Procedure: " << __app->errorRoutine[j] << endl;
            f << "Error:     " << __app->errorMessage[j] << endl;
        } 
        f << endl << endl;
        __app->errorMessage.clear();
        __app->errorRoutine.clear();
        //if (__app->stop_on_error) exit(-1);
    }

    void setError(string r, string msg) {
        appinit();
        __app->errorMessage.push_back(msg);
        __app->errorRoutine.push_back(r);
    }


    /**
     * Set the default directories for input files. Four directories are
     * added to the search path used by findInputFile. These are
     * 'data', 'data/inputs', 'data/thermo', and
     * 'data/transport'. These names are for convenience only -
     * findInputFile searches all of them, independent of the type of
     * file. The location of the 'data' directory depends on how
     * environment variables are set. If CANTERA_DATA_DIR is set, then
     * this will be used instead of 'data'. In addition, if
     * WIN_CANTERA_ROOT or CANTERA_ROOT are set, then 'data' is
     * assumed to be a top-level subdirectory. WIN_CANTERA_ROOT should
     * only be set on PCs, and should be in 'DOS' format, for example
     * 'C:\CANTERA'. CANTERA_ROOT, on the other hand, should be in
     * unix-like format ('/home/usr/cantera'). This allows Cantera to
     * be built on PCs using a unix-like environment (Cygwin) and
     * compiler (g++), as well as using Win32 compilers.
     */
    void setDefaultDirectories() {
        appinit();
        vector<string>& dirs = __app->inputDirs;

        // always look in the local directory first
        dirs.push_back(".");

#ifdef WIN32
        /*
         * Under Windows, the Cantera setup utility puts data files in
         * a directory 'Cantera\data' below the one the environment
         * variable COMMONPROGRAMFILES points to. (This is usually
         * C:\Program Files\Common Files.) If this environment
         * variable is defined, then this directory is assumed to
         * exist and is added to the search path.
         */
        const char* comfiles = getenv("COMMONPROGRAMFILES");
        if (comfiles != 0) {
            string cfiles = string(comfiles);

            // remove quotes if necessary
            if (cfiles[0] == '\'') 
                cfiles = cfiles.substr(1,1000);
            if (cfiles[cfiles.size()-1] == '\'') cfiles[cfiles.size()-1] = '\n';

            string datadir = string(comfiles) + "/Cantera/data";
            string tmpldir = string(comfiles) + "/Cantera/templates";
            dirs.push_back(datadir);
            dirs.push_back(tmpldir);
        }
#endif
        if (getenv("CANTERA_DATA") != 0) {
            string datadir = string(getenv("CANTERA_DATA"));
            dirs.push_back(datadir);
        }

        // CANTERA_ROOT is defined in file ctdir.h. This file is written
        // during the build process (unix), and points to the directory
        // specified by the 'prefix' option to 'configure', or else to
        // /usr/local/cantera. This will generally not be defined under
        // Windows.
#ifdef CANTERA_ROOT
	string datadir = string(CANTERA_ROOT) + "/data";
	dirs.push_back(datadir);
#endif
    }




    void addDirectory(string dir) {
        appinit();
        if (__app->inputDirs.size() == 0) setDefaultDirectories();
        __app->inputDirs.push_back(stripnonprint(dir));
    }

    /**
     *  findInputFile():
     *
     *    This routine will search for a file in the default
     *    locations specified for the application.
     *    See the routine setDefaultDirectories() listed above.
     *
     *    The default set of directories specified for the application
     *    will be searched if a '/' or an '\\' is not found in
     *    name. If either is found then a relative path name is
     *    presumed and the default directories are not searched.
     *
     *    The presence of the file is determined by whether the file
     *    can be opened for reading by the current user.
     *
     *    Return
     *    -------
     *      The absolute path name of the first matching
     *      file is returned. If a relative path name
     *      is indicated, the relative path name is returned.
     *  
     *      If the file is not found, a message is written to 
     *      stdout and  a CanteraError exception is thrown.
     */
    string findInputFile(string name) {
        appinit();
        int islash = name.find('/');
        int ibslash = name.find('\\');
        string inname;
        vector<string>& dirs = __app->inputDirs;
        if (dirs.size() == 0) setDefaultDirectories();

        int nd;
        if (islash < 0 && ibslash < 0) {
            nd = dirs.size();
            int i;
            inname = "";
            for (i = 0; i < nd; i++) {
                inname = dirs[i] + "/" + name;
                writelog("looking for file "+inname+"\n");
                ifstream fin(inname.c_str());
                if (fin) {
                    writelog("found it\n");
                    fin.close();
                    return inname;
                }
                else
                    writelog("not found.\n");
            }
            string msg;
            msg = "\nInput file " + name 
                 + " not found in director";
            msg += (nd == 1 ? "y " : "ies ");
            for (i = 0; i < nd; i++) {
                msg += "\n'" + dirs[i] + "'";
                if (i < nd-1) msg += ", ";
            }
            msg += "\n\n";
            msg += "To fix this problem, either:\n";
            msg += "    a) move the missing files into the local directory;\n";
            msg += "    b) define environment variable CANTERA_DATA to\n";
            msg += "         point to the directory containing the file.";
            throw CanteraError("findInputFile", msg);
            return "";
        }
        else {
            return name;
        }
    }

    //void setMatlabMode(bool m) { 
    //    appinit();
    //    __app->matlab = m; 
    //}

    //void write(const string& msg) {cout << msg;}
    //void write(const char* msg) {cout << msg;}

    void writelog(const char* msg) {writelog(string(msg));}
    //void getlog(string& s) {
    //    appinit();
    //    s = __app->msglog;
    //    //__app->msglog = "";
    //}
    //void clearlog() {
    //    __app->msglog = "";
    //}

    doublereal toSI(string unit) {
        doublereal f = Unit::units()->toSI(unit);
        if (f) return f;
        else return 1.0;
    }

    doublereal actEnergyToSI(string unit) {
        doublereal f = Unit::units()->actEnergyToSI(unit);
        if (f) return f;
        else return 1.0;
    }


    string canteraRoot() {
        char* ctroot = 0;
        ctroot = getenv("CANTERA_ROOT");
        if (ctroot != 0) { return string(ctroot); }
        else {
#ifdef CANTERA_ROOT
            return string(CANTERA_ROOT);
#else
            return "";
#endif
        }
    }
}


