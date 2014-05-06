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
#include "../../ext/libexecstream/exec-stream.h"

#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <windows.h>
#endif

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

void ct2ctml(const char* file, const int debug)
{
    string xml = ct2ctml_string(file);
    string out_name = file;
#ifdef _WIN32
    // For Windows, make the path POSIX compliant so code looking for directory
    // separators is simpler.  Just look for '/' not both '/' and '\\'
    std::replace_if(out_name.begin(), out_name.end(),
                    std::bind2nd(std::equal_to<char>(), '\\'), '/') ;
#endif
    size_t idir = out_name.rfind('/');
    if (idir != npos) {
        out_name = out_name.substr(idir+1, out_name.size());
    }
    size_t idot = out_name.rfind('.');
    if (idot != npos) {
        out_name = out_name.substr(0, idot) + ".xml";
    } else {
        out_name += ".xml";
    }
    std::ofstream out(out_name.c_str());
    out << xml;
}

std::string ct2ctml_string(const std::string& file)
{
#ifdef HAS_NO_PYTHON
    /*
     *  Section to bomb out if python is not
     *  present in the computation environment.
     */
    throw CanteraError("ct2ctml",
                       "python cti to ctml conversion requested for file, " + file +
                       ", but not available in this computational environment");
#endif

    string python_output, error_output;
    int python_exit_code;
    try {
        exec_stream_t python;
        python.set_wait_timeout(exec_stream_t::s_all, 1800000); // 30 minutes
        stringstream output_stream, error_stream;
        std::vector<string> args;
        args.push_back("-c");

        args.push_back(
                    "from __future__ import print_function\n"
                    "import sys\n"
                    "try:\n"
                    "    from cantera import ctml_writer\n"
                    "except ImportError:\n"
                    "    print('sys.path: ' + repr(sys.path) + '\\n', file=sys.stderr)\n"
                    "    raise\n"
                    "ctml_writer.convert(r'" + file + "', 'STDOUT')\n"
                    "sys.exit(0)\n");

        python.start(pypath(), args.begin(), args.end());
        std::string line;

        while (python.out().good()) {
            std::getline(python.out(), line);
            output_stream << line << std::endl;
        }

#ifdef _WIN32
        // Sleeping for 1 ms prevents a (somewhat inexplicable) deadlock while
        // reading from the stream.
        Sleep(1);
#endif
        while (python.err().good()) {
            std::getline(python.err(), line);
            error_stream << line << std::endl;
        }
        python.close();
        python_exit_code = python.exit_code();
        error_output = stripws(error_stream.str());
        python_output = output_stream.str();
    } catch (std::exception& err) {
        // Report failure to execute Python
        stringstream message;
        message << "Error executing python while converting input file:\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << err.what() << std::endl;
        throw CanteraError("ct2ctml_string", message.str());
    }

    if (python_exit_code != 0) {
        // Report a failure in the conversion process
        stringstream message;
        message << "Error converting input file \"" << file << "\" to CTML.\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << "The exit code was: " << python_exit_code << "\n";
        if (error_output.size() > 0) {
            message << "-------------- start of converter log --------------\n";
            message << error_output << std::endl;
            message << "--------------- end of converter log ---------------";
        } else {
            message << "The command did not produce any output." << endl;
        }
        throw CanteraError("ct2ctml_string", message.str());
    }

    if (error_output.size() > 0) {
        // Warn if there was any output from the conversion process
        stringstream message;
        message << "Warning: Unexpected output from CTI converter\n";
        message << "-------------- start of converter log --------------\n";
        message << error_output << std::endl;
        message << "--------------- end of converter log ---------------\n";
        writelog(message.str());
    }
    return python_output;
}


void ck2cti(const std::string& in_file, const std::string& thermo_file,
            const std::string& transport_file, const std::string& id_tag)
{
#ifdef HAS_NO_PYTHON
    /*
     *  Section to bomb out if python is not
     *  present in the computation environment.
     */
    string ppath = in_file;
    throw CanteraError("ct2ctml",
                       "python ck to cti conversion requested for file, " + ppath +
                       ", but not available in this computational environment");
#endif

    string python_output;
    int python_exit_code;
    try {
        exec_stream_t python;
        python.set_wait_timeout(exec_stream_t::s_all, 1800000); // 30 minutes
        python.start(pypath(), "-i");
        stringstream output_stream;

        ostream& pyin = python.in();
        pyin << "if True:\n" << // Use this so that the rest is a single block
                "    import sys\n" <<
                "    sys.stderr = sys.stdout\n" <<
                "    try:\n" <<
                "        from cantera import ck2cti\n" <<
                "    except ImportError:\n" <<
                "        print('sys.path: ' + repr(sys.path))\n" <<
                "        raise\n"
                "    ck2cti.Parser().convertMech(r'" << in_file << "',";
        if (thermo_file != "" && thermo_file != "-") {
            pyin << " thermoFile=r'" << thermo_file << "',";
        }
        if (transport_file != "" && transport_file != "-") {
            pyin << " transportFile=r'" << transport_file << "',";
        }
        pyin << " phaseName='" << id_tag << "',";
        pyin << " permissive=True,";
        pyin << " quiet=True)\n";
        pyin << "    sys.exit(0)\n\n";
        pyin << "sys.exit(7)\n";
        python.close_in();

        std::string line;
        while (python.out().good()) {
            std::getline(python.out(), line);
            output_stream << line << std::endl;;
        }
        python.close();
        python_exit_code = python.exit_code();
        python_output = stripws(output_stream.str());
    } catch (std::exception& err) {
        // Report failure to execute Python
        stringstream message;
        message << "Error executing python while converting input file:\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << err.what() << std::endl;
        throw CanteraError("ct2ctml", message.str());
    }

    if (python_exit_code != 0) {
        // Report a failure in the conversion process
        stringstream message;
        message << "Error converting input file \"" << in_file << "\" to CTI.\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << "The exit code was: " << python_exit_code << "\n";
        if (python_output.size() > 0) {
            message << "-------------- start of converter log --------------\n";
            message << python_output << std::endl;
            message << "--------------- end of converter log ---------------";
        } else {
            message << "The command did not produce any output." << endl;
        }
        throw CanteraError("ck2cti", message.str());
    }

    if (python_output.size() > 0) {
        // Warn if there was any output from the conversion process
        stringstream message;
        message << "Warning: Unexpected output from CTI converter\n";
        message << "-------------- start of converter log --------------\n";
        message << python_output << std::endl;
        message << "--------------- end of converter log ---------------\n";
        writelog(message.str());
    }
}

void get_CTML_Tree(Cantera::XML_Node* rootPtr, const std::string& file, const int debug)
{
    std::string ext = "";

    // find the input file on the Cantera search path
    std::string inname = findInputFile(file);
    writelog("Found file: "+inname+"\n", debug);

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
        string phase_xml = ctml::ct2ctml_string(inname);
        stringstream s(phase_xml);
        rootPtr->build(s);
        return;
    }

    writelog("Attempting to parse xml file " + inname + "\n", debug);
    ifstream fin(inname.c_str());
    if (!fin) {
        throw
        CanteraError("get_CTML_Tree",
                     "XML file " + inname + " not found");
    }
    rootPtr->build(fin);
    fin.close();
}

Cantera::XML_Node getCtmlTree(const std::string& file)
{
    Cantera::XML_Node root;
    get_CTML_Tree(&root, file);
    return root;
}

}
