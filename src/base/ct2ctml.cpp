/**
 * @file ct2ctml.cpp
 * Driver for the system call to the python executable that converts
 * cti files to ctml files (see \ref inputfiles).
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/base/ctml.h"
#include "cantera/base/stringUtils.h"
#include "../../ext/libexecstream/exec-stream.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <functional>

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std;

namespace Cantera
{

//! return the full path to the Python interpreter.
/*!
 * Use the environment variable PYTHON_CMD if it is set. If not, return
 * the string 'python'.
 *
 * Note, there are hidden problems here that really direct us to use a full
 * pathname for the location of python. Basically the system call will use the
 * shell /bin/sh, in order to launch python. This default shell may not be the
 * shell that the user is employing. Therefore, the default path to python may
 * be different during a system call than during the default user shell
 * environment. This is quite a headache. The answer is to always set the
 * PYTHON_CMD environmental variable in the user environment to an absolute path
 * to locate the python executable. Then this issue goes away.
 */
static string pypath()
{
    string s = "python";
    const char* py = getenv("PYTHON_CMD");

    if (py) {
        string sp = trimCopy(string(py));
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
                    std::bind2nd(std::equal_to<char>(), '\\'), '/');
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
    std::ofstream out(out_name);
    out << xml;
}

static std::string call_ctml_writer(const std::string& text, bool isfile)
{
    std::string file, arg;

    if (isfile) {
        file = text;
        arg = "r'" + text + "'";
    } else {
        file = "<string>";
        arg = "text=r'''" + text + "'''";
    }

    string python_output, error_output;
    int python_exit_code;
    try {
        exec_stream_t python;
        python.set_wait_timeout(exec_stream_t::s_all, 1800000); // 30 minutes
        stringstream output_stream, error_stream;
        python.start(pypath(), "");
        ostream& pyin = python.in();

        pyin << "from __future__ import print_function\n"
                "if True:\n"
                "    import sys\n"
                "    try:\n"
                "        from cantera import ctml_writer\n"
                "    except ImportError:\n"
                "        print('sys.path: ' + repr(sys.path) + '\\n', file=sys.stderr)\n"
                "        raise\n"
                "    ctml_writer.convert(";
        pyin << arg << ", outName='STDOUT')\n";
        pyin << "    sys.exit(0)\n\n";
        pyin << "sys.exit(7)\n";

        python.close_in();
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
        error_output = trimCopy(error_stream.str());
        python_output = output_stream.str();
    } catch (std::exception& err) {
        // Report failure to execute Python
        stringstream message;
        message << "Error executing python while converting input file:\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << err.what() << std::endl;
        throw CanteraError("call_ctml_writer", message.str());
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
        throw CanteraError("call_ctml_writer", message.str());
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

std::string ct2ctml_string(const std::string& file)
{
    return call_ctml_writer(file, true);
}

std::string ct_string2ctml_string(const std::string& cti)
{
    return call_ctml_writer(cti, false);
}

void ck2cti(const std::string& in_file, const std::string& thermo_file,
            const std::string& transport_file, const std::string& id_tag)
{
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
                "    _ = ck2cti.convertMech(r'" << in_file << "',";
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
        python_output = trimCopy(output_stream.str());
    } catch (std::exception& err) {
        // Report failure to execute Python
        stringstream message;
        message << "Error executing python while converting input file:\n";
        message << "Python command was: '" << pypath() << "'\n";
        message << err.what() << std::endl;
        throw CanteraError("ck2cti", message.str());
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
        message << "Warning: Unexpected output from ck2cti\n";
        message << "-------------- start of converter log --------------\n";
        message << python_output << std::endl;
        message << "--------------- end of converter log ---------------\n";
        writelog(message.str());
    }
}

}
