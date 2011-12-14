/**
 * @file makedsp.cpp
 *
 * Write a Visual Studio Project File with settings for Cantera.  This
 * program writes a project file to build a console application only.
 *
 */

// copyright 2001 California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void write_prog(string fname);

int main() {

    char buf[500];
    for (int j = 0; j < 500; j++) buf[j] = ' ';
    string line, projname, projdir, ctroot;

    // prompt for the project name
    cout << "Project name: ";
    cin >> projname;

    // prompt for the output directory
    cout << "Project directory (enter '.' for local directory): ";
    cin >> projdir;

    // get the Cantera root directory either from environment
    // variable CANTERA_ROOT (if set) or from user input

    const char* ctr;
    int iroot = 0;
    ctr = getenv("WIN_CANTERA_ROOT");
    if (ctr == 0) {
        ctr = getenv("CANTERA_ROOT");
        iroot = 1;
    }
    if (ctr != 0) {
        ctroot = ctr;
        cout << "\nCantera root directory: " << ctroot; 
        cout << (iroot == 1 ? " (CANTERA_ROOT)" : " (WIN_CANTERA_ROOT)")  << endl << endl;
    }
    else {
        iroot = -1;
        cout << "Cantera root directory: ";
        cin >> ctroot;
    }

#ifdef CVF
    ifstream proto((ctroot+"/tools/src/protocvf.dsp").c_str());
    if (!proto) {
        cout << "Error: file protocvf.dsp not found in " << ctroot << "/tools/src" << endl; 
        return 1; 
    }
#else
    ifstream proto((ctroot+"/tools/src/protocxx.dsp").c_str());
    if (!proto) {
        cout << "Error: file protocxx.dsp not found in " 
             << ctroot << "/tools/src" << endl; 
        return 1; 
    }
#endif

    string fname, progname;
    if (projdir != "." && projdir != "'.'") { 
        fname = projdir+"\\"+projname+".dsp";
        progname = projdir+"\\"+projname+".cpp";
    }
    else {
        fname = string(projname)+".dsp";
        progname = string(projname)+".cpp";
    }
    ofstream fout(fname.c_str());
    while (!proto.eof()) {
        proto.getline( buf,500);
        line = buf;
        int i;
        while (i = line.find("__PROJECT__"), i >=0) {
            line.replace(i,11,projname);
        }
        while (i = line.find("__CTROOT__"), i >=0) {
            line.replace(i,10,ctroot);
        }
        fout << line << endl;
    }
    proto.close();
    fout.close();

    // if main program file doesn't exist, create a prototype
    ifstream fprog(progname.c_str());
    bool wrote_prog = false;
    if (!fprog) {
        write_prog(progname);
        wrote_prog = true;
    }
    else
        fprog.close();

    cout << "created Developer Studio project file " 
         << fname;
    if (wrote_prog) 
        cout << " and main program file " << progname;
    cout << endl; 
    return 0;
}

void write_prog(string fname) {

    ofstream f(fname.c_str());
    f << endl
      << "#include \"Cantera.h\"" << endl
      << "// include additional header files here if needed" << endl
      << endl
      << "int main(int argc, char** argv) {" << endl
      << "  try {" << endl
      << "    // your code goes here" << endl
      << "    return 0;" << endl
      << "  }" << endl
      << "  catch (CanteraError) {" << endl
      << "    Application::showErrors(cerr);" << endl
      << "    cerr << \"program terminating.\" << endl;" << endl
      << "    return -1;" << endl
      << "  }" << endl
      << "}" << endl;
    f.close();
}


    

