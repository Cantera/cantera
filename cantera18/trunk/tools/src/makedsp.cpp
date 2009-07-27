/**
 * @file makedsp.cpp
 *
 * Write a Win32 Visual Studio Project File with settings for Cantera
 *
 */

// copyright (c) 2001 California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main() {

    char buf[500];
    for (int j = 0; j < 500; j++) buf[j] = ' ';
    string line, projname, ctroot;

    // prompt for the project name
    cout << "Project name: ";
    cin >> projname;


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

    ifstream proto((ctroot+"/tools/src/proto.dsp").c_str());
    if (!proto) {
        cout << "Error: file proto.dsp not found in " << ctroot << "/tools/src" << endl; 
        return 1; 
    }

    ofstream fout((projname+".dsp").c_str());
    while (!proto.eof()) {
        proto.getline( buf,500);
        line = buf;
        int i;
        while (i = line.find("__PROTO__"), i >=0) {
            line.replace(i,9,projname);
        }
        while (i = line.find("__CANTERAROOT__"), i >=0) {
            line.replace(i,15,ctroot);
        }
        fout << line << endl;
    }
    proto.close();
    fout.close();

    cout << "created Visual Studio project file " << projname << ".dsp" << endl;  
    return 0;
}


    

