//
// main program to run all C++ demos.
//
#include <cantera/Cantera.h>

#define CXX_DEMO
#include "rankine.cpp"
#include "flamespeed.cpp"
#include "kinetics1.cpp"

typedef int (*exfun)();
 
int run_example(int n, exfun f, int job = 2) {
    cout << "\n\n\n\n>>>>>  example " << n+1 << "\n\nDescription:  ";
    int i = f(job);
    showErrors(cout);
    return i;
}
   
// array of demo functions   
exfun fex[] = {kinetics1, rankine, flamespeed};
 
string demostr[] = {"zero-D kinetics        ", 
                    "open Rankine cycle     ",
                    "flamespeed             "};

#define NDEMOS 3

main() {
    int i, idemo;
    for (i = 0; i < NDEMOS; i++) {
        cout << i << ") " << demostr[i] << endl;
    }
    cout << i << ")" << "run all demos" << endl;
    cout << "Enter demo number: " << endl;

    cin >> idemo;

    try {
        if (idemo > 0 && idemo < NDEMOS) {
            return fex[idemo]();
        }
        else if (idemo == NDEMOS) {
            int iout = 0;
            for (i = 0; i < NDEMOS; i++) 
                iout = fex[i]();
            return iout;
        }
    }
    catch (CanteraError) {
        showErrors(cerr);
        appDelete();
        return -1;
    }
}
