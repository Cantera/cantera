//
// main program to run all C++ demos.
//
#include <cantera/Cantera.h>

#define CXX_DEMO
#include "rankine.cpp"
#include "flamespeed.cpp"
#include "kinetics1.cpp"

#include <time.h>

typedef int (*exfun)(int, void*);
 
// array of demo functions   
exfun fex[] = {kinetics1, openRankine, flamespeed};
 
string demostr[] = {"zero-D kinetics        ", 
                    "open Rankine cycle     ",
                    "flamespeed             "};

int np[] = {0, 0, 1};
double p[] = {0, 0, 0.9};

#define NDEMOS 3

int mainmenu() {
    int i, idemo;
    cout << "C++ Demo Programs " << endl;
    for (i = 0; i < NDEMOS; i++) {
        cout << "     " << i+1 << ") " << demostr[i] << endl;
    }
    cout << "     " << i+1 << ")" << " run all demos" << endl;
    cout << "Enter demo number (or 0 to quit): ";

    cin >> idemo;
    if (idemo <= 0) return -99;

    int iout = 0;
    try {
        if (idemo > 0 && idemo < NDEMOS+1) {
            clock_t t0 = clock();
            iout = fex[idemo-1](0, 0);
            clock_t t1 = clock();
            cout << endl << "elapsed time: " 
                 << 1.0*(t1 - t0)/CLOCKS_PER_SEC << " s " << endl;
            return idemo;
        }
        else if (idemo >= NDEMOS+1) {
            clock_t t0 = clock();
            for (i = 0; i < NDEMOS; i++) 
                iout = fex[i](np[i], &p[i]);
            clock_t t1 = clock();
            cout << "time: " << 1.0*(t1 - t0)/CLOCKS_PER_SEC << " s " << endl;
            return iout;
        }
        return 0;
    }
    catch (CanteraError) {
        showErrors(cerr);
        return -1;
    }
}

int main() {
    int i;
    while (1 > 0) {
        i = mainmenu();
        appdelete();
        if (i == -99) break;
    }
    return 0;
}
