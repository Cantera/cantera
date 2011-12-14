
#include <cantera/Cantera.h>
using namespace Cantera;
using namespace std;

//#include "ctexceptions.h"

// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#define NUM_EXAMPLES 6

int kinetics_example1(int job);
int kinetics_example2(int job);
int kinetics_example3(int job);
int transport_example1(int job);
int transport_example2(int job);
int equil_example1(int job);
int rxnpath_example1(int job);

typedef int (*exfun)(int n);
 
int run_example(int n, exfun f, int job = 2) {
    cout << "\n\n\n\n>>>>>  example " << n+1 << "\n\nDescription:  " << endl;
    int i = f(job);
    return i;
}
   
// array of example functions   
exfun fex[] = {kinetics_example1, kinetics_example2, kinetics_example3,
               equil_example1, 
               transport_example1, transport_example2}; //, rxnpath_example1};
 

// main program
int main(int argc, char** argv) {

    int example_num = 0;
    cout << endl
         << "-----------------------------------" << endl
         << "       Cantera C++ examples        " << endl
         << "-----------------------------------" << endl;

    if (argc > 1) {
        string v1 = string(argv[1]);
        if (v1 == "-h") {
            cout << "\nusage: examples <n> \n\nwhere <n> "
                "is the example number to run." << endl;
            cout << "if <n> is omitted, all examples are run." 
                 << endl << endl;
            cout << "Examples: " << endl;
            for (int i = 0; i < NUM_EXAMPLES; i++) {
                cout << "(" << i+1 << ") ";
                fex[i](0);
                cout << endl;
            }
            exit(0);
        }
        else
            example_num = atoi(argv[1]);
    }
    try {

        int i = 0;
        if (example_num == 0) {
            int j;
            for (j = 0; j < NUM_EXAMPLES; j++) {
                i = run_example(j, fex[j], 2);
            }
        }
        else if (example_num > 0 && example_num <= NUM_EXAMPLES)
            i = run_example(example_num-1, fex[example_num-1], 2);

        return 0;
    }
    catch (CanteraError) {
        showErrors(cerr);
        cerr << "program terminating." << endl;
        return -1;
    }
}
