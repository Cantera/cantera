/*
 * This is an example of a user-defined function that can be linked into Cantera.
 */

#include <iostream>

#include "ct_defs.h"
using namespace Cantera;

namespace User {

    void hello() {
        cout << "Hello!" << endl;
    }

}
