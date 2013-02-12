/**
 *  @file mixGasTransport.cpp
 *       test problem for mixture transport
 */

//  Example
//
// Test case for mixture transport in a gas
// The basic idea is to set up a gradient of some kind.
// Then the resulting transport coefficients out.
// Essentially all of the interface routines should be
// exercised and the results dumped out.
//
// A blessed solution test will make sure that the actual
// solution doesn't change as a function of time or
// further development.

// perhaps, later, an analytical solution could be added

#include "cantera/transport.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/thermo.h"

#include <iostream>
#include <string>
#include <vector>
#include <string>
#include <iomanip>

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    string infile = "diamond.xml";

    try {

        Cantera::setDefaultIndependentVars(INDVAR_TP_MOLEFRACTION_VECTOR);

	IdealGasPhase<doubleFAD> * ig_fad = new IdealGasPhase<doubleFAD>("diamond.xml", "");

	delete ig_fad;


    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return -1;
    }

    return 0;
}
/***********************************************************/
