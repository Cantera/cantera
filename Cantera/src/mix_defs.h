#ifndef CT_MIX_DEFS_H
#define CT_MIX_DEFS_H

namespace Cantera {

    const int cNone = 0;

    // species thermo types
    const int cNASA = 1;
    const int cShomate = 2;

    // equation of state types
    const int cIdealGas = 1;
    const int cIncompressible = 2;
    const int cSurf = 3;

    // kinetic manager types
    const int cGasKinetics = 2;
    const int cGRI30 = 3;
    const int cInterfaceKinetics = 4;
}

#endif
