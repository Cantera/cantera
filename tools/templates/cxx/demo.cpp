//
//     Replace this sample main program with your program
//
//

#include "Cantera.h"
#include "IdealGasMix.h"
using namespace Cantera;

main() {

    IdealGasMix gas("h2o2.xml");
    gas.setState_TPX_String(1200.0, OneAtm, 
        "H2:2, O2:1, OH:0.01, H:0.01, O:0.01");
    equilibrate(gas,"HP");

    printf("**** C++ Test Program ****\n\n");
    printf(
        "Temperature:    14.5g K\n"
        "Pressure:       14.5g Pa\n"
        "Density:        14.5g kg/m3\n"
        "Molar Enthalpy: 14.5g J/kmol\n"
        "Molar Entropy:  14.5g J/kmol-K\n"
        "Molar cp:       14.5g J/kmol-K\n",
        gas.temperature(), gas.pressure(), gas.density(),
        gas.enthalpy_mole(), gas.entropy_mole(), gas.cp_mole());



// c
// c     Reaction information
// c
//       irxns = nReactions()
//       call getFwdRatesOfProgress(qf)
//       call getRevRatesOfProgress(qr)
//       call getNetRatesOfProgress(q)
//       do i = 1,irxns
//          call getReactionEqn(i,eq)
//          write(*,20) eq,qf(i),qr(i),q(i)
//  20      format(a20,3g14.5,' kmol/m3/s')
//       end do
//       stop
//       end
}
      
