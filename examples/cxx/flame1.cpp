
#include "/Applications/Cantera/include/cantera/Cantera.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport.h"

int main() {

    // create the gas object

    IdealGasMix gas("gri30.cti");
    doublereal temp = 500.0;
    doublereal pres = 2.0*OneAtm;
    gas.setState_TPX(temp, pres, "CH4:1.0, O2:2.0, N2:7.52");

    // create a transport manager that implements
    // mixture-averaged transport properties

    Transport* tr = newTransportMgr("Mix", &gas);


    //=============  build each domain ========================


    //-------- step 1: create the stagnation flow -------------

    StFlow flow(&gas);

    // create an initial grid
    doublereal z[] = {0.0, 0.05, 0.1, 0.15, 0.2};
    flow.setupGrid(5, z);

    // specify the objects to use to compute kinetic rates and 
    // transport properties
    flow.setKinetics(&gas);
    flow.setTransport(&tr);

    flow.setPressure(0.05*OneAtm);



    //------- step 2: create the inlet  -----------------------

    Inlet1D inlet;
    inlet.setMoleFractions("CH4:1, O2:2, N2:7.52");
    inlet.setMdot(0.1);


    //------- step 3: create the surface  ---------------------

    Surf1D surf;

    //=================== create the container and insert the domains =====
    
    vector<Domain1D*> domains;
    domains.push_back(inlet);
    domains.push_back(flow);
    domains.push_back(surf);

    OneDim flamesim(domains);
}
