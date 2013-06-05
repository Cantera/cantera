/**
 *  @file statmech
 *       test problem for statistical mechanics in cantera
 */

//  Example
//
// Test case for the statistical mechanics in cantera
//

#include "cantera/transport.h"
#include "cantera/IdealGasMix.h"
#include "cantera/equil/equil.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{

    try {
        int k;
        IdealGasMix g("test_stat_trans.xml", "example");
        int nsp = g.nSpecies();
        double pres = 1.0E5;

        vector_fp Xset(nsp, 0.0);
        Xset[0] =  0.5 ;
        Xset[1] =  0.5;

        g.setState_TPX(1500.0, pres, DATA_PTR(Xset));
        equilibrate(g, "TP", -1);

        // init pecos transport
        int log_level = 0;
        Transport* tran = newTransportMgr("Pecos", &g, log_level=0);
        PecosTransport* tranMix = dynamic_cast<PecosTransport*>(tran);

        cout << "here";


        vector_fp cp_R(nsp, 0.0);
        g.getCp_R(DATA_PTR(cp_R));

        //for(int i=0;i<nsp;i++)
        //{
        //  std::cout.precision(10);
        //  std::cout << cp_R[i] << std::endl;
        //	}

        // error check-- exactly 2.5 for atoms
        if (cp_R[0] != 2.5) {
            std::cout << "Error for monotomic Species!\n";
            return 1;
        }

        // error check: analytical result is more complicated for
        // molecules. One species should suffice, lets try NO2, with
        // three vibrational modes:
        /// theta[0]: 1.07900e3
        /// theta[1]: 1.90000e3
        /// theta[2]: 2.32700e3
        // at T = 1500
        //
        // This is precisely: 6.655804161 (e.g. 5/2 + 2 + 3.1558..)
        //
        double theta[3];
        theta[0] = 1.07900e3;
        theta[1] = 1.90000e3;
        theta[2] = 2.32700e3;

        double T;
        T = 1500.0;

        double denom;
        double ctr = 0.0;
        double GasConstant = 1.0;

        for (int i = 0; i < 3; i++) {
            denom = exp(2*theta[i]/T) - 2* exp(theta[i]/T) + 1;
            ctr += GasConstant * theta[i] * (theta[i] * exp(theta[i]/T)/(T*T))/ (denom);
            //std::cout << "survey says: " << ctr << " and denom is: " << denom << std::endl;
        }
        //std::cout << "survey says: " << ctr << " and denom is: " << denom << std::endl;
        double sol = ctr + 5/2 + 2;
        double tol = 1e-9;

        if (abs(cp_R[3] - sol) >= tol) {
            double diff = cp_R[3]-sol;
            std::cout << "Error for Species NO2!\n";
            std::cout << "Diff was: " << diff << "\n";
            return 1;
        }

    } catch (CanteraError) {
        showErrors(cout);
        return 1;
    }

    // Mark it zero!
    return 0;

}
