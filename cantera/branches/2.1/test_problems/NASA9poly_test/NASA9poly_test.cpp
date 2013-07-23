/**
 *  @file NASA9poly_test
 *       test problem for NASA 9 coefficient formulation
 */
#include "cantera/transport.h"
#include "cantera/IdealGasMix.h"
#include "cantera/transport/TransportFactory.h"

#include <memory>
#include <cstdio>

using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        IdealGasMix g("gasNASA9.xml", "gri30_mix");
        size_t nsp = g.nSpecies();
        double pres = 1.0E5;
        vector_fp Xset(nsp, 0.0);
        Xset[0] = 0.5 ;
        Xset[1] = 0.5;

        g.setState_TPX(1500.0, pres, &Xset[0]);

        vector_fp cp_R(nsp, 0.0);
        vector_fp H_RT(nsp, 0.0);
        vector_fp S_R(nsp, 0.0);

        g.getEntropy_R(&S_R[0]);
        g.getCp_R(&cp_R[0]);
        g.getEnthalpy_RT(&H_RT[0]);

        printf("Comparisons of H2 calculated via several equivalent classes:\n");
        printf("1500 K and 1 atm:\n");
        printf("         NasaThermo   Nasa9   Nasa9_4reg \n");
        printf("  cp/R: %11.6g %11.6g %11.6g\n", cp_R[0], cp_R[1], cp_R[2]);
        printf("  H/RT: %11.6g %11.6g %11.6g\n", H_RT[0], H_RT[1], H_RT[2]);
        printf("  S/R: %11.6g %11.6g %11.6g\n", S_R[0], S_R[1], S_R[2]);

        std::auto_ptr<Transport> tran(newTransportMgr("Mix", &g));
        vector_fp Gvalues(nsp, 0.0);

        printf("Viscosity and thermal Cond vs. T\n");
        for (int k = 0; k < 40; k++) {
            double T1 = 400. + 200. * k;
            g.setState_TPX(T1, pres, &Xset[0]);
            g.getPureGibbs(&Gvalues[0]);
            double visc = tran->viscosity();
            double cond = tran->thermalConductivity();
            printf("    %13g %13.5g %13.5g\n", T1, visc, cond);
        }
    } catch (CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
    return 0;
}
