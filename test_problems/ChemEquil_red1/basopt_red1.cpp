#include "cantera/IdealGasMix.h"
#include "cantera/equilibrium.h"

using namespace std;
using namespace Cantera;

int main(int argc, char** argv)
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    try {
        IdealGasMix g("red1.xml", "gri30_mix");

#ifdef DEBUG_BASISOPTIMIZE
        Cantera::BasisOptimize_print_lvl = 0;
#endif
#ifdef DEBUG_CHEMEQUIL
        Cantera::ChemEquil_print_lvl = 0;
#endif

        double pres = 1.0E5;
        g.setState_TPX(2000.0, pres, "C2H2:0.9, CH:0.1");

        MultiPhase mphase;
        mphase.addPhase(&g, 10.0);
        mphase.init();
        int usedZeroedSpecies = 0;
        std::vector<size_t> orderVectorSpecies;
        std::vector<size_t> orderVectorElements;

        bool doFormMatrix = true;
        vector_fp formRxnMatrix;

        size_t nc = BasisOptimize(&usedZeroedSpecies, doFormMatrix,
                                  &mphase, orderVectorSpecies,
                                  orderVectorElements,
                                  formRxnMatrix);

        cout << "number of components = " << nc << endl;

        /*
         * The ChemEquil solver throws an error for this case.
         * The MultiPhaseEquil solver just gets the wrong result.
         */
        int it = equilibrate(g, "TP", -1);
        if (it != 1) {
            cerr << "incorrect number of iterations." << endl;
            return -1;
        }
        cout.unsetf(ios::floatfield);
        cout.precision(3);
        //cout << g;

        return 0;
    } catch (CanteraError& err) {
        std::cerr << err.what() << std::endl;
        cerr << "program terminating." << endl;
        return -1;
    }
}
