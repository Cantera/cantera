#include "cantera/thermo/WaterPropsIAPWSphi.h"
#include <new>

#include <cstdio>

using namespace std;
using namespace Cantera;

int main()
{
#ifdef _MSC_VER
    _set_output_format(_TWO_DIGIT_EXPONENT);
#endif
    WaterPropsIAPWSphi* phi = new WaterPropsIAPWSphi();

    phi->check1();
    phi->check2();

    delete phi;
    return 0;
}
