/*
 * $Id$
 */
#include "stdio.h"
#include "WaterPropsIAPWSphi.h"
#include <new>
using namespace std;

int main () {

    WaterPropsIAPWSphi *phi = new WaterPropsIAPWSphi();

    phi->check1();
    phi->check2();

    delete phi;
    return 0;
}
