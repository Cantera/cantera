/*
 * $Id$
 */
#include "WaterPropsIAPWSphi.h"
#include <new>

#include <cstdio>

using namespace std;

int main () {

    WaterPropsIAPWSphi *phi = new WaterPropsIAPWSphi();

    phi->check1();
    phi->check2();

    delete phi;
    return 0;
}
