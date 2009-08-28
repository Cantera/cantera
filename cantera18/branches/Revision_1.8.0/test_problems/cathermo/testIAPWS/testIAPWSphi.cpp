/*
 * $Id: testIAPWSphi.cpp,v 1.2 2008/12/17 17:31:13 hkmoffa Exp $
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
