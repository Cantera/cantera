/*
 * $Id: testIAPWSphi.cpp,v 1.1 2006/07/04 00:21:28 hkmoffa Exp $
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
