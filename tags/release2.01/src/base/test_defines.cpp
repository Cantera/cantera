#include "CanteraDefines.h"

int main() {

  Cantera::Real x;
  Cantera::array_fp a(7);
  int na = a.size();
  int n;

  for (n = 0; n < na; n++) {
    std::cout << n << "  " <<  na << "  " << a[n] << std::endl;
  }
		   

  
}
