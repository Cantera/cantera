#include "subs.h"
#include "utils.h"

namespace tpx {


Substance * GetSub(int isub) {
 	if (isub == 0)
                return new water;
 	else if (isub == 1)
 		return new nitrogen;
 	else if (isub == 2)
 		return new methane;
 	else if (isub == 3)
 		return new hydrogen;
 	else if (isub == 4)
 		return new oxygen;
// 	else if (isub == 5)
// 		return new HFC134a;
  	else
            return 0;
	}

}
