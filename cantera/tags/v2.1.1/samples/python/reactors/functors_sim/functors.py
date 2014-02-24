# This example shows how to create 'functors' - objects that evaluate
# functions. These are useful for specifying the expansion rate or
# heat flux at a wall.

from Cantera.Func import *

# create f1(t) = 4 + 6t + 8t^2 + t^3
f1 = Polynomial([4.0, 6.0, 8.0, 1.0])

# create sin(t)
f2 = Fourier(1.0, [(0.0, 0.0), (0.0, 1.0)])

# functors can be combined by +,*,or / to create new functors
f3 = f2*f2

xpts = 0.1*array(range(100))

for x in xpts:
    print x, f1(x), f2(x), f3(x)
