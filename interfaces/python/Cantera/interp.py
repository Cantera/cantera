
def interp(z0, z, f):
    """
    Linear interpolation.

    Sequences z and f must be of the same length,
    and the entries in z must be monotonically increasing.

    Example:
    >>> z = [0.0, 0.2, 0.5, 1.2, 2.1]
    >>> f = [3.0, 2.0, 1.0, 0.0, -1.0]
    >>>print interp(-2, z, f), interp(0.5, z, f), interp(6, z, f)
    3.0 7.0 -9.0
    """

    n = len(z)

    # if z0 is outside the range of z, then return the endpoint value,
    # instead of extrapolating.
    if z0 <= z[0]:
        return f[0]
    elif z0 > z[-1]:
        return f[-1]

    for i in range(1,n):
        if z0 <= z[i]:
            return f[i-1] + (f[i] - f[i-1])*(z0 - z[i-1])/(z[i] - z[i-1])

    # if this statement is reached, then there is an error.
    raise 'interpolation error!'


def quadInterp(z0, z, f):

    n = len(z)

    # if z0 is outside the range of z, then return the endpoint value,
    # instead of extrapolating.
    if z0 <= z[0]:
        return f[0]
    elif z0 > z[-1]:
        return f[-1]

    for i in range(1,n):
        if z0 <= z[i]:
            j = max(2,i)
            dx21 = z[j-1] - z[j-2]
            dx32 = z[j] - z[j-1]
            dx31 = dx21 + dx32
            dy32 = f[j] - f[j-1]
            dy21 = f[j-1] - f[j-2]
            a = (dx21*dy32 - dy21*dx32)/(dx21*dx31*dx32)
            return a*(z0 - z[j-2])*(z0 - z[j-1]) + (
                (dy21/dx21)*(z0 - z[j-1]) + f[j-1])

    # if this statement is reached, then there is an error.
    raise 'interpolation error!'

## if __name__ == '__main__':
##     z = [0.0, 0.2, 0.5, 1.2, 2.1]
##     f = [3.0, 5.0, 11.0, 0.0, -9.0]

##     print interp(-2, z, f), interp(0.3, z, f), interp(6, z, f)
##     print quadInterp(-2, z, f), quadInterp(0.3, z, f), quadInterp(6, z, f)
