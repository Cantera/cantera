import numpy as np

def compareTimeSeries(reference, sample, rtol=1e-5, atol=1e-12):
    """
    Compare two 2D arrays of time series data. Each data set should contain
    the time in the first column and data series to be compared in successive
    columns.

    The times in each data set do not need to be the same: The data from the
    second data set will be interpolated onto the times in the first data set
    before being compared. This means that the time range of the "sample" data
    set should be at least as long as the "reference" data set.

    After interpolation, each data point must satisfy either the relative
    error tolerance specified by `rtol` or the absolute error tolerance
    specified by `atol`.

    If the comparison succeeds, this function returns `None`. If the comparison
    fails, a formatted report of the differing elements is returned.
    """
    if isinstance(reference, str):
        reference = np.genfromtxt(reference, delimiter=',').T
    else:
        reference = np.asarray(reference).T

    sample = np.asarray(sample).T

    assert reference.shape[0] == sample.shape[0]

    nVars = reference.shape[0]
    nTimes = reference.shape[1]

    bad = []
    template = '{0:9.4e}  {1: 3d}   {2:14.8e}  {3:14.8e}  {4:9.3e}  {5:9.3e}'
    for i in range(1, nVars):
        comp = np.interp(reference[0], sample[0], sample[i])
        for j in range(nTimes):
            a = comp[j]
            b = reference[i,j]
            abserr = abs(a-b)
            relerr = 2 * abs(a-b) / (abs(a) + abs(b) + atol)
            if abserr > atol and relerr > rtol:
                bad.append(template.format(reference[0][j], i, a, b, abserr, relerr))

    if bad:
        header = ['Failed time series comparisons:',
                  '  time   component  reference       test value    abs. err.  rel. err',
                  '----------  ---   --------------  --------------  ---------  ---------']
        return '\n'.join(header + bad)
    else:
        return None
