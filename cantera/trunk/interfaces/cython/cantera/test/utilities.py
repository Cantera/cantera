import numpy as np
import unittest


class CanteraTest(unittest.TestCase):
    def assertNear(self, a, b, rtol=1e-8, atol=1e-12, msg=None):
        cmp = 2 * abs(a - b)/(abs(a) + abs(b) + 2 * atol / rtol)
        if cmp > rtol:
            message = ('AssertNear: %.14g - %.14g = %.14g\n' % (a, b, a-b) +
                       'Relative error of %10e exceeds rtol = %10e' % (cmp, rtol))
            if msg:
                message = msg + '\n' + message
            self.fail(message)

    def assertArrayNear(self, A, B, rtol=1e-8, atol=1e-12, msg=None):
        if len(A) != len(B):
            self.fail("Arrays are of different lengths ({0}, {1})".format(len(A), len(B)))
        A = np.asarray(A)
        B = np.asarray(B)

        for a,b in zip(A.flat, B.flat):
            self.assertNear(a,b, rtol, atol, msg)


def compareProfiles(reference, sample, rtol=1e-5, atol=1e-12, xtol=1e-5):
    """
    Compare two 2D arrays of spatial or time profiles. Each data set should
    contain the time or space coordinate in the first column and data series
    to be compared in successive columns.

    The coordinates in each data set do not need to be the same: The data from
    the second data set will be interpolated onto the coordinates in the first
    data set before being compared. This means that the range of the "sample"
    data set should be at least as long as the "reference" data set.

    After interpolation, each data point must satisfy a combined relative and absolute
    error criterion specified by `rtol` and `atol`.

    If the comparison succeeds, this function returns `None`. If the comparison
    fails, a formatted report of the differing elements is returned.
    """
    if isinstance(reference, str):
        reference = np.genfromtxt(reference, delimiter=',').T
    else:
        reference = np.asarray(reference).T

    if isinstance(sample, str):
        sample = np.genfromtxt(sample, delimiter=',').T
    else:
        sample = np.asarray(sample).T

    assert reference.shape[0] == sample.shape[0]

    nVars = reference.shape[0]
    nTimes = reference.shape[1]

    bad = []
    template = '{0:9.4e}  {1: 3d}   {2:14.7e}  {3:14.7e}  {4:9.3e}  {5:9.3e}  {6:9.3e}'
    for i in range(1, nVars):
        scale = max(max(abs(reference[i])), reference[i].ptp(),
                    max(abs(sample[i])), sample[i].ptp())
        slope = np.zeros(nTimes)
        slope[1:] = np.diff(reference[i]) / np.diff(reference[0]) * reference[0].ptp()

        comp = np.interp(reference[0], sample[0], sample[i])
        for j in range(nTimes):
            a = reference[i,j]
            b = comp[j]
            abserr = abs(a-b)
            relerr = abs(a-b) / (scale + atol)

            # error that can be accounted for by shifting the profile along
            # the time / spatial coordinate
            xerr = abserr / (abs(slope[j]) + atol)

            if abserr > atol and relerr > rtol and xerr > xtol:
                bad.append(template.format(reference[0][j], i, a, b, abserr, relerr, xerr))

    if bad:
        header = ['Failed series comparisons:',
                  'coordinate  comp.  reference val    test value    abs. err   rel. err    pos. err',
                  '----------  ---   --------------  --------------  ---------  ---------  ---------']
        return '\n'.join(header + bad)
    else:
        return None
