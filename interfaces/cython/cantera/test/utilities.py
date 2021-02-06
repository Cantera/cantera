import numpy as np
import sys
import os
import warnings
import shutil
import tempfile
import unittest
import errno
import cantera


class CanteraTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Create a working directory for output files. If this is
        # an in-source test, create the directory in the root
        # test/work directory. Otherwise, create a system level
        # temporary directory
        root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                '..', '..', '..', '..'))
        if os.path.exists(os.path.join(root_dir, 'SConstruct')):
            cls.test_work_dir = os.path.join(root_dir, 'test', 'work', 'python')
            try:
                os.makedirs(cls.test_work_dir)
            except OSError as e:
                if e.errno == errno.EEXIST:
                    pass
                elif e.errno == errno.EACCES:
                    cls.test_work_dir = tempfile.mkdtemp()
                else:
                    raise
        else:
            cls.test_work_dir = tempfile.mkdtemp()

        cantera.make_deprecation_warnings_fatal()
        cantera.add_directory(cls.test_work_dir)
        cls.test_data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'data'))
        cls.cantera_data = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', 'data'))

    @classmethod
    def tearDownClass(cls):
        # Remove the working directory after testing, but only if its a temp directory
        if tempfile.tempdir is not None:
            try:
                shutil.rmtree(cls.test_work_dir)
            except OSError:
                pass

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

        worst = 0, ''
        for i in np.ndindex(A.shape):
            a = A[i]
            b = B[i]
            cmp = 2 * abs(a - b)/(abs(a) + abs(b) + 2 * atol / rtol)
            if cmp > rtol:
                message = ('AssertNear: {:.14g} - {:.14g} = {:.14g}\n'.format(a, b, a-b) +
                           'Relative error for element {} of {:10e} exceeds rtol = {:10e}'.format(i, cmp, rtol))
                if msg:
                    message = msg + '\n' + message
                if cmp > worst[0]:
                    worst = cmp, message

        if worst[0]:
            self.fail(worst[1])

    def compare(self, data, reference_file, rtol=1e-8, atol=1e-12):
        """
        Compare an array with a reference data file, or generate the reference
        file if it does not exist.
        """
        data = np.array(data)
        if os.path.exists(reference_file):
            # Compare with existing output file
            ref = np.genfromtxt(reference_file)
            self.assertEqual(data.shape, ref.shape)
            for i in range(ref.shape[0]):
                self.assertArrayNear(ref[i], data[i], rtol, atol)
        else:
            # Generate the output file for the first time
            warnings.warn('Generating test data file:' +
                          os.path.abspath(reference_file))
            np.savetxt(reference_file, data, fmt='%.10e')


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
                bad.append((reference[0][j], i, a, b, abserr, relerr, xerr))

    footer = []
    maxrows = 10
    if len(bad) > maxrows:
        bad.sort(key=lambda row: -row[5])
        footer += ['Plus {0} more points exceeding error thresholds.'.format(len(bad)-maxrows)]
        bad = bad[:maxrows]

    if bad:
        header = ['Failed series comparisons:',
                  'coordinate  comp.  reference val    test value    abs. err   rel. err    pos. err',
                  '----------  ---   --------------  --------------  ---------  ---------  ---------']
        return '\n'.join(header + [template.format(*row) for row in bad] + footer)
    else:
        return None
