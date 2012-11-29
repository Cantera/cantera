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
            self.fail("Arrays are of different lengths ({0}, {0})".format(len(A), len(B)))
        A = np.asarray(A)
        B = np.asarray(B)

        for a,b in zip(A.flat, B.flat):
            self.assertNear(a,b, rtol, atol, msg)
