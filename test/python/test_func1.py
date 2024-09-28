import math
import numpy as np
import pytest

import cantera as ct

from .utilities import (
    assertNear,
    assertArrayNear
)

class TestFunc1:
    def test_function(self):
        f = ct.Func1(np.sin)
        assertNear(f(0), np.sin(0))
        assertNear(f(0.1), np.sin(0.1))
        assertNear(f(0.7), np.sin(0.7))

    def test_lambda(self):
        f = ct.Func1(lambda t: np.sin(t)*np.sqrt(t))
        assert f.type == "functor"
        for t in [0.1, 0.7, 4.5]:
            assertNear(f(t), np.sin(t)*np.sqrt(t))

    def test_callable(self):
        class Multiplier:
            def __init__(self, factor):
                self.factor = factor
            def __call__(self, t):
                return self.factor * t

        m = Multiplier(8.1)
        f = ct.Func1(m)
        assert f.type == "functor"
        for t in [0.1, 0.7, 4.5]:
            assertNear(f(t), 8.1*t)

    def test_constant(self):
        f = ct.Func1(5)
        for t in [0.1, 0.7, 4.5]:
            assertNear(f(t), 5)
        assert f.type == "constant"

    def test_sequence(self):
        f = ct.Func1([5])
        for t in [0.1, 0.7, 4.5]:
            assertNear(f(t), 5)

        with pytest.raises(TypeError):
            ct.Func1([3,4])

    def test_numpy(self):
        f = ct.Func1(np.array(5))
        assert f.type == "constant"
        g = ct.Func1(np.array([[5]]))
        assert g.type == "constant"
        for t in [0.1, 0.7, 4.5]:
            assertNear(f(t), 5)
            assertNear(g(t), 5)

        with pytest.raises(TypeError):
            ct.Func1(np.array([3,4]))

    def test_failure(self):
        def fails(t):
            raise ValueError('bad')

        f = ct.Func1(fails)
        with pytest.raises(ValueError):
            f(0.1)

    def test_unpicklable(self):
        import pickle
        f = ct.Func1(np.sin)
        with pytest.raises(NotImplementedError):
            pickle.dumps(f)

    def test_uncopyable(self):
        import copy
        f = ct.Func1(np.sin)
        with pytest.raises(NotImplementedError):
            copy.copy(f)

    def test_simple(self):
        with pytest.raises(NotImplementedError):
            ct.Func1("spam")
        functors = {
            'sin': math.sin,
            'cos': math.cos,
            'exp': math.exp,
            'log': math.log,
        }
        for name, fcn in functors.items():
            coeff = 2.34
            func = ct.Func1(name, coeff)
            assert func.cxx_type == f"Cantera::{name.capitalize()}1"
            assert func.type == name
            for val in [.1, 1., 10.]:
                assert name in func.write()
                assert func(val) == pytest.approx(fcn(coeff * val))

    def test_deprecated(self):
        with pytest.warns(DeprecationWarning, match="use alternative constructor"):
            ct.Func1.cxx_functor("cos")

    def test_compound(self):
        functors = {
            'sum': lambda x, y: x + y,
            'diff': lambda x, y: x - y,
            'product': lambda x, y: x * y,
            'ratio': lambda x, y: x / y,
        }
        f1 = ct.Func1('pow', 2)
        f2 = ct.Func1('sin')
        for name, fcn in functors.items():
            func = ct.Func1(name, f1, f2)
            assert func.type == name
            for val in [.1, 1., 10.]:
                assert name not in func.write()
                assert func(val) == pytest.approx(fcn(f1(val), f2(val)))
        f0 = 3.1415
        fcn = lambda x, y: x + y
        func1 = ct.Func1('sum', f0, f1)
        func2 = ct.Func1('sum', f2, f0)
        for val in [.1, 1., 10.]:
            assert func1(val) == pytest.approx(fcn(f0, f1(val)))
            assert func2(val) == pytest.approx(fcn(f2(val), f0))
        with pytest.raises(ValueError):
            ct.Func1('sum', f0, f0)
        with pytest.raises(ValueError):
            ct.Func1('sum', 'spam', 'eggs')

    def test_new_sum(self):
        f1 = ct.Func1('sin', 3.)
        f2 = f1 + 0
        assert f2.type == 'sin'
        assert f2(.5) == np.sin(.5*3)
        f2 = 0 + f1
        assert f2.type == 'sin'
        assert f2(.5) == np.sin(.5*3)
        f2 = f1 + ct.Func1(1.)
        assert f2.type == 'plus-constant'
        assert f2(.5) == np.sin(.5*3) + 1.
        f2 = f1 + f1
        assert f2.type == 'times-constant'
        assert f2(.5) == 2 * np.sin(.5*3)

    def test_new_diff(self):
        f1 = ct.Func1('sin', 3.)
        f2 = f1 - 0
        assert f2(.5) == np.sin(.5*3)
        assert f2.type == 'sin'
        f2 = 0 - f1
        assert f2.type == 'times-constant'
        assert f2(.5) == -np.sin(.5*3)
        f2 = f1 - ct.Func1(1.)
        assert f2.type == 'plus-constant'
        assert f2(.5) == np.sin(.5*3) - 1.
        f2 = f1 - f1
        assert f2.type == 'constant'
        assert f2(.5) == 0.

    def test_new_prod(self):
        f1 = ct.Func1('sin', 3.)
        f2 = f1 * ct.Func1(1.)
        assert f2.type == 'sin'
        assert f2(.5) == np.sin(.5*3)
        f2 = 2 * f1
        assert f2.type == 'times-constant'
        assert f2(.5) == 2 * np.sin(.5*3)
        f2 = f1 * 0
        assert f2.type == 'constant'
        assert f2(.5) == 0.

    def test_new_ratio(self):
        f1 = ct.Func1('sin', 3.)
        f2 = f1 / ct.Func1(1.)
        assert f2.type == 'sin'
        assert f2(.5) == np.sin(.5*3)
        f2 = 0 / f1
        assert f2.type == 'constant'
        assert f2(.5) == 0.
        f2 = f1 / 2
        assert f2.type == 'times-constant'
        assert f2(.5) == .5*np.sin(.5*3)
        f2 = ct.Func1( 0.) / 3.
        assert f2.type == 'constant'
        assert f2(.5) == 0.

    def test_modified(self):
        functors = {
            'plus-constant': lambda x, y: x + y,
            'times-constant': lambda x, y: x * y,
        }
        f1 = ct.Func1('sin')
        constant = 2.34
        for name, fcn in functors.items():
            with pytest.raises(ValueError):
                ct.Func1(name, constant, f1)
            func = ct.Func1(name, f1, constant)
            assert func.type == name
            for val in [.1, 1., 10.]:
                assert name not in func.write()
                assert func(val) == pytest.approx(fcn(f1(val), constant))

    def test_tabulated1(self):
        # this implicitly probes advanced functors
        arr = np.array([[0, 2], [1, 1], [2, 0]])
        time = arr[:, 0]
        fval = arr[:, 1]
        fcn0 = ct.Tabulated1(time, fval)
        fcn1 = ct.Func1("tabulated-linear", time, fval)
        assert fcn0.type == "tabulated-linear"
        assert fcn1.type == "tabulated-linear"
        for t, f in zip(time, fval):
            assert fcn0(t) == pytest.approx(f)
            assert fcn1(t) == pytest.approx(f)

    def test_tabulated2(self):
        time = [0, 1, 2]
        fval = [2, 1, 0]
        fcn = ct.Tabulated1(time, fval)
        assert fcn.type == "tabulated-linear"
        for t, f in zip(time, fval):
            assertNear(f, fcn(t))

    def test_tabulated3(self):
        time = 0, 1, 2,
        fval = 2, 1, 0,
        fcn = ct.Tabulated1(time, fval)
        assertNear(fcn(-1), fval[0])
        assertNear(fcn(3), fval[-1])

    def test_tabulated4(self):
        time = np.array([0, 1, 2])
        fval = np.array([2, 1, 0])
        fcn = ct.Tabulated1(time, fval)
        tt = .5*(time[1:] + time[:-1])
        ff = .5*(fval[1:] + fval[:-1])
        for t, f in zip(tt, ff):
            assertNear(f, fcn(t))

    def test_tabulated5(self):
        time = [0, 1, 2]
        fval = [2, 1, 0]
        fcn = ct.Tabulated1(time, fval, method='previous')
        assert fcn.type == "tabulated-previous"
        val = np.array([fcn(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]])
        assertArrayNear(val, np.array([2.0, 2.0, 2.0, 1.0, 0.0, 0.0]))

    def test_tabulated_failures(self):
        with pytest.raises(ct.CanteraError, match="even number of entries"):
            ct.Tabulated1(range(2), range(3))
        with pytest.raises(ct.CanteraError, match="at least 4 entries"):
            ct.Tabulated1([], [])
        with pytest.raises(ct.CanteraError, match="monotonically"):
            ct.Tabulated1((0, 1, 0.5, 2), (2, 1, 1, 0))
        with pytest.raises(ct.CanteraError, match="No such type"):
            ct.Tabulated1((0, 1, 1, 2), (2, 1, 1, 0), method='spam')
