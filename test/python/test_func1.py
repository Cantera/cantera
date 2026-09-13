import math
import numpy as np
import pytest
from pytest import approx

import cantera as ct


class TestFunc1:
    def test_function(self):
        f = ct.Func1(np.sin)
        assert f(0) == approx(np.sin(0))
        assert f(0.1) == approx(np.sin(0.1))
        assert f(0.7) == approx(np.sin(0.7))

    def test_lambda(self):
        f = ct.Func1(lambda t: np.sin(t)*np.sqrt(t))
        assert f.type == "functor"
        for t in [0.1, 0.7, 4.5]:
            assert f(t) == approx(np.sin(t)*np.sqrt(t))

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
            assert f(t) == approx(8.1*t)

    def test_constant(self):
        f = ct.Func1(5)
        for t in [0.1, 0.7, 4.5]:
            assert f(t) == approx(5)
        assert f.type == "constant"

    def test_sequence(self):
        f = ct.Func1([5])
        for t in [0.1, 0.7, 4.5]:
            assert f(t) == approx(5)

        with pytest.raises(TypeError):
            ct.Func1([3,4])

    def test_numpy(self):
        f = ct.Func1(np.array(5))
        assert f.type == "constant"
        g = ct.Func1(np.array([[5]]))
        assert g.type == "constant"
        for t in [0.1, 0.7, 4.5]:
            assert f(t) == approx(5)
            assert g(t) == approx(5)

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
                assert func(val) == approx(fcn(coeff * val))

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
                assert func(val) == approx(fcn(f1(val), f2(val)))
        f0 = 3.1415
        fcn = lambda x, y: x + y
        func1 = ct.Func1('sum', f0, f1)
        func2 = ct.Func1('sum', f2, f0)
        for val in [.1, 1., 10.]:
            assert func1(val) == approx(fcn(f0, f1(val)))
            assert func2(val) == approx(fcn(f2(val), f0))
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
                assert func(val) == approx(fcn(f1(val), constant))

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
            assert fcn0(t) == approx(f)
            assert fcn1(t) == approx(f)

    def test_tabulated2(self):
        time = [0, 1, 2]
        fval = [2, 1, 0]
        fcn = ct.Tabulated1(time, fval)
        assert fcn.type == "tabulated-linear"
        for t, f in zip(time, fval):
            assert f == approx(fcn(t))

    def test_tabulated3(self):
        time = 0, 1, 2,
        fval = 2, 1, 0,
        fcn = ct.Tabulated1(time, fval)
        assert fcn(-1) == approx(fval[0])
        assert fcn(3) == approx(fval[-1])

    def test_tabulated4(self):
        time = np.array([0, 1, 2])
        fval = np.array([2, 1, 0])
        fcn = ct.Tabulated1(time, fval)
        tt = .5*(time[1:] + time[:-1])
        ff = .5*(fval[1:] + fval[:-1])
        for t, f in zip(tt, ff):
            assert f == approx(fcn(t))

    def test_tabulated5(self):
        time = [0, 1, 2]
        fval = [2, 1, 0]
        fcn = ct.Tabulated1(time, fval, method='previous')
        assert fcn.type == "tabulated-previous"
        val = np.array([fcn(v) for v in [-0.5, 0, 0.5, 1.5, 2, 2.5]])
        assert val == approx(np.array([2.0, 2.0, 2.0, 1.0, 0.0, 0.0]))

    def test_tabulated_failures(self):
        with pytest.raises(ct.CanteraError, match="even number of entries"):
            ct.Tabulated1(range(2), range(3))
        with pytest.raises(ct.CanteraError, match="at least 4 entries"):
            ct.Tabulated1([], [])
        with pytest.raises(ct.CanteraError, match="monotonically"):
            ct.Tabulated1((0, 1, 0.5, 2), (2, 1, 1, 0))
        with pytest.raises(ct.CanteraError, match="No such type"):
            ct.Tabulated1((0, 1, 1, 2), (2, 1, 1, 0), method='spam')


class TestFunc1UseAfterFree:
    """Regression tests for GitHub issue #2161.

    Python-callable Func1 objects (lambdas, functions) must remain alive as
    long as any compound or modified functor that depends on them is alive.
    Previously, the composed functor did not hold a Python reference to its
    base(s), so the garbage collector could free the base while the C++ layer
    still held a raw pointer into its memory, causing a segmentation fault.
    """

    def _lambda_func(self, t):
        return 2.0 * t

    def test_mul_scalar_keeps_lambda_alive(self):
        """f = lambda_func * scalar should survive gc.collect()."""
        import gc
        f = ct.Func1(lambda t: 2.0 * t) * 0.5
        gc.collect()
        assert f(4.0) == approx(4.0)

    def test_rmul_scalar_keeps_lambda_alive(self):
        """f = scalar * lambda_func should survive gc.collect()."""
        import gc
        f = 0.5 * ct.Func1(lambda t: 2.0 * t)
        gc.collect()
        assert f(4.0) == approx(4.0)

    def test_add_scalar_keeps_lambda_alive(self):
        """f = lambda_func + scalar should survive gc.collect()."""
        import gc
        f = ct.Func1(lambda t: t) + 1.0
        gc.collect()
        assert f(3.0) == approx(4.0)

    def test_radd_scalar_keeps_lambda_alive(self):
        """f = scalar + lambda_func should survive gc.collect()."""
        import gc
        f = 1.0 + ct.Func1(lambda t: t)
        gc.collect()
        assert f(3.0) == approx(4.0)

    def test_sub_scalar_keeps_lambda_alive(self):
        """f = lambda_func - scalar should survive gc.collect()."""
        import gc
        f = ct.Func1(lambda t: t) - 1.0
        gc.collect()
        assert f(5.0) == approx(4.0)

    def test_rsub_scalar_keeps_lambda_alive(self):
        """f = scalar - lambda_func should survive gc.collect()."""
        import gc
        f = 10.0 - ct.Func1(lambda t: t)
        gc.collect()
        assert f(3.0) == approx(7.0)

    def test_truediv_scalar_keeps_lambda_alive(self):
        """f = lambda_func / scalar should survive gc.collect()."""
        import gc
        f = ct.Func1(lambda t: 4.0 * t) / 2.0
        gc.collect()
        assert f(3.0) == approx(6.0)

    def test_rtruediv_scalar_keeps_lambda_alive(self):
        """f = scalar / lambda_func should survive gc.collect()."""
        import gc
        f = 12.0 / ct.Func1(lambda t: t)
        gc.collect()
        assert f(3.0) == approx(4.0)

    def test_string_compound_keeps_lambda_alive(self):
        """Func1('sum', lambda_func, other) should survive gc.collect()."""
        import gc
        f = ct.Func1("sum", ct.Func1(lambda t: t), ct.Func1(lambda t: 2.0 * t))
        gc.collect()
        assert f(3.0) == approx(9.0)

    def test_string_product_keeps_lambda_alive(self):
        """Func1('product', lambda_func, other) should survive gc.collect()."""
        import gc
        f = ct.Func1("product", ct.Func1(lambda t: t), ct.Func1(lambda t: 2.0))
        gc.collect()
        assert f(3.0) == approx(6.0)

    def test_string_composite_keeps_lambda_alive(self):
        """Func1('composite', f, g) computes f(g(t)); both survive gc."""
        import gc
        f = ct.Func1("composite", ct.Func1(lambda t: 2.0 * t),
                     ct.Func1(lambda t: t + 1.0))
        gc.collect()
        # composite(f, g)(t) = f(g(t)) = 2*(t+1)
        assert f(3.0) == approx(8.0)

    def test_string_modified_times_constant_keeps_lambda_alive(self):
        """Func1('times-constant', lambda_func, c) should survive gc.collect()."""
        import gc
        f = ct.Func1("times-constant", ct.Func1(lambda t: t), 3.0)
        gc.collect()
        assert f(4.0) == approx(12.0)

    def test_string_modified_periodic_keeps_lambda_alive(self):
        """Func1('periodic', lambda_func, period) should survive gc.collect()."""
        import gc
        f = ct.Func1("periodic", ct.Func1(lambda t: t), 5.0)
        gc.collect()
        # periodic wraps t into [0, period); t=3 is already in [0,5)
        assert f(3.0) == approx(3.0)

    def test_deep_chain_keeps_all_lambdas_alive(self):
        """(f * g) + h where all three are lambdas; all survive gc."""
        import gc
        product = ct.Func1(lambda t: t) * ct.Func1(lambda t: 2.0)
        f = product + ct.Func1(lambda t: 1.0)
        del product
        gc.collect()
        # f(t) = t*2 + 1
        assert f(3.0) == approx(7.0)

    def test_callable_class_survives_gc(self):
        """Func1 wrapping a callable class instance survives gc."""
        import gc

        class Scale:
            def __init__(self, factor):
                self.factor = factor
            def __call__(self, t):
                return self.factor * t

        f = ct.Func1(Scale(3.0)) * 2.0
        gc.collect()
        assert f(5.0) == approx(30.0)

    def test_method_survives_gc(self):
        """Func1 wrapping a bound method survives gc."""
        import gc
        f = ct.Func1(self._lambda_func) * 0.5
        gc.collect()
        assert f(4.0) == approx(4.0)

    def test_reactor_net_survives_gc(self):
        """Regression for the original crash from issue #2161.

        MassFlowController mdot set to (base_func * scalar); base_func goes
        out of scope before ReactorNet.advance() is called.
        """
        import gc

        def build():
            gas = ct.Solution('gri30.yaml')
            gas.TPX = 300, 5e5, 'CH4:1'
            r1 = ct.IdealGasReactor(gas, volume=0.005)
            r2 = ct.IdealGasReactor(gas, volume=0.005)
            fuel = ct.Reservoir(ct.Solution('gri30.yaml'))
            ox = ct.Reservoir(ct.Solution('gri30.yaml'))
            env = ct.Reservoir(ct.Solution('air.yaml'))

            # base Func1 objects deliberately not stored after build() returns
            base_fuel = ct.Func1(lambda t: 1.0)
            base_ox = ct.Func1(lambda t: 40.0)
            ct.MassFlowController(fuel, r1, mdot=base_fuel * 0.7)
            ct.MassFlowController(fuel, r2, mdot=base_fuel * 0.3)
            ct.MassFlowController(ox, r1, mdot=base_ox * 0.3)
            ct.MassFlowController(ox, r2, mdot=base_ox * 0.7)
            ct.Valve(r1, r2, K=1.0)
            ct.MassFlowController(r2, env, mdot=ct.Func1(lambda t: 1.0))
            return ct.ReactorNet([r1, r2])

        nets = [build() for _ in range(4)]
        gc.collect()  # would previously free base_fuel / base_ox
        for net in nets:
            net.advance(0)
            net.advance(1e-4)  # previously crashed with SIGSEGV here

