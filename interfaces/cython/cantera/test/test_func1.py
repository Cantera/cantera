import numpy as np

import cantera as ct
from . import utilities
from .utilities import unittest

class TestFunc1(utilities.CanteraTest):
    def test_function(self):
        f = ct.Func1(np.sin)
        self.assertNear(f(0), np.sin(0))
        self.assertNear(f(0.1), np.sin(0.1))
        self.assertNear(f(0.7), np.sin(0.7))

    def test_lambda(self):
        f = ct.Func1(lambda t: np.sin(t)*np.sqrt(t))
        for t in [0.1, 0.7, 4.5]:
            self.assertNear(f(t), np.sin(t)*np.sqrt(t))

    def test_callable(self):
        class Multiplier:
            def __init__(self, factor):
                self.factor = factor
            def __call__(self, t):
                return self.factor * t

        m = Multiplier(8.1)
        f = ct.Func1(m)
        for t in [0.1, 0.7, 4.5]:
            self.assertNear(f(t), 8.1*t)

    def test_constant(self):
        f = ct.Func1(5)
        for t in [0.1, 0.7, 4.5]:
            self.assertNear(f(t), 5)

    def test_sequence(self):
        f = ct.Func1([5])
        for t in [0.1, 0.7, 4.5]:
            self.assertNear(f(t), 5)

        with self.assertRaises(TypeError):
            ct.Func1([3,4])

    def test_numpy(self):
        f = ct.Func1(np.array(5))
        g = ct.Func1(np.array([[5]]))
        for t in [0.1, 0.7, 4.5]:
            self.assertNear(f(t), 5)
            self.assertNear(g(t), 5)

        with self.assertRaises(TypeError):
            ct.Func1(np.array([3,4]))

    def test_failure(self):
        def fails(t):
            raise ValueError('bad')

        f = ct.Func1(fails)
        with self.assertRaises(ValueError):
            f(0.1)

    def test_unpicklable(self):
        import pickle
        f = ct.Func1(np.sin)
        with self.assertRaises(NotImplementedError):
            pickle.dumps(f)

    def test_uncopyable(self):
        import copy
        f = ct.Func1(np.sin)
        with self.assertRaises(NotImplementedError):
            copy.copy(f)
