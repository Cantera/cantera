from os.path import join as pjoin
import os
# from pathlib import Path
# import numpy as np
# import gc

import numpy as np
import warnings

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

import cantera as ct
from . import utilities


class TestModels(utilities.CanteraTest):

    def test_load_thermo_models(self):
        yml_file = pjoin(self.test_data_dir, "thermo-models.yaml")
        with open(yml_file, 'rt', encoding="utf-8") as stream:
            yml = yaml.safe_load(stream)

        for ph in yml['phases']:

            ph_name = ph['name']
            try:
                sol = ct.Solution(yml_file, ph_name)

                T0, p0 = sol.TP
                z = sol.state # calls Phase::saveState
                sol.TP = 300, 2*ct.one_atm
                sol.state = z # calls Phase::restoreState
                self.assertEqual(sol.T, T0)
                self.assertEqual(sol.P, p0)

                if sol._default_state == 2:
                    # stoich phase (fixed composition)
                    self.assertEqual(sol.n_species, 1)
                    self.assertEqual(len(z), 2)
                else:
                    self.assertEqual(len(z), 2 + sol.n_species)

            except Exception as inst:

                # raise meaningful error message without breaking test suite
                # ignore deprecation warnings originating in C++ layer
                # (converted to errors in test suite)
                if 'Deprecated' not in str(inst):

                    msg = "Error in processing of phase with type '{}'"
                    raise TypeError(msg.format(ph['thermo'])) from inst

    def test_restore_thermo_models(self):

        def check(a, b):
            self.assertArrayNear(a.T, b.T)
            self.assertArrayNear(a.P, b.P)
            self.assertArrayNear(a.X, b.X)

        yml_file = pjoin(self.test_data_dir, "thermo-models.yaml")
        with open(yml_file, 'rt', encoding="utf-8") as stream:
            yml = yaml.safe_load(stream)

        for ph in yml['phases']:

            skipped = ['ions-from-neutral-molecule',
                       'pure-fluid']
            if ph['thermo'] in skipped:
                continue

            ph_name = ph['name']

            try:
                sol = ct.Solution(yml_file, ph_name)

                a = ct.SolutionArray(sol, 10)

                # assign some state
                T = 373.15 + 100*np.random.rand(10)
                if sol.n_species > 1:
                    X = a.X
                    X[:, 1] = .01
                    X = np.diag(X.sum(axis=1)).dot(X)
                    a.TPX = T, np.linspace(1., 2., 10), X
                else:
                    a.TP = T, np.linspace(1., 2., 10)

                # default columns
                data, labels = a.collect_data()
                b = ct.SolutionArray(sol)
                b.restore_data(data, labels)
                check(a, b)

            except Exception as inst:

                # raise meaningful error message without breaking test suite
                # ignore deprecation warnings originating in C++ layer
                # (converted to errors in test suite)
                if 'Deprecated' not in str(inst):

                    msg = "Error in processing of phase with type '{}'"
                    raise TypeError(msg.format(ph['thermo'])) from inst


class TestRestoreIdealGas(utilities.CanteraTest):
    """ Test restoring of the IdealGas class """
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.gas = ct.Solution('h2o2.xml')

    def test_restore_gas(self):

        def check(a, b, atol=None):
            if atol is None:
                self.assertArrayNear(a.T, b.T)
                self.assertArrayNear(a.P, b.P)
                self.assertArrayNear(a.X, b.X)
            else:
                self.assertArrayNear(a.T, b.T, atol=atol)
                self.assertArrayNear(a.P, b.P, atol=atol)
                self.assertArrayNear(a.X, b.X, atol=atol)

        # test ThermoPhase
        a = ct.SolutionArray(self.gas)
        for i in range(10):
            T = 300 + 1800*np.random.random()
            P = ct.one_atm*(1 + 10*np.random.random())
            X = np.random.random(self.gas.n_species)
            X[-1] = 0.
            X /= X.sum()
            a.append(T=T, P=P, X=X)

        data, labels = a.collect_data()

        # basic restore
        b = ct.SolutionArray(self.gas)
        b.restore_data(data, labels)
        check(a, b)

        # skip concentrations
        b = ct.SolutionArray(self.gas)
        b.restore_data(data[:, :2], labels[:2])
        self.assertTrue(np.allclose(a.T, b.T))
        self.assertTrue(np.allclose(a.density, b.density))
        self.assertFalse(np.allclose(a.X, b.X))

        # wrong data shape
        b = ct.SolutionArray(self.gas)
        with self.assertRaises(TypeError):
            b.restore_data(data.ravel(), labels)

        # inconsistent data
        b = ct.SolutionArray(self.gas)
        with self.assertRaises(ValueError):
            b.restore_data(data, labels[:-2])

        # inconsistent shape of receiving SolutionArray
        b = ct.SolutionArray(self.gas, 9)
        with self.assertRaises(ValueError):
            b.restore_data(data, labels)

        # incomplete state
        b = ct.SolutionArray(self.gas)
        with self.assertRaises(ValueError):
            b.restore_data(data[:,1:], labels[1:])

        # add extra column
        t = np.arange(10, dtype=float)[:, np.newaxis]

        # auto-detection of extra
        b = ct.SolutionArray(self.gas)
        b.restore_data(np.hstack([t, data]), ['time'] + labels)
        check(a, b)

        # explicit extra
        b = ct.SolutionArray(self.gas, extra=('time',))
        b.restore_data(np.hstack([t, data]), ['time'] + labels)
        check(a, b)
        self.assertTrue((b.time == t.ravel()).all())

        # wrong extra
        b = ct.SolutionArray(self.gas, extra=('xyz',))
        with self.assertRaises(KeyError):
            b.restore_data(np.hstack([t, data]), ['time'] + labels)

        # missing extra
        b = ct.SolutionArray(self.gas, extra=('time'))
        with self.assertRaises(KeyError):
            b.restore_data(data, labels)

        # inconsistent species
        labels[-1] = 'Y_invalid'
        b = ct.SolutionArray(self.gas)
        with self.assertRaises(ValueError):
            b.restore_data(data, labels)

        # incomplete species info (using threshold)
        data, labels = a.collect_data(threshold=1e-6)

        # basic restore
        b = ct.SolutionArray(self.gas)
        b.restore_data(data, labels)
        check(a, b, atol=1e-6)

        # skip calculated properties
        cols = ('T', 'P', 'X', 'gibbs_mass', 'forward_rates_of_progress')
        data, labels = a.collect_data(cols=cols, threshold=1e-6)

        b = ct.SolutionArray(self.gas)
        b.restore_data(data, labels)
        check(a, b)
        self.assertTrue(len(b._extra) == 0)


class TestRestorePureFluid(utilities.CanteraTest):
    """ Test restoring of the PureFluid class """
    @classmethod
    def setUpClass(cls):
        utilities.CanteraTest.setUpClass()
        cls.water = ct.Water()

    def test_restore_water(self):

        def check(a, b):
            self.assertArrayNear(a.T, b.T)
            self.assertArrayNear(a.P, b.P)
            self.assertArrayNear(a.X, b.X)

        # benchmark
        a = ct.SolutionArray(self.water, 10)
        a.TX = 373.15, np.linspace(0., 1., 10)

        # complete data
        cols = ('T', 'P', 'X')
        data, labels = a.collect_data(cols=cols)
        b = ct.SolutionArray(self.water)
        b.restore_data(data, labels)
        check(a, b)

        # partial data
        cols = ('T', 'X')
        data, labels = a.collect_data(cols=cols)
        b = ct.SolutionArray(self.water)
        b.restore_data(data, labels)
        check(a, b)

        # default columns
        data, labels = a.collect_data()
        self.assertEqual(labels, ['T', 'density'])
        b = ct.SolutionArray(self.water)
        b.restore_data(data, labels)
        check(a, b)

        # default state plus Y
        cols = ('T', 'density', 'Y')
        data, labels = a.collect_data(cols=cols)
        b = ct.SolutionArray(self.water)
        b.restore_data(data, labels)
        check(a, b)
