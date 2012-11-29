import unittest
import numpy as np

import cantera as ct
from . import utilities

class TestTransport(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.Solution('h2o2.xml')
        self.phase.X = [0.1, 1e-4, 1e-5, 0.2, 2e-4, 0.3, 1e-6, 5e-5, 0.4]
        self.phase.TP = 800, 2*ct.OneAtm

    def test_scalar_properties(self):
        self.assertTrue(self.phase.viscosity > 0.0)
        self.assertTrue(self.phase.thermalConductivity > 0.0)

    def test_mixtureAveraged(self):
        self.assertEqual(self.phase.transportModel, 'Mix')
        Dkm1 = self.phase.mixDiffCoeffs
        Dkm1b = self.phase.mixDiffCoeffsMole
        Dkm1c = self.phase.mixDiffCoeffsMass
        Dbin1 = self.phase.binaryDiffCoeffs

        self.phase.transportModel = 'Multi'
        Dkm2 = self.phase.mixDiffCoeffs
        Dkm2b = self.phase.mixDiffCoeffsMole
        Dkm2c = self.phase.mixDiffCoeffsMass
        Dbin2 = self.phase.binaryDiffCoeffs
        self.assertArrayNear(Dkm1, Dkm2)
        self.assertArrayNear(Dkm1b, Dkm2b)
        self.assertArrayNear(Dkm1c, Dkm2c)
        self.assertArrayNear(Dbin1, Dbin2)
        self.assertArrayNear(Dbin1, Dbin1.T)

    def test_multiComponent(self):
        self.assertRaises(Exception,
                          lambda: self.phase.MultiDiffCoeffs)

        self.assertArrayNear(self.phase.thermalDiffCoeffs,
                             np.zeros(self.phase.nSpecies))

        self.phase.transportModel = 'Multi'
        self.assertTrue(all(self.phase.multiDiffCoeffs.flat >= 0.0))
        self.assertTrue(all(self.phase.thermalDiffCoeffs.flat != 0.0))


class TestDustyGas(utilities.CanteraTest):
    def setUp(self):
        self.phase = ct.DustyGas('h2o2.xml')
        self.phase.porosity = 0.2
        self.phase.tortuosity = 0.3
        self.phase.meanPoreRadius = 1e-4
        self.phase.meanParticleDiameter = 5e-4
        self.Dref = self.phase.multiDiffCoeffs

    def test_porosity(self):
        self.phase.porosity = 0.4
        D = self.phase.multiDiffCoeffs
        self.assertArrayNear(self.Dref * 2, D)

    def test_tortuosity(self):
        self.phase.tortuosity = 0.6
        D = self.phase.multiDiffCoeffs
        self.assertArrayNear(self.Dref * 0.5, D)

    # The other parameters don't have such simple relationships to the diffusion
    # coefficients, so we can't test them as easily
