import unittest

import Cantera as ct

class ImportTest(unittest.TestCase):
    def testCtiImport(self):
        gas = ct.importPhase('test/data/air-no-reactions.cti', 'air')

        self.assertAlmostEqual(gas.temperature(), 300)
