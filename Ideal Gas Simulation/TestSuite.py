"""Test suite for ideal gas simulation project.
"""
import unittest, unittest.mock

from tests.test_gui import TestIdealGasGUI
from tests.test_idealgas import TestIdealGas, TestParticle

def suite():
    suite = unittest.TestSuite()
    match 2:
        case 0:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestIdealGasGUI))
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestIdealGas))
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestParticle))
        case 1:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestIdealGasGUI))
        case 2:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestIdealGas))
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestParticle))
    return suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite())