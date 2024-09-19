"""Test suite for ideal gas simulation project."""

import os
import unittest, unittest.mock

import sys
sys.path.append(os.path.dirname(r"")) #TODO: Insert parent directory here.

from test_gui import TestGUI
from test_idealgas import TestIdealGas

def suite():
    suite = unittest.TestSuite()
    match 0:
        case 0:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestGUI))
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestIdealGas))
        case 1:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestGUI))
        case 2:
            suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(TestIdealGas))
    return suite

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=1)
    runner.run(suite())