import sys
from unittest import TestCase

from PyQt6 import QtWidgets

from gui import IdealGasGUI

class TestIdealGasGUI(TestCase):
    
    @classmethod
    def setUpClass(self) -> None:
        self.qapp = QtWidgets.QApplication(sys.argv)
        
    @classmethod
    def tearDownClass(self) -> None:
        self.qapp.quit()
        
    def setUp(self) -> None:
        self.app = IdealGasGUI(show=False)
        
    def tearDown(self) -> None:
        self.app.close()
    
    def test_init(self):
        raise NotImplementedError("Test not implemented.")