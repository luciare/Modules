# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 10:52:38 2020

@author: Lucia
"""

import numpy as np
import time
from PyQt5 import Qt


class GeneralTimer(Qt.QThread):
    TimerDone = Qt.pyqtSignal()

    def __init__(self):
        super(GeneralTimer, self).__init__()
        self.InitTime = time.time()
        self.ElapsedTime = 0.00
        
    def run(self):
        while True:
            NewTime = time.time()
            self.ElapsedTime = NewTime - self.InitTime
            self.TimerDone.emit()
            Qt.QThread.msleep(1000)

