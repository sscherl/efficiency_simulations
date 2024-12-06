#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Generate map of hit rates with Oscar's fit
#
#	Author: Sigrid Scherl
#
#	Created: March 2023
#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import csv

class MTHitRates:

    def __init__(self):
        self.tracker_width  = 5000 #mm
        self.tracker_height = 1200 #mm
        self.sensor_width   = 20 #mm
        self.sensor_height  = 20 #mm

        self.q1_p0  = -0.0042579
        self.q1_p1  = -2.63322e-07
        self.q1_p2  = -0.0109647
        self.q1_p3  = 1.79803e-06
        self.q1_p4  = 1.1804e-05
        self.q1_p5  = 5.84583

    def hitratefit(self, x, y):
        p0 = self.q1_p0
        p1 = self.q1_p1
        p2 = self.q1_p2
        p3 = self.q1_p3
        p4 = self.q1_p4
        p5 = self.q1_p5

        z = np.exp(p0 * x + p1 * x**2 + p2 * y + p3 * y**2 + p4 * x * y + p5)

        return z
    

    def hitratefit_donal(self, x, y):
        p0 = self.q1_p0
        p1 = self.q1_p1
        p2 = self.q1_p2
        p3 = self.q1_p3
        p4 = self.q1_p4
        p5 = self.q1_p5

        z = p0 * x + p1 * x**2 + p2 * y + p3 * y**2 + p4 * x * y + p5

        return z
    

hitrates = MTHitRates()

print(hitrates.hitratefit(154, 12))