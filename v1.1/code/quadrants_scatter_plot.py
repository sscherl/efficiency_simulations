#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Plot to illustrate how data from different quadrants is overlapped
#
#	Author: Sigrid Scherl
#
#	Created: January 2022
#

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import axes
import matplotlib
import random

from classes.simulation_data_plot import MyColours

###### Main #####################################################

class PlotQuadrants:

    def __init__(self, q1 = True, q2 = True, q3 = True, q4 = True, mirror_q2 = False, mirror_q3 = False, mirror_q4 = False, arrow_q2to1 = False, arrow_q3to1 = False, arrow_q4to1 = False):

        col = MyColours()

        matplotlib.rcParams.update({'font.size': 14})

        # Make randoms stay the same
        random.seed(10)

        self.data_points = 100
        self.data_range = 20

        self.mu = 0
        self.sigma = 8

        self.q1 = q1
        self.q2 = q2
        self.q3 = q3
        self.q4 = q4
        
        self.mirror_q2 = mirror_q2
        self.mirror_q3 = mirror_q3
        self.mirror_q4 = mirror_q4

        self.arrow_q2to1 = arrow_q2to1
        self.arrow_q3to1 = arrow_q3to1
        self.arrow_q4to1 = arrow_q4to1

        self.plot_dir = '/Users/sigridscherl/Documents/PhD/Writings/EfficiencySimulation/Figures/'
        self.plot_name = ''


    def plot_quadrants(self):

        for i in range(self.data_points):
            x = (random.gauss(self.mu, self.sigma))
            y = (random.gauss(self.mu, self.sigma))

            # Q1: (+x,+y)
            if x > 0 and y > 0: clr = 'tab:blue'
            # Q2: (-x,+y)
            elif x < 0 and y > 0: clr = 'tab:orange'
            # Q3 = (-x,-y)
            elif x < 0 and y < 0: clr = 'tab:green'
            # Q4 = (+x,-y)
            elif x > 0 and y < 0: clr = 'tab:red'
            else: print('ERROR: Wonky (x,y)!')

            if not self.q1 and y > 0 and x > 0: continue
            if not self.q2 and y > 0 and x < 0: continue
            if not self.q3 and y < 0 and x < 0: continue
            if not self.q4 and y < 0 and x > 0: continue

            if (self.mirror_q2 and x < 0 and y > 0) or (self.mirror_q3 and x < 0 and y < 0): x = -x
            if (self.mirror_q4 and y < 0 and x > 0) or (self.mirror_q3 and y < 0): y = -y

            plt.scatter(x,y,color=clr, s = 15)
            
            if self.arrow_q2to1: plt.arrow(-5, 5, 10, 0)
            if self.arrow_q3to1: plt.arrow(-5, -5, 10, 10)
            if self.arrow_q4to1: plt.arrow(5, -5, 0, 10)

        x1, y1 = [- self.data_range, self.data_range], [0,0]
        x2, y2 = [0,0], [- self.data_range, self.data_range]

        plt.plot(x1, y1, x2, y2, color='black', linewidth = 1)
        plt.xlabel('x')
        plt.ylabel('y')

        plt.tight_layout()

        plt.savefig(self.plot_dir+self.plot_name+'.pdf')
        plt.show()


plot_all = PlotQuadrants(q1 = True, q2 = True, q3 = True, q4 = True, mirror_q2 = False, mirror_q3 = False, mirror_q4 = False)
plot_all.plot_name = 'QAll'
plot_all.plot_quadrants()

plot_all = PlotQuadrants(q1 = True, q2 = True, q3 = True, q4 = True, mirror_q2 = False, mirror_q3 = False, mirror_q4 = False)
plot_all.plot_name = 'QAll_'
plot_all.plot_quadrants()

plot_q1 = PlotQuadrants(q1 = True, q2 = False, q3 = False, q4 = False, mirror_q2 = False, mirror_q3 = False, mirror_q4 = False)
plot_q1.plot_name = 'Q1'
plot_q1.plot_quadrants()

plot_q2 = PlotQuadrants(q1 = True, q2 = True, q3 = False, q4 = False, mirror_q2 = True, mirror_q3 = False, mirror_q4 = False)
plot_q2.plot_name = 'Q1to2'
plot_q2.plot_quadrants()

plot_q2 = PlotQuadrants(q1 = True, q2 = True, q3 = True, q4 = True, mirror_q2 = True, mirror_q3 = True, mirror_q4 = True)
plot_q2.plot_name = 'Q1to4'
plot_q2.plot_quadrants()

