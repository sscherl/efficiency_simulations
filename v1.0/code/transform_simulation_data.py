#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Output simulation data for chosen input data in proper format for pixel matrix model
#
#	Author: Sigrid Scherl
#
#	Created: November 2021
#

import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

from simulation_data_transformation import DataTransformation
from toolkit import Chip, Helper

###### Run options ##############################################

# Mute SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'

###### Settings #################################################

helper = Helper()
mightypix = Chip()
# mightypix = Chip(name = 'MightyPix', cols = 128, rows = 400, col_width = 0.150, row_width = 0.05)

data_types = ['scifi', 'random']
quadrants = [1,2,3,4] # Single to quadruple data rate for simulation data
rates = [1.7, 3.4, 5.1, 6.8] # Single to quadruple data rate for random data
secs = [True]
clstrs = [False]
tof_prctgs = [100, 99]
events = [3000]

# Run simulation data:
for quad in quadrants:
    for sec in secs:
        for clstr in clstrs:
            for tof_prctg in tof_prctgs:

                trafo = DataTransformation('scifi', mightypix, events = 500, layers = 6, quadrants = quad)
                trafo.add_output_path('../data/input_data/')
                # trafo.add_output_path('../../../mightypix/verification/simulation_data/')

                # Conversion factor of energy in MeV to ToT in us, if set to 0 fixed tot value will be used
                trafo.mev_to_us = 0
                trafo.fixed_tot = 2 #us

                # Run options
                trafo.keep_secondaries = sec
                trafo.clusters = clstr

                trafo.tof_plot_percentage = tof_prctg

                ###### Run ######################################################

                trafo.transform(save_me = True)
                trafo.plots(save_me = True, show_me = False)

# Run random data:
for rate in rates:
    for clstr in clstrs:
        for event in events:

            trafo = DataTransformation('random', mightypix, events = event, layers = 1, quadrants = 1, rate = rate)
            trafo.add_output_path('../data/input_data/')
            # trafo.add_output_path('../../../mightypix/verification/simulation_data/')

            # Conversion factor of energy in MeV to ToT in us, if set to 0 fixed tot value will be used
            trafo.mev_to_us = 0
            trafo.fixed_tot = 2 #us

            # Run options
            trafo.keep_secondaries = True
            trafo.clusters = clstr

            trafo.tof_plot_percentage = 100

            ###### Run ######################################################

            trafo.transform(save_me = True)
            trafo.plots(save_me = True, show_me = False)
