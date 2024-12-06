#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class to generate Monte Carlo simulation data as a complimentary test to the simulation data created with SciFi geometry by Zurich group
#
#   Note: Rate = Hits per event (25 ns) and area (chip or pixel area in mm^2)
#
#	Author: Sigrid Scherl
#
#	Created: March 2022
#

import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import random
from random import randint
import sys

from classes.simulation_data_plot import Plot
from classes.toolkit import Chip, File, LHCParameters

###### Classes ##################################################

class DataGenerator:

    def __init__(self, chip, events, rate):
        self.chip = chip
        self.events = events
        self.full_chip_rate = rate
        self.full_data_rate = rate / (20*19.2) * (20*20)   # Hits per 20 mm x 19.2 mm chip and event (25 ns)

        self.full_data_length = 20 #mm
        self.full_data_width = 20 #mm
        self.full_data_size = self.full_data_length * self.full_data_width

        self.lhc = LHCParameters()

        self.data_file = File(path = '../data/rnd_input_data/', suffix = '.csv')
        self.data_file.generate_name(prefix = 'rnd_data', chip = self.chip, first_event = 1, last_event = self.events, rate = self.full_data_rate)

        self.set_parameters()

    def set_parameters(self):
        # Rate per pixel (Divide rate per hit by pixels)
        self.chip_rate = self.full_data_rate / self.full_data_size * self.chip.size
        self.pixel_rate = self.full_data_rate / self.chip.pixels

    def nr_of_hits(self):
        return np.random.poisson(self.full_data_rate)

    def get_x_coord(self): # in [-9.6, 9.6] as chip is 19.2 mm wide
        return round(random.uniform(- self.full_data_width/2, self.full_data_width/2), 4)

    def get_y_coord(self): # in [-10.0, 10.0] as chip is 20 mm long
        return round(random.uniform(- self.full_data_length/2, self.full_data_length/2), 4)

    def get_tof(self):
        # return round(random.uniform(self.lhc.min_tof_theory, self.lhc.min_tof_theory + self.lhc.bx_period), 4)
        return round(random.uniform(self.lhc.min_tof_theory, self.lhc.max_tof), 4)

    def get_energy(self):
        return round(random.uniform(0.1, 2.0), 4)

    def get_col(self):
        return randint(0, self.chip.cols)

    def get_row(self):
        return randint(0, self.chip.rows)

    def generate_data(self):
        # Create dataframe in same format as scifi raw data
        self.df = pd.DataFrame(columns = ['Event', 'Layer', 'Quadrant', 'X', 'Y', 'ToF', 'Energy'])
        events, xs, ys, tofs, energies = [], [], [], [], []

        # Generate the data
        for new_event in range(1, self.events+1):
            for new_hit in range(self.nr_of_hits()):
                events.append(new_event)
                xs.append(self.get_x_coord())
                ys.append(self.get_y_coord())
                tofs.append(self.get_tof())
                energies.append(self.get_energy())

        empty = [0] * len(events)
        one = [1] * len(events)

        # Fill dataframe
        self.df.Event = events
        self.df.Layer = one
        self.df.Quadrant = one
        self.df.X = xs
        self.df.Y = ys
        self.df.ToF = tofs
        self.df.Energy = energies

    def get_data(self):
        return self.df

    def save_data(self):
        self.df.to_csv(self.data_file.whole, index = False)
