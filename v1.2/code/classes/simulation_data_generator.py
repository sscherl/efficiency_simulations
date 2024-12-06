#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class to generate random data as a complimentary test to the simulation data created with SciFi geometry by Zurich group
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
import yaml

from classes.simulation_data_plot import Plot
from classes.toolkit import Chip, File, LHCParameters, AstroPixParameters

###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


###### Classes ##################################################

class DataGenerator:

    def __init__(self, events = '', rate = '', chips = ''):

        #self.lhc        = LHCParameters()
        self.params      = AstroPixParameters()
        
        self.chip       = Chip()                            # Default initialisation via yaml config input

        if chips == ''  : self.chips = config['chip']['last'] - config['chip']['first'] + 1
        else            : self.chips = chips

        if events == '' : self.events   = config['events']['last']
        else            : self.events   = events

        if rate == ''   : self.rate     = config['random']['rate']          # Particle hit rate in MHz/cm^2
        else            : self.rate     = rate
    
        self.data_file  = File(path = config['directories']['data']['input'], prefix = 'rnd_data', suffix = '.csv')

        self.calculate_chip2x2_rate()
        #self.generate_data(chip_id = chip_id)
        #self.get_data()


    def calculate_chip2x2_rate(self):
        self.chip2x2_rate = self.rate * 4 * self.params.event_length # cm^2                      # Data simulated for 2 cm x 2 cm area, like scifi data

    def get_number_of_hits(self):
        return np.random.poisson(self.chip2x2_rate)
    
    def get_coordinate(self): # in [-9.6, 9.6] as chip is 19.2 mm wide
        return round(random.uniform(-10, 10), 4)

    def get_tof(self):
        # return round(random.uniform(self.lhc.min_tof_theory, self.lhc.min_tof_theory + self.lhc.bx_period), 4)
        return round(random.uniform(self.params.min_tof_theory, self.params.max_tof), 4)
    
    def get_tot(self):
        return round(random.uniform(1, 2), 4)


    def generate_data(self):                                        # Create dataframe in same format as scifi raw data
        self.df = pd.DataFrame(columns = ['ChipID', 'Event', 'Layer', 'Quadrant', 'X', 'Y', 'ToF', 'Energy'])
        ids, events, xs, ys, tofs, tots = [], [], [], [], [], []

        for id in range(0,self.chips):

            # Generate the data
            for new_event in range(1, self.events+1):

                hits = self.get_number_of_hits()

                for new_hit in range(hits):
                    ids.append(id)
                    events.append(new_event)
                    xs.append(self.get_coordinate())
                    ys.append(self.get_coordinate())
                    tofs.append(self.get_tof())
                    tots.append(self.get_tot())

        empty = [0] * len(events)
        one = [1] * len(events)

        # Fill dataframe
        self.df.ChipID = ids
        self.df.Event = events
        self.df.Layer = one
        self.df.Quadrant = one
        self.df.X = xs
        self.df.Y = ys
        self.df.ToF = empty
        self.df.Energy = tots

    def get_data(self):
        return self.df

    def save_data(self):
        self.df.to_csv(self.data_file.whole, index = False)
