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
import yaml
from tqdm import tqdm
from multiprocessing import Pool

from classes.simulation_data_transformation import DataTransformation
from classes.toolkit import Chip, Helper

###### Run options ##############################################

pd.options.mode.chained_assignment = None  # Mute SettingWithCopyWarning, default='warn'

###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Functions ################################################

def run_trafo(sensors, events, layers, quadrants, rates, save_trafo = True, save_plots = False, show_plots = False):

    for l in tqdm(layers, desc = 'Layers', leave = False):
        for q in tqdm(quadrants, desc = 'Quadrants', leave = False):
            for s in tqdm(sensors, desc = 'Sensors', leave = False):
                for e in tqdm(events, desc = 'Events', leave = False):
                    for r in tqdm(rates, desc = 'Rates', leave = False):

                        trafo = DataTransformation(s, e, l, q, r)
                        trafo.add_output_path(config['directories']['data']['input'])
                        trafo.transform(save_me = save_trafo)
                        if save_plots or show_plots: trafo.plots(save_me = save_plots, show_me = show_plots)

def generate_list(name):
    if config[name]['list'] != []:  list = config[name]['list']
    else:                           list = [*range(config[name]['first'], config[name]['last']+1)]
    return list


###### Run ######################################################

sensors     = [739]                             #generate_list('sensors')
layers      = generate_list('layers')
quadrants   = generate_list('quadrants')
events      = [500]                             #generate_list('events')
rates       = ['']                                #generate_list('rates')



run_trafo(sensors, events, layers, quadrants, rates)