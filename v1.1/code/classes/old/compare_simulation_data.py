#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Compare the data going into and coming out of MightyPix pixel matrix model
#
#	Author: Sigrid Scherl
#
#	Created: January 2022
#

import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import yaml
from tqdm import tqdm

from classes.simulation_data_comparison import DataComparison
from classes.toolkit import Chip, Helper, File, ProgressBar

###### Run options ##############################################

# Mute SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'

###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Run ######################################################

def generate_list(name):
    if config[name]['list'] != []:  list = config[name]['list']
    else:                           list = [*range(config[name]['first'], config[name]['last']+1)]
    return list


sensors     = [739]       #generate_list('sensors')
layers      = generate_list('layers')
quadrants   = generate_list('quadrants')
events      = 500                            #generate_list('events')
rates       = ['']                #generate_list('rates')



for l in tqdm(layers, desc = 'Layers', leave = False):
    for q in tqdm(quadrants, desc = 'Quadrants', leave = False):
        for s in tqdm(sensors, desc = 'Sensors', leave = False):
            for r in tqdm(rates, desc = 'Rates', leave = False):
                #try:
                comp = DataComparison(sensor = s, events = events, layers = l, quadrants = q, rate = r)
                comp.run()
                #except: print('Skipping...')

print('Done!')
