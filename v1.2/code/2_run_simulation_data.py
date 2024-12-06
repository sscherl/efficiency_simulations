#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Run transformation or comparison of simulation data - generate the input and analyse the output files
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

from classes.simulation_data_transformation import DataTransformation
from classes.simulation_data_comparison import DataComparison
from classes.toolkit import Helper

###### Run options ##############################################

# Mute SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'

###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Run ######################################################

help = Helper()

sim         = config['sim']
data        = config['data']['type']

if config['sensors']['use_list']:
    sensors = config['sensors']['list']
else:
    sensors = help.generate_list(config['sensors']['first'],config['sensors']['last'])

if data == 'random':

    events = config['events']['list']

    if config['rates']['use_list']:
        rates = config['rates']['list']
    else:
        rates = help.generate_list(config['rates']['first'],config['rates']['last'])

    layers = [1]
    quadrants = [1]

elif data == 'scifi':

    events = [500]
    rates = ['']

    if config['layers']['use_list']:
        layers = config['layers']['list']
    else:
        layers = help.generate_list(config['layers']['first'],config['layers']['last'])

    if config['quadrants']['use_list']:
        quadrants = config['quadrants']['list']
    else:
        quadrants = help.generate_list(config['quadrants']['first'],config['quadrants']['last'])

else: print('Invalid data type!')

repetitions = config['repetitions']

tots        = config['tot']['fixed']
secondaries = config['secondaries']
clusters    = config['clusters']['generate']
extra       = config['extra']

save_plots = config['plots']['save']

dis = not config['progress_bars']



for s in tqdm(sensors, desc = 'Sensors', leave = False, disable=dis):
    for e in tqdm(events, desc = 'Events', leave = False, disable=dis):
        for r in tqdm(rates, desc = 'Rates', leave = False, disable=dis):
            for l in tqdm(layers, desc = 'Layers', leave = False, disable=dis):
                for q in tqdm(quadrants, desc = 'Quadrants', leave = False, disable=dis):
                    for secs in tqdm(secondaries, desc = 'Secondaries', leave = False, disable=dis):
                        for cl in tqdm(clusters, desc = 'Clusters', leave = False, disable=dis):
                            for tot in tqdm(tots, desc = 'ToT', leave = False, disable=dis):
                                for ex in tqdm(extra, desc = 'Extra', leave = False, disable=dis):
                                    for rep in tqdm(repetitions, desc = 'Reps', leave = False, disable=dis):
                                    
                                        if sim == 'trafo':

                                            trafo = DataTransformation(data = data, sensor = s, events = e, rate = r, layers = l, quadrants = q, secondaries = secs, clusters = cl, tot = tot, extra = ex, repetition = rep, save_plots = save_plots)
                                            trafo.add_output_path(config['directories']['data']['input'])
                                            trafo.transform(save_me = True)

                                        elif sim == 'comp':

                                            # try:
                                            comp = DataComparison(data = data, sensor = s, events = e, rate = r, layers = l, quadrants = q, secondaries = secs, clusters = cl, tot = tot, extra = ex, repetition = rep, save_plots = save_plots)
                                            comp.run()
                                            # except:
                                            #     pass

                                        else:
                                            print('Simulation type not valid!')
print('Done!')
