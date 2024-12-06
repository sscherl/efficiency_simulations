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

from classes.simulation_data_comparison import DataComparison
from classes.toolkit import Chip, Helper, File

###### Run options ##############################################

# Mute SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'

###### Functions ################################################

helper = Helper()
mightypix = Chip()

def comp_sim(quadrant, secondaries, clusters, clockspeed, save_plots, show_plots):

    if not secondaries: secs = '_wosecs'
    else: secs = ''

    if clusters: clstr = '_clusters'
    else: clstr = ''

    comp = DataComparison(mightypix, 'C4785x17600um2_P165x55um2_E1to500_L1to6_app_Q'+quadrant+secs+clstr, str(clockspeed)+'MHz')

    comp.compare()
    comp.plots(save_me = save_plots, show_me = show_plots)
    comp.end()


def comp_rnd(event, rate, clusters, clockspeed, save_plots, show_plots):

    if clusters: clstr = '_clusters'
    else: clstr = ''

    comp = DataComparison(mightypix, 'rnd_C4785x17600um2_P165x55um2_E1to'+event+'_L1_Q1_R'+rate+clstr, str(clockspeed)+'MHz')

    comp.compare()
    comp.plots(save_me = save_plots, show_me = show_plots)
    comp.end()


def comp_all():

    for clockspeed in [40, 160]:
        for quadrant in ['1', '1to2', '1to3', '1to4']:
            for secondaries in [True, False]:
                for clusters in [True, False]:
                    try: comp_sim(quadrant, secondaries, clusters, clockspeed, save_plots = True, show_plots = False)
                    except: pass

    for clockspeed in [40, 160]:
        for rate in ['1-7', '3-4', '5-1', '6-8']:
            for event in ['3000', '10000', '100000']:
                for clusters in [True, False]:
                    try: comp_rnd(event, rate, clusters, clockspeed, save_plots = True, show_plots = False)
                    except: pass

###### Run ######################################################
# comp_sim('1', True, False, 40, False, False)
# comp_rnd('3000', '1-7', False, 40, False, False)

print('Done!')
