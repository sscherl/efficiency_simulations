#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Process efficiency measurements
#
#	Author: Sigrid Scherl
#
#	Created: March 2024
#

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path
import random
import sys
import time
import math
import csv

from classes.simulation_data_plot import Plot
from classes.toolkit import File


class PixelMap:

    def __init__(self):

        self.data_path = ''
        self.data_files = []

        # Constants
        self.nr_of_cols = 29
        self.nr_of_rows = 320
        self.active_are = self.nr_of_cols * 165e-6 * self.nr_of_rows * 55e-6 # in m^2

        self.clock_period = 25e-9 #s

        self.data_set_counter = 0


    # Set paramters
    def set_parameters(self, inj_pixels, inj_freq_clk_cy = '', inj_hit_rate = ''):

        # Input parameters
        self.inj_pixels = inj_pixels
        self.nr_of_inj_pixels = len(self.inj_pixels)

        self.inj_freq_clk = inj_freq_clk_cy # in 40 MHz clock cycles
        self.hit_rate_MHz_cm = inj_hit_rate
        if self.hit_rate_MHz_cm != '':
            self.hit_rate_Hz_m = self.to_Hz_m(self.hit_rate_MHz_cm)

        if self.inj_freq_clk != '':
            self.calc_inj_hit_rate()
        elif inj_hit_rate != '':
            self.calc_inj_freq()
        
        if type(self.inj_freq_clk) != int:
            print('*** WARNING:\t Injection frequency ', self.inj_freq_clk, ' (clk) has to be multiple of clock period ', self.clock_period,'(s). Rounding to nearest integer.')
            self.inj_freq_clk = round(self.inj_freq_clk)
            print('*** INFO:\t New injection frequency is ', self.inj_freq_clk, '(clk).')
            self.calc_inj_hit_rate()
            print('*** INFO:\t New injection hit rate is', self.hit_rate_MHz_cm,'(MHz/cm^2).')
        else:
            print('*** INFO:\t New injection frequency is ', self.inj_freq_clk, '(clk).')
            print('*** INFO:\t New injection hit rate is', self.hit_rate_MHz_cm,'(MHz/cm^2).')



    # Get parameters
    def calc_inj_hit_rate(self):
        if self.inj_freq_clk == '':
            print('ERROR: Injection Frequency not provided.')
        self.inj_freq_s = self.inj_freq_clk * self.clock_period
        self.hit_rate_Hz_m = self.nr_of_inj_pixels / (self.active_are * self.inj_freq_s) # Hz/m
        self.hit_rate_MHz_cm = self.to_MHz_cm(self.hit_rate_Hz_m)

    def calc_inj_freq(self):
        if self.hit_rate_Hz_m == '':
            print('ERROR: Hit rate not provided.')
        self.inj_freq_s = self.nr_of_inj_pixels / (self.active_are * self.hit_rate_Hz_m)
        self.inj_freq_clk = self.inj_freq_s / self.clock_period

    def to_MHz_cm(self,inp):
        return inp / 1e10
    
    def to_Hz_m(self,inp):
        return inp * 1e10
    


    # Get data
    def add_dataset(self, path, basename = '', file_date = '', file_time = ''):
        self.data_set_counter += 1
        self.path = path
        if basename == '':
            if file_date != '' and file_time != '':
                self.basename = file_date + '_' + file_time + '_injection_'
            else: 
                print('ERROR: No file specified. Please provide whole name or data and time.')
        else:
            self.basename = basename
        self.import_data()
    
    def import_hit_data(self):
        cols = ['event_number', 'col', 'row', 'ts', 'ts_n', 'ts2', 'ts3', 'timestamp', 'err']

        dataframes = []

        for f in self.data_files:
            filename = self.data_path + f + 'hit_data.csv'
            print(filename)

            df = pd.read_csv(filename, index_col=None, header=0, names=cols)
            dataframes.append(df)

        self.hit_data = pd.concat(dataframes, axis=0, ignore_index=True)

    def import_ts_data(self):
        ts_file = File(path=self.path, base_name = self.basename)
        ts_file.file_name += 'ts_data'
        ts_file.suffix = '.csv'
        ts_file.set_parameters()
        cols = ['event_number', 'timestamp', 'err']
        self.ts_data = pd.DataFrame(pd.read_csv(ts_file.whole, header=None, names=cols))
        
    def import_event_data(self):
        event_file = File(path=self.path, base_name = self.basename)
        event_file.file_name += 'event_data'
        event_file.suffix = '.csv'
        event_file.set_parameters()
        cols = ['event_number', 'col', 'row', 'toa', 'tot', 'tdc', 'err', 'timestamp']
        self.event_data = pd.DataFrame(pd.read_csv(event_file.whole, header=None, names=cols))
    
    def import_data(self):
        self.import_hit_data()
        #self.import_ts_data()
        #self.import_event_data()




###### Run ##############################
        
pm = PixelMap()

# Set injection parameters

pixels = [[5,100], [10,200], [15, 300], [20, 210], [25,110]] # [[col_pixel0,row_pixel0],..., [col_pixelN, row_pixelN]]
# For 5 pixels: 24 clk cycles = 9.9 MHz, 12 cc = 19.8 MHz, 8 cc = 29.7 MHz, 6 cc = 39.6 MHz
pm.set_parameters(inj_pixels=pixels, inj_freq_clk_cy = 20, inj_hit_rate = '')


# Specify data set

five_faulty_pixels_data = ['20240328_151150_injection_', '20240328_151337_injection_', '20240328_151508_injection_', '20240328_151933_injection_', '20240328_152148_injection_']
noisy_data = ['20240404_151328_injection_','20240404_150827_injection_']

# pm.data_path = '/Users/sigridscherl/Documents/GitRepos/mightypix_measurements/gecco-daq/output/selected_injection_scans/'
# pm.data_files = ['20240404_151220_injection_', '20240404_151117_injection_', '20240404_150718_injection_', '20240404_150618_injection_']
pm.data_path = '/Users/sigridscherl/Documents/GitRepos/MightyTracker/verification_framework/simulation_data/v1.1/code/extra/h5data/'
pm.data_files = ['20240617_151643_injection_']
pm.import_data()



# Make plot

plot_path = 'extra/h5data/'
file_date = '20240617_'
file_time = '151643_'

plt_file = File(path=plot_path, prefix = 'hitmap', base_name=(file_date+file_time+'injection'), suffix='.pdf') # Need to use file function for plot
plot = Plot(file=plt_file, save=True, show=True)

cols = pm.hit_data.col
rows = pm.hit_data.row
# cols = [i[0] for i in pixels]
# rows = [i[1] for i in pixels]

plot.hitmap_side_histos(x=cols,y=rows,xbins=np.arange(0,29), ybins=np.arange(0,320), vmax = '', col_max = '', row_max = '', ticks = [], cmap = 'viridis_r', vmin = 1, cmin = 1)

# col_vals = pm.hit_data.col.value_counts().to_frame().reset_index()
# col_vals.columns = ['Column', 'Count']
# print(col_vals)
# plt.scatter(col_vals.Column, col_vals.Count)
# plt.plot([0,30], [10000,10000], color = 'tab:orange')
# plt.xlim(0,30)
# plt.ylim(9000,11000)
# plt.xlabel('Column')
# plt.ylabel('Hits')
# plt.tight_layout()
# plt.show()