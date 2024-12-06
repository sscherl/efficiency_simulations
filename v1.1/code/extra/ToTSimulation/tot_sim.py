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


class Plot:

    def __init__(self):

        self.plot_name       = 'MightyPix1_SimulatedToT'
        self.measurement_dir = ''
        self.data_dir        = '' 
        self.data_file       = 'mp1_recout.csv' 

        # Plot settings
        self.plot_lines      = False
        self.plot_logy       = False

        # Font size for plot labels
        matplotlib.rcParams.update({'font.size': 12})

        self.marker_counter = 0
        self.markers = ['v', '^', '<', '>', 'o']

        self.colour_counter = 0
        self.colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']


    # Get data of chip
    def get_data(self):

        self.df_data = pd.DataFrame(pd.read_csv(self.data_file))

        self.df_data.columns = ['3p5e-7_X', '3p5e-7_Y', '3p5e-8_X','3p5e-8_Y', '8p5e-8_X', '8p5e-8_Y']
        print(self.df_data)


    def plot(self):
        
        plt.plot(self.df_data['3p5e-8_X'] * 1e6, self.df_data['3p5e-8_Y'], label = '  2185 e-')
        plt.plot(self.df_data['8p5e-8_X'] * 1e6, self.df_data['8p5e-8_Y'], label = '  5305 e-')
        plt.plot(self.df_data['3p5e-7_X'] * 1e6, self.df_data['3p5e-7_Y'], label = '21850 e-')

        plt.yticks(np.arange(0, 2, step=0.3))
        plt.xlim(-0.1,3)


        plt.xlabel('Time (μs)')
        plt.ylabel('Voltage (V)')



    # Plot all added IVs
    def wrapup(self):

        # plt.tight_layout()
        plt.legend()

        if self.plot_logy: self.plot_name += '_log'

        plt.savefig(self.plot_name+'.pdf')
        plt.show()
            

    # Run all functions
    def run(self):

        self.get_data()
        self.plot()
        self.wrapup()



  

###### IVs with PCB without chip for background - 5 measurements per value ######
        
thisplot = Plot()

thisplot.run()



