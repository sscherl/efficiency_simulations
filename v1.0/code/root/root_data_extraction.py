#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class to extract data from .root file, write to .csv file, plot stuff
#
#	Author: Sigrid Scherl
#
#	Created: June 2022
#

import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas

import csv
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
from pathlib import Path
import random
from random import randint
import sys
import time

from toolkit import ProgressBar, Constants
from simulation_data_plot import RootHisto

###### Class ####################################################

class RootData:

    def __init__(self, data_folder, histo_folder, plot_folder, data_file, file_name):
        self.data_folder = data_folder
        self.histo_folder = histo_folder
        self.plot_folder = plot_folder
        self.data_file = data_file
        self.file_name = file_name
        
        self.const = Constants()

        if self.data_file[-4:] == '.csv': self.get_csv_data()
        elif self.data_file[-5:] == '.root': self.get_root_data()
        self.get_sensors()

    
    def get_csv_data(self):

        data = pd.read_csv(self.data_folder+self.data_file)
        self.df = pd.DataFrame(data)


    def get_root_data(self):

        file = TFile(self.data_folder+self.data_file)
        tree = file.Get('SensorHits')
        
        self.entries = tree.GetEntries()

        evt_num, layer, quadrant, run_num, sensor_num, tof, unq_evt_num, x, x0, x_symloc, y, y0, y_symloc, z = [], [], [], [], [], [], [], [], [], [], [], [], [], []

        progress_bar = ProgressBar(self.entries)
        
        print('Extracting data from tree...')

        for i in range(0, self.entries):
            tree.GetEntry(i)

            # Check if new sensor:
            if i == 0:
                sensor_num.append(0)
            else:
                if tree.x0 == x0[-1] and tree.y0 == y0[-1]:
                    sensor_num.append(sensor_num[-1])
                else: 
                    sensor_num.append(sensor_num[-1] + 1)

            # Fill in original values
            run_num.append(tree.run_num)      # collection of simulated events
            evt_num.append(tree.evt_num)      # identifier of event within a given run
            layer.append(tree.layer)          # layer number of hit, only x-aligned layers included, 1-6 in ascending order with z position
            x.append(tree.x)                  # x position of hit in SciFi in global detector coordinate system
            y.append(tree.y)                  # y position of hit in SciFi in global detector coordinate system
            z.append(tree.z)                  # z position of hit in SciFi in global detector coordinate system
            x_symloc.append(tree.x_symloc)    # symmetrised local sensor coordinates, described above in preamble
            y_symloc.append(tree.y_symloc)    # symmetrical local sensor coordinates, described above in preamble
            x0.append(tree.x0)                # sensor x centre position
            y0.append(tree.y0)                # sensor y centre position

            # Calculate time of flight 
            tof.append(tree.z/(1000*self.const.c))              # time of flight is distance to detector (z) for speed of light

            # Determine quadrant
            quadrant.append(self.get_quadrant(tree.x, tree.y))

            # Calculate unique event number
            unq_evt_num.append(int(evt_num[-1] + 5 * (run_num[-1] - run_num[0]) + 5 * 100 * sensor_num[-1]))

            progress_bar.PrintBar(i)

        # Fill into dataframe
        self.df = pd.DataFrame(columns = ['Sensor', 'X0', 'Y0', 'UniqueEvent','Run', 'Event', 'Layer', 'X', 'Y', 'Z', 'X_symloc', 'Y_symloc', 'Quadrant', 'ToF'])

        self.df.Sensor = sensor_num
        self.df.X0 = x0
        self.df.Y0 = y0
        self.df.UniqueEvent = unq_evt_num
        self.df.Run = run_num
        self.df.Event = evt_num
        self.df.Layer = layer
        self.df.X = x
        self.df.Y = y
        self.df.Z = z
        self.df.X_symloc = x_symloc
        self.df.Y_symloc = y_symloc
        self.df.Quadrant = quadrant
        self.df.ToF = tof


    def get_data_rate():

        # Calculate hits per event
        df_count_events = self.df['Event'].value_counts().to_frame().reset_index()
        df_count_events.columns = ['Event', 'Count']

        # Omit contributions
        # df_count_events.drop(df_count_events.index[0], axis=0, inplace=True)

        self.data_rate = df_count_events.Count.sum()/math.ceil(self.df['Event'].iloc[-1])


    def get_quadrant(self, x, y, definition = 'lhcb'):

        if definition == 'lhcb': # Definition from LHCb SciFi Simulation Data: (x,y) : (+,+) := Quad 1 , (-,+) := Quad 2 , (-,-) := Quad 3 , (+,-) := Quad 4.
            if   x > 0 and y > 0: q = 1
            elif x < 0 and y > 0: q = 2
            elif x < 0 and y < 0: q = 3
            elif x > 0 and y < 0: q = 4
            else: print('Error! Quadrant not defined.')

            return q


    def get_sensors(self):

        self.df_sensors = self.df[['X0', 'Y0']].drop_duplicates()
    

    def get_sensor(self, nr):

        df_sensor = self.df[self.df.X0 == self.df_sensors.iloc[nr]['X0']]
        df_sensor = df_sensor[df_sensor.Y0 == self.df_sensors.iloc[nr]['Y0']]

        return df_sensor


    def plot_root_histos(self):
        
        folder = self.histo_folder
        data_type = 'LHCb Simulation Data'

        print('Creating histograms...')

        RootHisto(dir = folder, plot_name = 'histo_run_num', plot_title = data_type, x_axis_title = 'Run Number', x_values = self.df.Run, bins = 100, min = 4990, max = 5110)
        RootHisto(dir = folder, plot_name = 'histo_evt_num', plot_title = data_type, x_axis_title = 'Event Number', x_values = self.df.Event, bins = 60, min = 0, max = 6)
        RootHisto(dir = folder, plot_name = 'histo_layer', plot_title = data_type, x_axis_title = 'Layer', x_values = self.df.Layer, bins = 70, min = 0, max = 7)
        
        RootHisto(dir = folder, plot_name = 'histo_x', plot_title = data_type, x_axis_title = 'X (mm)', x_values = self.df.X, bins = 100, min = -3000, max = 3000)
        RootHisto(dir = folder, plot_name = 'histo_y', plot_title = data_type, x_axis_title = 'Y (mm)', x_values = self.df.Y, bins = 100, min = -700, max = 700)
        RootHisto(dir = folder, plot_name = 'histo_z', plot_title = data_type, x_axis_title = 'Z (mm)', x_values = self.df.Z, bins = 100, min = 7700, max = 9500)
        
        RootHisto(dir = folder, plot_name = 'histo_x_symloc', plot_title = data_type, x_axis_title = 'X_symloc (mm)',x_values = self.df.X_symloc, bins = 100, min = -15, max = 15)
        RootHisto(dir = folder, plot_name = 'histo_y_symloc', plot_title = data_type, x_axis_title = 'Y_symloc (mm)',x_values = self.df.Y_symloc, bins = 100, min = -15, max = 15)
        
        RootHisto(dir = folder, plot_name = 'histo_x0', plot_title = data_type, x_axis_title = 'X0 (mm)', x_values = self.df.X0, bins = 100, min = 0, max = 2500)
        RootHisto(dir = folder, plot_name = 'histo_y0', plot_title = data_type, x_axis_title = 'Y0 (mm)', x_values = self.df.Y0, bins = 100, min = 0, max = 700)
        

    def plot_sensor_positions(self, save = True, show = True):

        plt.scatter(self.df_sensors.X0, self.df_sensors.Y0, marker = '.')
        plt.xlabel('Sensor x centre position (mm)')
        plt.ylabel('Sensor y centre position (mm)')

        if save: plt.savefig(self.plot_folder+'sensor_centre_positions.pdf')
        if show: plt.show()


    def plot_sensor_xy(self, sensor, save = True, show = True): # Plot hits of single sensor - all original quadrants

        sensor_x0 = int(sensor.iloc[0]['X0'])
        sensor_y0 = int(sensor.iloc[0]['Y0'])

        x_min, x_max = sensor_x0 - 13, sensor_x0 + 13
        y_min, y_max = sensor_y0 - 13, sensor_y0 + 13

        fig, ax = plt.subplots(2, 2) # [row, col]
        ax[0,0].scatter(sensor.X, sensor.Y, marker = '.')
        ax[1,0].scatter(sensor.X, sensor.Y, marker = '.')
        ax[0,1].scatter(sensor.X, sensor.Y, marker = '.')
        ax[1,1].scatter(sensor.X, sensor.Y, marker = '.')

        ax[0,0].set_xlim(-x_max, -x_min)
        ax[0,0].set_ylim(y_min, y_max)
        ax[1,0].set_xlim(-x_max, -x_min)
        ax[1,0].set_ylim(-y_max, -y_min)
        ax[0,1].set_xlim(x_min, x_max)
        ax[0,1].set_ylim(y_min, y_max)
        ax[1,1].set_xlim(x_min, x_max)
        ax[1,1].set_ylim(-y_max, -y_min)

        ax[0,0].spines['bottom'].set_visible(False)
        ax[0,0].spines['right'].set_visible(False)
        ax[0,1].spines['bottom'].set_visible(False)
        ax[0,1].spines['left'].set_visible(False)
        ax[1,0].spines['top'].set_visible(False)
        ax[1,0].spines['right'].set_visible(False)
        ax[1,1].spines['top'].set_visible(False)
        ax[1,1].spines['left'].set_visible(False)

        ax[0,0].xaxis.tick_top()
        ax[0,1].xaxis.tick_top()
        ax[0,0].tick_params(labeltop = False)
        ax[0,1].tick_params(left = False, labeltop = False, labelleft = False)
        ax[1,1].tick_params(left = False, labelleft = False)
        ax[1,0].xaxis.tick_bottom()
        ax[1,1].xaxis.tick_bottom()
        
        d = .5  # proportion of vertical to horizontal extent of the slanted line
        kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12, linestyle="none", color='k', mec='k', mew=1, clip_on=False)
        ax[0,0].plot([0, 1], [0, 1], transform=ax[0,0].transAxes, **kwargs)
        ax[0,1].plot([0, 1], [1, 0], transform=ax[0,1].transAxes, **kwargs)
        ax[1,0].plot([0, 1], [1, 0], transform=ax[1,0].transAxes, **kwargs)
        ax[1,1].plot([0, 1], [0, 1], transform=ax[1,1].transAxes, **kwargs)

        fig.text(0.5, 0.02, 'X (mm)', ha = 'center')
        fig.text(0.02, 0.5, 'Y (mm)', va = 'center', rotation='vertical')
        fig.text(0.5, 0.92, 'Hits for sensor at X0 = '+str(sensor_x0)+' mm, Y0 = '+str(sensor_y0)+' mm', ha = 'center')

        if save: plt.savefig(self.plot_folder+'sensorX0'+str(sensor_x0)+'Y0'+str(sensor_y0)+'_xy.pdf')
        if show: plt.show()


    def plot_sensor_xyz(self, sensor, save = True, show = True):  #3D plot hits of single sensor - mirrored quadrants

        sensor_x0 = int(sensor.iloc[0]['X0'])
        sensor_y0 = int(sensor.iloc[0]['Y0'])

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        ax.scatter3D(sensor.X_symloc, sensor.Z, sensor.Y_symloc, c=sensor.Layer, cmap='jet')
        ax.set_xlabel('X_symloc (mm)')
        ax.set_ylabel('Z (mm)')
        ax.set_zlabel('Y_symloc (mm)')
        fig.text(0.5, 0.92, '3D hits for sensor at X0 = '+str(sensor_x0)+' mm, Y0 = '+str(sensor_y0)+' mm', ha = 'center')

        if save: plt.savefig(self.plot_folder+'sensorX0'+str(sensor_x0)+'Y0'+str(sensor_y0)+'_xyz.pdf')
        if show: plt.show()


    def plot_tof(self, l, save = True, show = True):

        if l == all:
            plt.hist(self.df.ToF*1e9, bins= 1000) # Want time of flight in ns, before given in s
            plt.title('All layers')
            plot_title = 'tof.pdf'
        else:
            layers = {'Number': [1,2,3,4,5,6], 'Min': [26, 26.5, 28, 29, 30.5, 31], 'Max': [26.5, 27, 28.5, 29.5, 31, 31.5]}
            df_layers = pd.DataFrame.from_dict(layers)
            df_tof_layer = self.df[self.df['ToF'].between(df_layers.iloc[l]['Min'], df_layers.iloc[l]['Max'])]
            plt.hist(df_tof_layer.ToF*1e9, bins= 100) # Want time of flight in ns, before given in s
            plt.title('Layer '+str(l))
            plot_title = 'tof_layer'+str(l)+'.pdf'

        plt.xlabel('Time of flight (ns)')
        plt.ylabel('Counts')

        if save: plt.savefig(self.plot_folder+plot_title)
        if show: plt.show()


    def save_data(self):

        print('Saving data file...')

        df2 = pd.DataFrame(columns = ['Event', 'Layer', 'Quadrant', 'X', 'Y', 'ToF', 'Energy'])

        df2.Event = self.df.UniqueEvent
        df2.Layer = self.df.Layer
        df2.Quadrant = self.df.Quadrant
        df2.X = self.df.X_symloc
        df2.Y = self.df.Y_symloc
        df2.ToF = self.df.ToF
        df2.Energy = np.zeros(self.entries)

        df2.to_csv(self.data_folder+self.file_name+'_Data.csv', index = False)


    def save_all_data(self):

        print('Saving debug data file...')

        self.df.to_csv(self.data_folder+self.file_name+'_DataAll.csv', index = False)


###### Run ######################################################

data_folder = '../data_v2/raw_data/'
histo_folder = '../plots_v2/root_data_histograms/'
plot_folder = '../plots_v2/root_data_plots/'
root_file = 'MagUp_bs2phiphi_violaine_energytime_oscar_files_1p5e34.root_Sorted_Sensors.root'
csv_file = 'MagUp_bs2phiphi_violaine_energytime_oscar_files_1p5e34_Sorted_Sensors_DataAll.csv'
file_name = 'MagUp_bs2phiphi_violaine_energytime_oscar_files_1p5e34_Sorted_Sensors'

data = RootData(data_folder = data_folder, histo_folder = histo_folder, plot_folder = plot_folder, data_file = root_file, file_name = file_name)

data.plot_sensor_xyz(sensor = data.get_sensor(1), save = 'False')