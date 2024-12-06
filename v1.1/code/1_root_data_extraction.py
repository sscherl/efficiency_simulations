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
from matplotlib.pyplot import figure
import matplotlib
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
from pathlib import Path
import random
from random import randint
import sys
import time
import math
import yaml

from classes.toolkit import ProgressBar, Constants
from classes.simulation_data_plot import RootHisto

###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Class ####################################################

class RootData:

    def __init__(self, data_source):

        self.data_source = data_source

        if data_source == 'csv':
            self.data_file = config['files']['raw_data']['csv']
            self.get_csv_data()

        elif data_source == 'root':
            self.data_file = config['files']['raw_data']['root']
            self.get_root_data()

        self.get_sensors()

        matplotlib.rcParams.update({'font.size': 19})
        # matplotlib.rcParams['text.usetex'] = True

    # Internal use only
    def get_csv_data(self):

        data = pd.read_csv(config['directories']['data']['raw']+self.data_file)
        self.df = pd.DataFrame(data)

    # Internal use only
    def get_root_data(self):

        # Read file and get tree
        file = TFile(config['directories']['data']['raw']+self.data_file)
        tree = file.Get('SensorHits')  
        self.entries = tree.GetEntries()

        # Initialise varaibles for leaves
        self.evt_num, self.layer, self.quadrant, self.run_num, self.sensor_num, self.tof, self.unq_evt_num, self.x, self.x0, self.x_symloc, self.y, self.y0, self.y_symloc, self.z = [], [], [], [], [], [], [], [], [], [], [], [], [], []

        progress_bar = ProgressBar(self.entries)
        
        # Get tree content
        print('Extracting data from tree...')
        for i in range(0, self.entries):
            tree.GetEntry(i)

            # Check if new sensor:
            if i == 0:
                self.sensor_num.append(0)
            else:
                if tree.x0 == self.x0[-1] and tree.y0 == self.y0[-1]:
                    self.sensor_num.append(self.sensor_num[-1])
                else: 
                    self.sensor_num.append(self.sensor_num[-1] + 1)

            # Fill in original values
            self.run_num.append(tree.run_num)      # collection of simulated events
            self.evt_num.append(tree.evt_num)      # identifier of event within a given run
            self.layer.append(tree.layer)          # layer number of hit, only x-aligned layers included, 1-6 in ascending order with z position
            self.x.append(tree.x)                  # x position of hit in SciFi in global detector coordinate system
            self.y.append(tree.y)                  # y position of hit in SciFi in global detector coordinate system
            self.z.append(tree.z)                  # z position of hit in SciFi in global detector coordinate system
            self.x_symloc.append(tree.x_symloc)    # symmetrised local sensor coordinates, described above in preamble
            self.y_symloc.append(tree.y_symloc)    # symmetrical local sensor coordinates, described above in preamble
            self.x0.append(tree.x0)                # sensor x centre position
            self.y0.append(tree.y0)                # sensor y centre position
            self.tof.append(tree.t_coll)           # time of collision = tof

            # Calculate time of flight 
            #self.tof.append(tree.z*1e6/(self.const.c))              # time of flight is distance to detector (z) for speed of light - given in ns here

            # Determine quadrant
            self.quadrant.append(self.get_quadrant(tree.x, tree.y))

            # Calculate unique event number (my new version) - all sensors have same events (500 events)
            self.unq_evt_num.append(int(self.evt_num[-1] + 5 * (self.run_num[-1] - self.run_num[0])))

            
            # Calculate unique event number (my old version) - count new sensor as new event to maximise number of events
            #self.unq_evt_num.append(int(self.evt_num[-1] + 5 * (self.run_num[-1] - self.run_num[0]) + 5 * 100 * self.sensor_num[-1]))

            # # Calculate unique event number (violaine's version)
            # self.unq_evt_num.append(int(self.evt_num[-1] + 100 * (self.run_num[-1])))

            progress_bar.PrintBar(i)

        # Fill into dataframe
        self.df = pd.DataFrame(columns = ['Sensor', 'X0', 'Y0', 'UniqueEvent','Run', 'Event', 'Layer', 'X', 'Y', 'Z', 'X_symloc', 'Y_symloc', 'Quadrant', 'ToF'])

        self.df.Sensor = self.sensor_num
        self.df.X0 = self.x0
        self.df.Y0 = self.y0
        self.df.UniqueEvent = self.unq_evt_num
        self.df.Run = self.run_num
        self.df.Event = self.evt_num
        self.df.Layer = self.layer
        self.df.X = self.x
        self.df.Y = self.y
        self.df.Z = self.z
        self.df.X_symloc = self.x_symloc
        self.df.Y_symloc = self.y_symloc
        self.df.Quadrant = self.quadrant
        self.df.ToF = self.tof

    # DOES NOT WORK YET
    def get_data_rate(self):

        # Calculate hits per event
        df_count_events = self.df['Event'].value_counts().to_frame().reset_index()
        df_count_events.columns = ['Event', 'Count']

        # Omit contributions
        # df_count_events.drop(df_count_events.index[0], axis=0, inplace=True)

        self.data_rate = df_count_events.Count.sum()/math.ceil(self.df['Event'].iloc[-1])
        print(self.data_rate)

    # Internal use only
    def get_quadrant(self, x, y, definition = 'lhcb'):

        if definition == 'lhcb': # Definition from LHCb SciFi Simulation Data: (x,y) : (+,+) := Quad 1 , (-,+) := Quad 2 , (-,-) := Quad 3 , (+,-) := Quad 4.
            if   x > 0 and y > 0: q = 1
            elif x < 0 and y > 0: q = 2
            elif x < 0 and y < 0: q = 3
            elif x > 0 and y < 0: q = 4
            else: print('Error! Quadrant not defined.')

            return q

    # Internal use only
    def get_sensors(self):

        self.df_sensors = self.df[['Sensor', 'X0', 'Y0']].drop_duplicates()


    # Get centre coordinate of sensor by position on list
    def get_sensor(self, nr):

        df_sensor = self.df[self.df.Sensor == nr]
        df_sensor = df_sensor.reset_index(drop=True)

        return df_sensor


    # pyplot histos of the raw root data
    def plot_histo(self, plot_name, x, x_label, bins, y_label = 'Counts', x_lim = '', y_lim = ''):

        fig = plt.figure(figsize=(13,6)) #
        fig.set_tight_layout(True)

        plt.hist(x, bins = bins)
        plt.ticklabel_format(style='sci', axis='y', scilimits=(3,3))

        plt.xlabel(x_label)
        plt.ylabel(y_label)

        plt.yscale('log')

        if x_lim != '': plt.xlim(x_lim)
        if y_lim != '': plt.ylim(y_lim)

        plt.tight_layout()

        plt.savefig(config['directories']['plots']['root_histos']+plot_name+'.pdf')
        plt.show()
        plt.clf()
        plt.close('all')

        try: fig.clf()
        except: pass

    def plot_new_histos(self):

        self.df_sel = self.df[self.df['Layer'] == 2] #1843
        self.df_sel_min = self.df_sel[self.df_sel['ToF'] < 26.2]
        print (self.df_sel_min)
        #print(self.df_sel.ToF.min())
        #print(self.df_sel.ToF.max())
        #self.df_sel = self.df_sel[self.df_sel['ToF'] > 20]
        #df_selected = self.df_sel['Y0'].value_counts().to_frame().reset_index()
        #df_hits_per_evt.columns = ['UniqueEvent', 'Hits']

        #self.plot_histo(plot_name = 'wide_histo_events', x = self.df.UniqueEvent, x_label = 'Event',  bins = 500, y_label = 'Hits')
        #self.plot_histo(plot_name = 'wide_histo_z',      x = self.df.Z,           x_label = 'Z (mm)', bins = 1700, y_label = 'Hits', x_lim = (7750,9450))
        #self.plot_histo(plot_name = 'wide_histo_x',      x = self.df.X,           x_label = 'X (mm)', bins = 4768, y_label = 'Hits')
        #self.plot_histo(plot_name = 'wide_histo_y',      x = self.df.Y,           x_label = 'Y (mm)', bins = 1204, y_label = 'Hits')
        #self.plot_histo(plot_name = 'wide_histo_tof',      x = self.df.ToF,        x_label = 'Time of Flight (ns)', bins = 1200, y_label = 'Hits')

    
    def plot_histos(self):
        self.plot_histo(plot_name = 'histo_run_num',  x = self.df.Run,      x_label = 'Run Number',    bins = 100)
        self.plot_histo(plot_name = 'histo_evt_num',  x = self.df.Event,    x_label = 'Event Number',  bins = 60)
        self.plot_histo(plot_name = 'histo_layer',    x = self.df.Layer,    x_label = 'Layer',         bins = 70)
        self.plot_histo(plot_name = 'histo_x',        x = self.df.X,        x_label = 'X (mm)',        bins = 500)
        self.plot_histo(plot_name = 'histo_y',        x = self.df.Y,        x_label = 'Y (mm)',        bins = 500)
        self.plot_histo(plot_name = 'histo_z',        x = self.df.Z,        x_label = 'Z (mm)',        bins = 70)
        self.plot_histo(plot_name = 'histo_x_symloc', x = self.df.X_symloc, x_label = 'X_symloc (mm)', bins = 100)
        self.plot_histo(plot_name = 'histo_y_symloc', x = self.df.Y_symloc, x_label = 'Y_symloc (mm)', bins = 100)
        self.plot_histo(plot_name = 'histo_x0',       x = self.df.X0,       x_label = 'X0 (mm)',       bins = 120)
        self.plot_histo(plot_name = 'histo_y0',       x = self.df.Y0,       x_label = 'Y0 (mm)',       bins = 120)

    # Plot the centre positions of all sensors
    def plot_sensor_positions(self, sensor = '', save = True, show = True):

        plt.scatter(self.df_sensors.X0, self.df_sensors.Y0, marker = '.')
        df_y_beampipe = self.df_sensors[self.df_sensors.X0 == 14.0]
        print(df_y_beampipe.Y0.min())

        if sensor != '':
            starred_sensor = self.get_sensor(sensor)
            starred_X0 = int(starred_sensor.iloc[0]['X0'])
            starred_Y0 = int(starred_sensor.iloc[0]['Y0'])
            plt.scatter(starred_X0, starred_Y0, marker = '*', s = 50, color = 'black')

        plt.xlabel('Sensor x centre position (mm)')
        plt.ylabel('Sensor y centre position (mm)')
        plt.tight_layout()

        if save and sensor != '': plt.savefig(config['directories']['plots']['root_plots']+'sensor_centre_positions_X0'+str(starred_X0)+'Y0'+str(starred_Y0)+'.pdf')
        elif save and sensor == '': plt.savefig(config['directories']['plots']['root_plots']+'sensor_centre_positions.pdf')
        if show: plt.show()

    # Plot the centre positions of all sensors in all four quadrants
    def plot_sensor_positions_4Q(self, save = True, show = True, zoom = False):

        df_sensors_Q2 = self.df_sensors.copy()
        df_sensors_Q2.X0 =  df_sensors_Q2.X0 * -1

        df_sensors_Q3 = self.df_sensors.copy()
        df_sensors_Q3.X0 = df_sensors_Q3.X0 * -1
        df_sensors_Q3.Y0 = df_sensors_Q3.Y0 * -1

        df_sensors_Q4 = self.df_sensors.copy()
        df_sensors_Q4.Y0 =  df_sensors_Q4.Y0 * -1

        if zoom:
            matplotlib.rcParams.update({'font.size': 16})
            figure(figsize=(4.8, 4.5), dpi=100)
            size = 70
        else:
            figure(figsize=(13, 4.5), dpi=100)
            size = 10

        plt.scatter(self.df_sensors.X0, self.df_sensors.Y0, marker = 's', s = size)
        plt.scatter(df_sensors_Q2.X0, df_sensors_Q2.Y0, marker = 's', s = size)
        plt.scatter(df_sensors_Q3.X0, df_sensors_Q3.Y0, marker = 's', s = size)
        plt.scatter(df_sensors_Q4.X0, df_sensors_Q4.Y0, marker = 's', s = size)

        plt.xlabel('Sensor position in x (mm)')
        plt.ylabel('Sensor position in y (mm)')

        if zoom:
            plt.xlim(-220,220)
            plt.ylim(-220,220)
        else:
            plt.xlim(-2500,2500)
            plt.ylim(-650,650)

        plt.tight_layout()


        if save:
            if zoom: plt.savefig(config['directories']['plots']['root_plots']+'sensor_centre_positions_4Q_zoom.pdf')
            else: plt.savefig(config['directories']['plots']['root_plots']+'sensor_centre_positions_4Q.pdf')
        if show: plt.show()


    # Plot hits of single sensor - all four original quadrants
    def plot_sensor_xy(self, sensor, save = True, show = True):

        sensor_id = int(sensor.iloc[0]['Sensor'])
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

        plt.tight_layout()

        if save: plt.savefig(config['directories']['plots']['root_plots']+'sensor'+str(sensor_id)+'_xy.pdf')
        if show: plt.show()


    #3D plot hits of single sensor - mirrored quadrants
    def plot_sensor_xyz(self, sensor, save = True, show = True):

        sensor_id = int(sensor.iloc[0]['Sensor'])
        sensor_x0 = int(sensor.iloc[0]['X0'])
        sensor_y0 = int(sensor.iloc[0]['Y0'])

        fig = plt.figure()
        ax = plt.axes(projection='3d')

        ax.scatter3D(sensor.X_symloc, sensor.Z, sensor.Y_symloc, c=sensor.Layer, cmap='jet')
        ax.set_xlabel('X_symloc (mm)')
        ax.set_ylabel('Z (mm)')
        ax.set_zlabel('Y_symloc (mm)')
        fig.text(0.5, 0.92, '3D hits for sensor '+str(sensor_id)+' at X0 = '+str(sensor_x0)+' mm, Y0 = '+str(sensor_y0)+' mm', ha = 'center')

        plt.tight_layout()

        if save: plt.savefig(config['directories']['plots']['root_plots']+'sensor'+str(sensor_id)+'_xyz.pdf')
        if show: plt.show()


    # Plot the time of flight of single layer (1 to 6) or all layers
    def plot_tof(self, l = 'all', save = True, show = True):

        fig = plt.figure()
        fig.set_tight_layout(True)

        if l == 'all':
            plt.hist(self.df.ToF*1e9, bins= 1000) # Want time of flight in ns, before given in s
            plt.title('All layers')
            plot_title = 'tof.pdf'
        else:
            layers = {'Number': [1,2,3,4,5,6], 'Min': [26e-9, 26.5e-9, 28e-9, 29e-9, 30.5e-9, 31e-9], 'Max': [26.5e-9, 27e-9, 28.5e-9, 29.5e-9, 31e-9, 31.5e-9]}
            df_layers = pd.DataFrame.from_dict(layers)
            df_tof_layer = self.df[self.df['ToF'].between(df_layers.iloc[l-1]['Min'], df_layers.iloc[l-1]['Max'])]
            plt.hist(df_tof_layer.ToF, bins= 100)
            plt.title('Layer '+str(l))
            plot_title = 'tof_layer'+str(l)+'.pdf'

        plt.xlabel('Time of flight (ns)')
        plt.ylabel('Counts')
        plt.tight_layout()

        if save: plt.savefig(config['directories']['plots']['root_plots']+plot_title)
        if show: plt.show()


    def print_secs(self, sensor = '', layer = 1, quadrant = 1):

        if sensor != '': df = self.get_sensor(sensor)
        else: df = self.df

        total_nr = len(df)
        df_tof = df[df['ToF'] > 35]
        df_l1_tof = df_tof[df_tof['Layer'] == layer]
        df_q1_l1_tof = df_l1_tof[df_l1_tof['Quadrant'] == quadrant]
        secs_nr = len(df_q1_l1_tof)

        for val in df_l1_tof['ToF']:
            if val > 35: print(val)

        print('Total hits: {}, Secondaries: {}'.format(total_nr, secs_nr))


    # Save the data to a csv file, to simplify further usage
    def save_selected_data_to_csv(self):

        print('Saving data file...')

        df2 = pd.DataFrame(columns = ['Event', 'Layer', 'Quadrant', 'X', 'Y', 'ToF', 'Energy'])

        df2.Event = self.df.UniqueEvent
        df2.Layer = self.df.Layer
        df2.Quadrant = self.df.Quadrant
        df2.X = self.df.X_symloc
        df2.Y = self.df.Y_symloc
        df2.ToF = self.df.ToF
        df2.Energy = np.zeros(len(df2))

        df2.to_csv(config['directories']['data']['raw']+config['files']['raw_data']['selected'], index = False)
    

    # Save the data for one sensor to a csv file, to simplify further usage
    def save_sensor_data(self, nr):

        sensor = self.get_sensor(nr)
        sensor_id = int(sensor.iloc[0]['Sensor'])

        print('Saving data file for sensor ' +str(sensor_id)+'...')
        
        df2 = pd.DataFrame(columns = ['Event', 'Layer', 'Quadrant', 'X', 'Y', 'ToF', 'Energy'])

        df2.Event = (sensor.UniqueEvent - sensor.UniqueEvent.iloc[0] + 1).astype(int)
        df2.Layer = sensor.Layer.astype(int)
        df2.Quadrant = sensor.Quadrant.astype(int)
        df2.X = sensor.X_symloc
        df2.Y = sensor.Y_symloc
        df2.ToF = sensor.ToF
        df2.Energy = np.zeros(len(df2))

        df2.to_csv(config['directories']['data']['raw_sensor']+config['files']['raw_data']['name']+'_Sensor'+str(sensor_id)+'.csv', index = False)


    # Save data of all sensors
    def save_all_sensor_data(self):

        sensor_list = self.df_sensors.Sensor.tolist()

        for sensor in sensor_list:
            self.save_sensor_data(sensor)


    # Save all the data available to a csv file
    def save_root_data_to_csv(self):

        if self.data_source == 'root':
            print('Saving data file...')
            self.df.to_csv(config['directories']['data']['raw']+config['files']['raw_data']['csv'], index = False)
        else:
            print('ERROR: Can only create csv file for all sensors from root file!')


    def save_sensor_list(self):

        df_sensors_sorted = self.df_sensors.sort_values(['Y0', 'X0'], ascending = [True,True]) 
        df_sensors_sorted.to_csv(config['directories']['data']['raw']+config['files']['raw_data']['sensor_list'], index = False)

    # Generate list of sensors to be filled with total and missing hits a la 'Sensor_Efficiencies_LXX_QXX.csv' in data/log_files
    def save_empty_sensor_efficiencies_list(self):

        df_sensors_sorted = self.df_sensors.sort_values(['Y0', 'X0'], ascending = [True,True])
        df_sensors_sorted['TotalHits'] = [''] * len(df_sensors_sorted)
        df_sensors_sorted['MissingHits'] = [''] * len(df_sensors_sorted)
        df_sensors_sorted.to_csv(config['directories']['data']['efficiencies']+config['files']['efficiencies']['empty'], index = False)


###### Run ######################################################

data = RootData(data_source = 'csv')
#data.save_root_data_to_csv()
#data.plot_new_histos()
#sensor = data.get_sensor(739)
#data.plot_sensor_positions()
data.plot_sensor_positions_4Q(zoom = True)
