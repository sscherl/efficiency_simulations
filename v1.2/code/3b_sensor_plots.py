#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Graphs for all sensors in Mighty Tracker region
#
#	Author: Sigrid Scherl
#
#	Created: November 2021
#

import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import yaml
from tqdm import tqdm
import math


###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Class ####################################################

class SensorPlots:

    def __init__(self, layer, quadrant, newfsm, save_plots = True, show_plots = True):
        self.layer = layer
        self.quadrant = quadrant
        self.newfsm = newfsm
        self.save_plots = save_plots
        self.show_plots = show_plots

        # For file naming
        if self.layer == 1: self.l = str(self.layer)
        else: self.l = '1to'+str(self.layer)+'_app'

        if self.quadrant == 1: self.q = str(self.quadrant)
        else: self.q = '1to'+str(self.quadrant)

        # Plot labels
        self.hitrate_label = 'Hit Rate (MHz/cm²)'
        self.efficiency_label = 'Efficiency (%)'
        self.inefficiency_label = 'Inefficiency (%)'
        self.distance_label = 'Distance to IP (mm)'
        self.hits_label = 'Total Hits'

        # Read data
        data = pd.read_csv('../data/efficiencies/efficiencies_C4785x17600um2_P165x55um2_E1to500_L'+self.l+'_Q'+self.q+'_wosecs.csv')
        df = pd.DataFrame(data)
        self.dim = len(df)

        # Generate numbered postitions
        self.add_new_x(df)
        self.add_new_y(df)
        self.df=self.new_xy(df)

        # Size of data
        self.x_dim = int(max(self.df.XPos))
        self.y_dim = int(max(self.df.YPos))

        # Generate hit rates
        self.add_hit_rates()
        self.import_rnd_eff()
        self.add_sensor_distances()
        self.add_sensor_efficiencies()

###### Calculate sensor positions ###############################

    def new_x(self,old_x):
        new_x = self.df_new_x.index[self.df_new_x['X0'] == old_x].to_list()
        if len(new_x) == 1: return new_x[0]

    def new_y(self, old_y):
        new_y = self.df_new_y.index[self.df_new_y['Y0'] == old_y].to_list()
        if len(new_y) == 1: return new_y[0]

    def add_new_x(self, df):
        df_new_x = df.X0.value_counts().to_frame().reset_index()
        df_new_x.columns = ['X0', 'Count']
        df_new_x = df_new_x.sort_values(by = ['X0']).reset_index(drop=True).drop(columns=['Count'])
        self.df_new_x = df_new_x

    def add_new_y(self, df):
        df_new_y = df.Y0.value_counts().to_frame().reset_index()
        df_new_y.columns = ['Y0', 'Count']
        df_new_y = df_new_y.sort_values(by = ['Y0']).reset_index(drop=True).drop(columns=['Count'])
        self.df_new_y = df_new_y

    def new_xy(self, df):
        df['XPos'] = df['X0']
        df['YPos'] = df['Y0']
        for j in range(self.dim):
            i = float(j)
            df = df.replace({'XPos': {i: self.new_x(i)}})
            df = df.replace({'YPos': {i: self.new_y(i)}})
        return df


###### Calculate hit rates ######################################

    def hit_rate(self, hits):
        time = 500 * self.layer * 25e-9 #s -- time = events * layers * BX
        area = 29*0.0165*320*0.0055 #cm²  -- sensing area = pixel size * rows * cols
        return hits/(time * area * 1e6) #MHz/cm²

    def add_hit_rates(self):
        rates = []
        for i in range(self.dim):
            rates.append(self.hit_rate(self.df.iloc[i]['TotalHits']))
        self.df['HitRate'] = rates


###### Calculate sensor distance to interaction point ###########

    def add_sensor_distances(self):
        dist = []
        for i in range(self.dim):
            dist.append(np.sqrt(float(self.df.iloc[i]['X0'])**2 + float(self.df.iloc[i]['Y0'])**2))
        self.df['Distance'] = dist

###### Calculate sensor efficiencies ############################

    def add_sensor_efficiencies(self):
        eff = []
        rnd_eff = []
        rnd_ineff = []
        for i in range(self.dim):
            eff_i = float((1. - (self.df.loc[i, 'MissingHits'] / self.df.loc[i, 'TotalHits'])) * 100.)
            eff.append(eff_i)
            rate_i = self.df.loc[i, 'HitRate']
            if np.isnan(rate_i):
                rnd_eff.append((rate_i))
                rnd_ineff.append((rate_i))
            else:
                rnd_eff.append(self.rnd_eff(rate_i))
                rnd_ineff.append(100.0 - self.rnd_eff(rate_i))
        self.df['Efficiency'] = eff
        self.df['RandomEfficiency'] = rnd_eff
        self.df['RandomInefficiency'] = rnd_ineff

        print('Max ineff = ', max(rnd_ineff))

###### Add sensor efficiencies caluclated from random data ######

    def import_rnd_eff(self):
        if self.newfsm == False:
            rnd_file = '../data/rates/rates_rnd_C4785x17600um2_P165x55um2_S739_E1to500000_L1_Q1_wosecs.csv'
        else:
            rnd_file = '../data/rates/rates_rnd_C4785x17600um2_P165x55um2_S739_E1to500000_L1_Q1_wosecs_newfsm.csv'
        self.df_rnd_eff = pd.DataFrame(pd.read_csv(rnd_file)) #old: '../data/efficiencies/rnd_efficiencies_C4785x17600um2_P165x55um2_S739_E1to4000_L1_Q1_wosecs_rep1to10.csv'
        self.df_rnd_eff = self.df_rnd_eff.set_index(['Rate'])

    def rnd_eff(self, rate):
        rounded_rate = float(int(math.ceil(rate)))
        rnd_tot = float(self.df_rnd_eff.loc[rounded_rate, 'TotalHits'])
        rnd_mis = float(self.df_rnd_eff.loc[rounded_rate, 'MissingHits'])
        rnd_eff = 100.* (rnd_tot - rnd_mis)/rnd_tot
        return rnd_eff
    
###### Define mapping function ##################################

    def empty_map(self, y_dim, x_dim):
        map = np.empty((y_dim+1, x_dim+1,))
        map[:] = np.nan
        return map


    def heatmap(self, map, name, label, bar_max = '', bar_min = ''):

        fig, ax = plt.subplots(figsize=(10,6),dpi=100)

        # masked_map = np.ma.masked_where(map == 0, map)
        # cmap = matplotlib.cm.viridis
        # cmap.set_bad(color='white')
        # plt.imshow(masked_map, cmap=cmap, interpolation='nearest', aspect='auto', origin='lower')

        if bar_max != '' and bar_min != '':
            plt.imshow(map, cmap='viridis_r', interpolation='nearest', vmax = bar_max, vmin = bar_min, aspect='auto', origin='lower')
        else: 
            plt.imshow(map, cmap='viridis_r', interpolation='nearest', aspect='auto', origin='lower')

        cbar = plt.colorbar()
        cbar.set_label(label)
        
        plt.xlabel('Sensor x position (mm)')
        plt.ylabel('Sensor y position (mm)')

        x_pos = [0, 20, 40, 60, 80]
        x_labels = []
        for pos in x_pos:
            x_labels.append(str(int(self.df[self.df['XPos'] == pos].X0.iloc[0])))
        ax.set_xticks(x_pos)
        ax.set_xticklabels(x_labels)

        y_pos = [0, 5, 10, 15, 20, 25]
        y_labels = []
        for pos in y_pos:
            y_labels.append(str(int(self.df[self.df['YPos'] == pos].Y0.iloc[0])))
        ax.set_yticks(y_pos)
        ax.set_yticklabels(y_labels)
        
        plt.tight_layout()
        if self.save_plots: plt.savefig('../plots/sensor_maps_thesis/viridis_r/sensor_'+name+'_map_L'+self.l+'_Q'+self.q+'.pdf') 
        if self.show_plots: plt.show()

        plt.clf()
        plt.close('all')
        try: fig.clf()
        except: pass


    def create_map(self, col_name, map_name, map_ax_label, fill_holes = True, bar_max = '', bar_min = ''):

        map = self.empty_map(self.y_dim, self.x_dim)

        for i in range(len(self.df)):

            x_i = int(self.df.loc[i, 'XPos'])
            y_i = int(self.df.loc[i, 'YPos'])

            var_i = self.df.loc[i, col_name]
            if fill_holes:
                if not var_i > 0: var_i = 0

            map[y_i][x_i] = var_i

        self.heatmap(map, map_name, map_ax_label, bar_max, bar_min)


###### Generate individual maps #################################

    def distance_map(self):
        self.create_map('Distance', 'distance', self.distance_label, fill_holes = False)

    def hit_map(self):
        self.create_map('TotalHits', 'hit', self.hits_label, fill_holes = False)

    def hitrate_map(self):
        self.create_map('HitRate', 'hitrate', self.hitrate_label, fill_holes=False, bar_max = 21, bar_min = 0)

    def efficiency_map(self):
        self.create_map('Efficiency', 'efficiency', self.efficiency_label, fill_holes = False, bar_max=100, bar_min = 90)

    def rnd_efficiency_map(self):
        if self.newfsm == False:
            rnd_label = 'rnd_efficiency'
        else:
            rnd_label = 'rnd_newfsm_efficiency'
        self.create_map('RandomEfficiency', rnd_label, self.efficiency_label, fill_holes=False, bar_max = 100, bar_min = 99.47)

    def rnd_inefficiency_map(self):
        if self.newfsm == False:
            rnd_label = 'rnd_inefficiency'
        else:
            rnd_label = 'rnd_newfsm_inefficiency'
        self.create_map('RandomInefficiency', rnd_label, self.inefficiency_label, fill_holes=False, bar_max = 0.52, bar_min = 0.0) #


###### Generate individual plots ################################

    def scatter_plot(self, name, x, y, x_label, y_label, log = False):

        plt.scatter(x, y, s = 4)
        plt.xlabel(x_label)
        plt.ylabel(y_label)

        if log:
            plt.yscale('log')
            name = name + '_log'

        plt.tight_layout()
        if self.save_plots: plt.savefig('../plots/sensor_plots/sensor_'+name+'_plot_L'+self.l+'_Q'+self.q+'.pdf')
        if self.show_plots: plt.show()

        plt.clf()
        plt.close('all')
        try: fig.clf()
        except: pass

    def position_vs_rate(self, log = False):
        self.scatter_plot('distance_vs_hitrate', self.df.Distance, self.df.HitRate, self.distance_label, self.hitrate_label, log)

    def rate_vs_efficiency(self, log = False):
        self.scatter_plot('hitrate_vs_efficiency', self.df.HitRate, self.df.Efficiency, self.hitrate_label, self.efficiency_label, log)

    def distance_vs_efficiency(self, log = False):
        self.scatter_plot('distance_vs_efficiency', self.df.Distance, self.df.Efficiency, self.distance_label, self.efficiency_label, log)

###### Run multiple functions together ##########################

    def create_maps(self):
        matplotlib.rcParams.update({'font.size': 18})
        #self.distance_map()
        self.hit_map()
        self.hitrate_map()
        self.efficiency_map()
        self.rnd_efficiency_map()
        self.rnd_inefficiency_map()

    def create_plots(self, log = False):
        matplotlib.rcParams.update({'font.size': 14})
        self.position_vs_rate(log)
        self.rate_vs_efficiency(log)
        self.distance_vs_efficiency(log)

    def create_all(self):
        self.create_maps()
        self.create_plots()

###### Run ######################################################

layers = [6]
quadrants = [1]

for l in tqdm(layers, desc = 'Layers'):
    for q in tqdm(quadrants, desc = 'Quadrants'):
        plots = SensorPlots(layer = l, quadrant = q, newfsm = True, save_plots = True, show_plots = True)
        matplotlib.rcParams.update({'font.size': 18})
        plots.rnd_inefficiency_map()