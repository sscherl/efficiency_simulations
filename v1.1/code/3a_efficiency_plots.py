#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Plots for random (MC) data testing of pixel matrix model
#
#	Author: Sigrid Scherl
#
#	Created: February 2023
#

import csv
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys
import yaml
from tqdm import tqdm
from statistics import mean, stdev
from scipy.stats import beta

from classes.toolkit import File, LHCParameters
from classes.simulation_data_plot import MyColours

###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Class ####################################################

class RandomPlots:

    def __init__(self, type, events, new_fsm = False, layer = 1, quadrant = 1, min_x = '', max_x = '', min_y = '', max_y = '', log = False, save_plots = True, show_plots = True):
        
        self.type       = type
        self.events     = events
        self.layer      = layer
        self.quadrant   = quadrant
        self.min_x      = min_x
        self.max_x      = max_x
        self.min_y      = min_y
        self.max_y      = max_y
        self.log        = log
        self.save_plots = save_plots
        self.show_plots = show_plots
        self.new_fsm    = new_fsm

        # Plot labels
        self.hitrate_label      = 'Hit Rate (MHz/cm²)'
        self.efficiency_label   = 'Efficiency (%)'
        self.distance_label     = 'Distance to IP (mm)'
        self.hits_label         = 'Total Hits'

        self.plot_count = 0
        self.plot_marker = ['o', 's', '^', 'd']
        self.plot_col = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']

        self.plot_readout_limit = True
        self.plot_line = True
        self.plot_error_band = True
        self.plot_error_bars = True

        # Plot settings
        matplotlib.rcParams.update({'font.size': 14}) # For JINST Proceeding: 19, Thesis rate plots probablz: 14

        self.col = MyColours()
        self.lhc = LHCParameters()

        # Plot rate vs efficiency for random data
        if type == 'random':
            self.rates_file = File(path = config['directories']['data']['rates'], prefix= 'rates_rnd', suffix = '.csv')
            self.rates_file.last_event = self.events
            if self.new_fsm: self.rates_file.extra = 'newfsm'
            self.rates_file.generate_name()
            self.df = pd.DataFrame(pd.read_csv(self.rates_file.whole))
            self.dim = len(self.df)

            self.add_sensor_efficiencies()
            self.add_sensor_corr_efficiencies()
            self.add_sensor_readout_efficiencies()
            self.add_sensor_corr_readout_efficiencies()
            self.add_sensor_corr_tot_efficiencies()
            self.calc_clopperpearson()

        elif type == 'scifi':
            self.eff_file = File(path = config['directories']['data']['rates'], prefix= 'rates_rnd', suffix = '.csv')


###### Calculate hit rates ######################################

    def hit_rate(self, hits):
        time = self.events * self.layer * 25e-9 #s -- time = events * layers * BX
        area = 29*0.0165*320*0.0055 #cm²  -- sensing area = pixel size * rows * cols
        return hits/(time * area * 1e6) #MHz/cm²


###### Calculate sensor efficiencies ############################

    def add_sensor_efficiencies(self):
        eff = []
        for i in range(self.dim):
            det = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits']
            tot = self.df.loc[i, 'TotalHits']
            eff.append(float(det / tot) * 100.)
        self.df['Efficiency'] = eff

    def add_sensor_corr_efficiencies(self):
        eff = []
        for i in range(self.dim):
            det = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] - self.df.loc[i, 'OverMaxROTime']
            tot = self.df.loc[i, 'TotalHits']
            eff.append(float(det / tot) * 100.)
        self.df['CorrEfficiency'] = eff

    def add_sensor_readout_efficiencies(self):
        ro_eff = []
        for i in range(self.dim):
            det = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits']
            tot = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'WithinPrevToT']
            ro_eff.append(float(det / tot) * 100.)
        self.df['ReadoutEfficiency'] = ro_eff

    def add_sensor_corr_readout_efficiencies(self):
        ro_eff = []
        for i in range(self.dim):
            det = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] - self.df.loc[i, 'OverMaxROTime']
            tot = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'WithinPrevToT']
            ro_eff.append(float(det / tot) * 100.)
        self.df['CorrReadoutEfficiency'] = ro_eff

    def add_sensor_corr_tot_efficiencies(self):
        tot_eff = []
        for i in range(self.dim):
            det = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] - self.df.loc[i, 'OverMaxROTime']
            tot = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] + self.df.loc[i, 'WithinPrevToT']
            tot_eff.append(float(det / tot) * 100.)
        self.df['CorrToTEfficiency'] = tot_eff

###### Calculate clopper pearson confidence interval ###############

    def cloppear(self, total, detected):
        k = detected
        n = total
        alpha = 0.05
        low, high = beta.ppf([alpha / 2, 1 - alpha / 2], [k, k + 1], [n - k + 1, n - k])
        
        if np.isnan(low): low = 0
        if np.isnan(high): high = 1
        
        return low, high
    
    def calc_clopperpearson(self):
        lows, highs = [], []
        ro_lows, ro_highs = [], []
        corr_lows, corr_highs = [], []
        corr_ro_lows, corr_ro_highs = [], []
        corr_tot_lows, corr_tot_highs = [], []

        for i in range(self.dim):
            low, high = self.cloppear(total = self.df.loc[i, 'TotalHits'], detected = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'])
            corr_low, corr_high = self.cloppear(total = self.df.loc[i, 'TotalHits'], detected = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] - self.df.loc[i, 'OverMaxROTime'])
            ro_low, ro_high = self.cloppear(total = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'WithinPrevToT'], detected = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'])
            corr_ro_low, corr_ro_high = self.cloppear(total = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'WithinPrevToT'], detected = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] - self.df.loc[i, 'OverMaxROTime'])
            corr_tot_low, corr_tot_high = self.cloppear(total = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] + self.df.loc[i, 'WithinPrevToT'], detected = self.df.loc[i, 'TotalHits'] - self.df.loc[i, 'MissingHits'] - self.df.loc[i, 'OverMaxROTime'])

            lows.append(low*100)
            highs.append(high*100)
            corr_lows.append(corr_low*100)
            corr_highs.append(corr_high*100)
            ro_lows.append(ro_low*100)
            ro_highs.append(ro_high*100)
            corr_ro_lows.append(corr_ro_low*100)
            corr_ro_highs.append(corr_ro_high*100)
            corr_tot_lows.append(corr_tot_low*100)
            corr_tot_highs.append(corr_tot_high*100)

        
        self.df['CPLow'] = lows
        self.df['CPHigh'] = highs

        self.df['Corr_CPLow'] = corr_lows
        self.df['Corr_CPHigh'] = corr_highs

        self.df['RO_CPLow'] = ro_lows
        self.df['RO_CPHigh'] = ro_highs

        self.df['Corr_RO_CPLow'] = corr_ro_lows
        self.df['Corr_RO_CPHigh'] = corr_ro_highs

        self.df['Corr_ToT_CPLow'] = corr_tot_lows
        self.df['Corr_ToT_CPHigh'] = corr_tot_highs


###### Calculate mean of sensor efficiencies ####################

    def calc_mean_efficiencies(self):

        self.mean_df = pd.DataFrame(columns = ['Rate', 'EfficiencyMean', 'EfficiencyStdev'])

        rates = self.df.Rate.to_numpy()
        effs = self.df.Efficiency.to_numpy()

        temp_rates, temp_mean_eff, temp_stdev_eff = [], [], []

        for val in rates:
            temp_list = []

            for i in range(len(rates)):
                if rates[i] == val:
                    temp_list.append(effs[i])

            temp_rates.append(val)
            temp_mean_eff.append(mean(temp_list))
            
            if len(temp_list) > 1:
                temp_stdev_eff.append(stdev(temp_list))
            else:
                temp_stdev_eff.append(0)

        self.mean_df.Rate = temp_rates
        self.mean_df.EfficiencyMean = temp_mean_eff
        self.mean_df.EfficiencyStdev = temp_stdev_eff


###### Define plotting functions ################################

    def conf_plot(self, subdir, prefix, x, y, x_label, y_label, yerr_low = '', yerr_high = '', last = True, label = ''):

        yerr = np.array(list(zip(y - yerr_low, yerr_high - y))).T

        if self.plot_line:
            lw = 1
        else:
            lw = 0

        if self.plot_error_bars:
            plt.errorbar(x, y, yerr= yerr, color = self.plot_col[self.plot_count], label = label, marker = self.plot_marker[self.plot_count], markersize = 4, linewidth=lw, elinewidth = 1)
        else:
            plt.plot(x, y, color = self.plot_col[self.plot_count], label = label, marker = self.plot_marker[self.plot_count], markersize = 4, linewidth=lw)

        if self.plot_error_band:
            plt.fill_between(x, yerr_low, yerr_high, alpha = 0.2, color = self.plot_col[self.plot_count])

        plt.xlabel(x_label)
        plt.ylabel(y_label)

        if self.plot_readout_limit:
            if self.new_fsm:
                limit = self.lhc.readout_limit_mp2
                col = 'tab:orange'#'tab:blue'
            else:
                limit = self.lhc.readout_limit_mp1
                col = 'tab:blue'
            plt.plot([limit, limit],[-10, 110], linewidth = 1.5, linestyle = '--', color = col)

        if self.plot_count == 0:
            self.plot_name = prefix
        else:
            self.plot_name += ('_'+prefix)

        self.plot_count += 1

        if last:

            if self.min_x == '': min_x = x.min() - 0.1
            else: min_x = self.min_x
            if self.max_x == '': max_x = x.max() + 0.1
            else: max_x = self.max_x
            plt.xlim(min_x, max_x)

            if self.min_y == '': min_y = y.min() - 0.1
            else: min_y = self.min_y
            if self.max_y == '': max_y = y.max() + 0.1
            else: max_y = self.max_y
            plt.ylim(min_y, max_y)

            plt.grid(color = 'lightgrey', linestyle = '--', linewidth = 1)

            plt.tight_layout()

            # Plot name
            self.plot_file = File(path = config['directories']['plots'][subdir], prefix = 'plot_rnd_rate_vs_' + self.plot_name)
            self.plot_file.suffix = '.pdf'
            self.plot_file.last_event = self.events
            
            if self.log:
                plt.yscale('log')
                self.plot_file.suffix = '_log' + self.plot_file.suffix 
            if self.min_x != '' or self.max_x != '':
                self.plot_file.suffix = '_zoom' + str(self.max_x) + self.plot_file.suffix
            if self.new_fsm: self.plot_file.extra = 'newfsm'
            self.plot_file.generate_name()

            if label != '':
                plt.legend(loc= 'lower left')

            if self.save_plots: plt.savefig(self.plot_file.whole)
            if self.show_plots: plt.show()

            plt.clf()
            plt.close('all')
            try: fig.clf()
            except: pass


###### Save data to file ########################################

    def save_data(self):
        data_file = File(path = config['directories']['data']['efficiencies'], prefix = 'rnd_efficiencies')
        data_file.suffix = '.csv'
        data_file.last_event = self.events
        if self.new_fsm: data_file.extra = 'newfsm'

        data_file.generate_name()
        self.df.to_csv(data_file.whole)


###### Generate individual plots ################################

    def rate_vs_efficiency(self, last, label = 'Efficiency'):

        self.conf_plot(subdir = 'rate_vs_eff',
                       prefix = 'eff',
                       x = self.df.Rate,
                       y = self.df.Efficiency,
                       x_label = self.hitrate_label,
                       y_label = self.efficiency_label,
                       yerr_low = self.df.CPLow,
                       yerr_high = self.df.CPHigh,
                       last = last,
                       label = label
                    )
        
    def rate_vs_corr_efficiency(self, last, label = 'Corrected Efficiency'):

        self.conf_plot(subdir = 'rate_vs_eff',
                       prefix = 'corr_eff',
                       x = self.df.Rate,
                       y = self.df.CorrEfficiency,
                       x_label = self.hitrate_label,
                       y_label = self.efficiency_label,
                       yerr_low = self.df.Corr_CPLow,
                       yerr_high = self.df.Corr_CPHigh,
                       last = last,
                       label = label
                    )

    def rate_vs_readout_efficiency(self, last, label = 'Readout Efficiency'):

        self.conf_plot(subdir = 'rate_vs_eff',
                       prefix = 'ro_eff',
                       x = self.df.Rate,
                       y = self.df.ReadoutEfficiency,
                       x_label = self.hitrate_label,
                       y_label = self.efficiency_label,
                       yerr_low = self.df.RO_CPLow,
                       yerr_high = self.df.RO_CPHigh,
                       last = last,
                       label = label
                    )
        
    def rate_vs_corr_readout_efficiency(self, last, label = 'Corrected Readout Efficiency'):

        self.conf_plot(subdir = 'rate_vs_eff',
                       prefix = 'corr_ro_eff',
                       x = self.df.Rate,
                       y = self.df.CorrReadoutEfficiency,
                       x_label = self.hitrate_label,
                       y_label = self.efficiency_label,
                       yerr_low = self.df.Corr_RO_CPLow,
                       yerr_high = self.df.Corr_RO_CPHigh,
                       last = last,
                       label = label
                    )

    def rate_vs_corr_tot_efficiency(self, last, label = 'Corrected ToT Efficiency'):

        self.conf_plot(subdir = 'rate_vs_eff',
                       prefix = 'corr_tot_eff',
                       x = self.df.Rate,
                       y = self.df.CorrToTEfficiency,
                       x_label = self.hitrate_label,
                       y_label = self.efficiency_label,
                       yerr_low = self.df.Corr_ToT_CPLow,
                       yerr_high = self.df.Corr_ToT_CPHigh,
                       last = last,
                       label = label
                    )

###### Run multiple functions together ##########################

    def create_plots(self):
        self.rate_vs_efficiency()
        self.rate_vs_corr_efficiency()
        self.rate_vs_readout_efficiency()
        self.rate_vs_corr_readout_efficiency()
        self.rate_vs_corr_tot_efficiency()

    def create_all(self):
        self.create_plots()

###### Run ######################################################

save = True
show = True
min_x = 0.9 # For zoom: 0.9, Full: 0.9
max_x = 35 # For zoom: 25, Full: 40; NewFSM: Zoom: 35
min_y = 99.9 # For zoom: 99, Full: 50
max_y = 100.001 # For zoom: 100.01, Full: 101

plots1 = RandomPlots(type = 'random', events = 500000, new_fsm = False, save_plots = save, show_plots = show, min_x = min_x, max_x = max_x, min_y = min_y, max_y = max_y)
plots2 = RandomPlots(type = 'random', events = 500000, new_fsm = True, save_plots = save, show_plots = show, min_x = min_x, max_x = max_x, min_y = min_y, max_y = max_y)

# plots1.plot_error_band = True
# plots1.plot_line = True
plots1.plot_error_bars = False
plots2.plot_error_bars = False

#plots.rate_vs_corr_efficiency(last = True, label = '') # Readout Efficiency
#plots.rate_vs_corr_readout_efficiency(last = True, label = 'FSM Efficiency')
#plots.rate_vs_corr_tot_efficiency(last = True, label = 'Readout Efficiency (ToT)')

#plots2.plot_name = plots.plot_name
plots1.plot_count = 0
plots1.rate_vs_corr_readout_efficiency(last = False, label = 'MightyPix1')

plots2.plot_name = plots1.plot_name
plots2.plot_count = 1
plots2.rate_vs_corr_readout_efficiency(last = True, label = 'MightyPix2')

# plots2.plot_name = plots.plot_name
# plots2.plot_count = 3
# plots2.rate_vs_corr_readout_efficiency(last = True, label = 'Corrected Readout Efficiency MP2')

