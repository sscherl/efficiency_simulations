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


class PlotIVs:

    def __init__(self, plot_multiple = False, plot_name = 'ivs', readings = 1):

        self.plot_name       = plot_name + '_chip'
        self.measurement_dir = '/Users/sigridscherl/Documents/GitRepos/mightypix_measurements/gecco-daq/output/' #'/Users/sigridscherl/Documents/PhD/Projects/MightyPix1/Measurements/'
        self.data_dir        = 'iv_scan/' #'2024-01_FIBedMP1_InitialTests/FirstIVs/'

        # Ranges for plot axes
        self.current_range   = 'nA' # 'pA', 'nA', 'uA'
        self.voltage_range   = 'V'

        # Plot settings
        self.plot_lines      = False
        self.plot_logy       = False

        # Remove outliers
        self.rem_outliers = False
        self.rem_pos_curr = False

        # Plot errors or not
        self.plot_errorbars  = False
        self.error_source    = '' # 'keithley', 'repeat', 'constant'

        # For Keithley errors
        self.keithley_volt_range = '200V'
        self.keithley_curr_range = '100nA' # '10nA', '100nA', '1uA', '10uA', '100uA'

        # For constant errors
        self.const_volt_error = None
        self.const_curr_error = None

        # For repeat measurements
        self.repeat_readings = readings

        # Font size for plot labels
        matplotlib.rcParams.update({'font.size': 14})

        if plot_multiple: 
            self.marker_counter = 0
        else:
            self.marker_counter = -1

        self.markers = ['v', '^', '<', '>', 'o']

        self.colour_counter = 0
        self.colours = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']



    # Add a new chip to be plottet
    def add_chip(self, data_name, chip_name, extra_label = '', marker = '', col = ''):

        self.data_name = data_name
        self.chip_name = chip_name

        self.data_file = self.measurement_dir + self.data_dir + self.data_name + '_chip' + self.chip_name + '.csv'

        if marker != '': self.marker_counter = marker
        if    col != '': self.colour_counter = col

        self.extra_label = extra_label

        self.plot_name += ('_' + self.chip_name)

        self.run()


    # Get data of chip
    def get_data(self):

        self.df_raw_data = pd.DataFrame(pd.read_csv(self.data_file))

        self.df_raw_data['AbsSetVoltage'] = self.df_raw_data.SetVoltage * -1
        self.df_raw_data['AbsMeasuredCurrent'] = self.df_raw_data.MeasuredCurrent * -1

        
        self.mean_voltages, self.mean_currents, self.std_voltages, self.std_currents = [], [], [], []

        temp_volt = self.df_raw_data.AbsSetVoltage.tolist()
        temp_curr = self.df_raw_data.AbsMeasuredCurrent.tolist()

        while len(temp_volt) > 0:
            list_volt, list_curr = [], []

            for i in range(self.repeat_readings):

                list_volt.append(temp_volt[0])
                list_curr.append(temp_curr[0])

                temp_volt.pop(0)
                temp_curr.pop(0)

            self.mean_voltages.append(np.mean(list_volt))
            self.mean_currents.append(np.mean(list_curr))

            self.std_voltages.append(np.std(list_volt))
            self.std_currents.append(np.std(list_curr))

        self.df_data = pd.DataFrame(columns=['Voltage', 'StdVoltage', 'Current', 'StdCurrent'])
        self.df_data.Voltage = self.mean_voltages
        self.df_data.StdVoltage = self.std_voltages
        self.df_data.Current = self.mean_currents
        self.df_data.StdCurrent = self.std_currents

    # Remove data points below outlier_min and/or above outlier_max
    def remove_outliers(self):

        if self.rem_outliers:

            try:    self.df_data = self.df_data[self.df_data.Current < self.outlier_max]
            except: print('No maximum for values found. No high outliers removed')

            try:    self.df_data = self.df_data[self.df_data.Current > self.outlier_min]
            except: print('No minimum for values found. No low outliers removed')

    def remove_positive_currents(self):
        if self.rem_pos_curr:
            self.df_data = self.df_data[self.df_data.Current > 0]


    def add_errors(self):
        if   self.error_source == 'keithley': self.add_keithley_error()
        elif self.error_source == 'constant': self.add_constant_error()
        elif self.error_source == 'repeat':   self.add_repeat_error()
        else: pass


    # Keithley errors
    def add_keithley_error(self):
        self.df_data['VoltErr'] = self.df_data['Voltage'].apply(self.keithley_volt_err)
        self.df_data['CurrErr'] = self.df_data['Current'].apply(self.keithley_curr_err)

    # Voltage accuracy for the Keithley 2470 SMU (as source)
    def keithley_volt_err(self, val):

        range = self.keithley_volt_range

        if range == '200V': accuracy = (val * 0.00015) + 24e-3 # in V
        else: print('ERROR: Voltage range not implemented!')

        return accuracy

    # Current accuracy for the Keithley 2470 SMU (as measurement)
    def keithley_curr_err(self, val):

        range = self.keithley_curr_range

        if   range ==  '10nA': accuracy = (val * 0.001)   + 250e-12 # in A
        elif range == '100nA': accuracy = (val * 0.0006)  + 300e-12 # in A
        elif range ==   '1uA': accuracy = (val * 0.00025) + 300e-12 # in A
        elif range ==  '10uA': accuracy = (val * 0.00025) + 700e-12 # in A
        elif range == '100uA': accuracy = (val * 0.0002)  + 6e-9 # in A
        else: print('ERROR: Current range not implemented!')

        return accuracy

    # Set error to constant value
    def add_constant_error(self):
        if self.const_volt_error != None: self.df_data['VoltErr'] = [self.error_volt] * len(self.df_data)
        else: print('WARNING: Constant error on voltage set to NONE. No error added.')

        if self.const_curr_error != None: self.df_data['CurrErr'] = [self.error_curr] * len(self.df_data)
        else: print('WARNING: Constant error on current set to NONE. No error added.')

    # Average over repeat readings to get mean and standard deviation as error
    def add_repeat_error(self):

        self.mean_voltages, self.mean_currents, self.std_voltages, self.std_currents = [], [], [], []

        temp_volt = self.df_data.Voltage.tolist()
        temp_curr = self.df_data.Current.tolist()

        while len(temp_volt) > 0:
            list_volt, list_curr = [], []

            for i in range(self.repeat_readings):

                list_volt.append(temp_volt[0])
                list_curr.append(temp_curr[0])

                temp_volt.pop(0)
                temp_curr.pop(0)

            self.mean_voltages.append(np.mean(list_volt))
            self.mean_currents.append(np.mean(list_curr))

            self.std_voltages.append(np.std(list_volt))
            self.std_currents.append(np.std(list_curr))


    def add_iv(self):

        if   self.current_range == 'nA': current_multi = 1e9
        elif self.current_range == 'uA': current_multi = 1e6
        elif self.current_range == 'pA': current_multi = 1e12
        else:
            current_multi = 1
            print('Unit not implemented. Setting current unit to "A".')

        if self.chip_name == 'none':
            label = 'Background'
        else: 
            label = 'Chip '+self.chip_name.replace('_', ' ')

        if self.extra_label != '':
            label += ', '+self.extra_label

        if self.plot_errorbars and self.error_source == 'repeat':
            x = self.mean_voltages
            y = [i * current_multi for i in self.mean_currents]
            xerr = self.std_voltages
            yerr = [i * current_multi for i in self.std_currents]

        else:
            x = self.df_data.Voltage
            y = self.df_data.Current * current_multi

            try:    yerr = self.df_data.CurrErr * current_multi
            except:
                yerr = None

            try:    xerr = self.df_data.VoltErr
            except:
                xerr = None

        if self.plot_lines: linewidth = 1
        else:               linewidth = 0

        plt.errorbar(x, y, xerr = xerr, yerr = yerr,
                     linewidth = linewidth, elinewidth = 1,
                     marker = self.markers[self.marker_counter], markersize = 4,
                     label = label,
                     color = self.colours[self.colour_counter]
                     )
        
        print(max(y))
        # x1 = np.arange(-10, 250)
        # plt.plot(x1, self.slope1*current_multi*x1+self.intercept1*current_multi, linewidth=2)

        # x2 = np.arange(230, 255)
        # plt.plot(x2, self.slope2*current_multi*x2+self.intercept2*current_multi, linewidth=2)

        plt.xlabel('Applied Voltage (- '+self.voltage_range+')')
        if self.current_range == 'uA':
            plt.ylabel('Current (- μA)')
        else:
            plt.ylabel('Current (- '+self.current_range+')')

        if self.plot_logy: plt.yscale('log')

        try: plt.xlim(self.x_min,self.x_max)
        except: pass

        try: plt.ylim(self.y_min,self.y_max)
        except: pass

        self.marker_counter += 1
        self.colour_counter += 1


    # Plot all added IVs
    def plot_ivs(self):

        plt.tight_layout()
        plt.legend(loc = 'upper left')

        if self.plot_logy: self.plot_name += '_log'

        plt.savefig(self.plot_name+'.pdf')
        plt.show()
            

    # Run all functions
    def run(self):

        self.get_data()
        self.remove_outliers()
        self.remove_positive_currents()
        self.add_errors()
        self.add_iv()


        #self.linear_fit()
        #self.breakdown_fit()




    def linear_fit(self):
        slope, intercept = np.polyfit(self.df_data.SetVoltage, self.df_data.MeasuredCurrent, 1)
        print('Linear fit (Chip ', self.chip_name, '): slope = ', slope, ', intercept = ', intercept, '(R = ', round(1/slope,0), 'Ohm)')

    def quadratic_fit(self):
        pass

    def breakdown_fit(self):

        x1 = self.df_data.SetVoltage[self.df_data.SetVoltage < 200]
        y1 = self.df_data.MeasuredCurrent[self.df_data.SetVoltage < 200]

        self.slope1, self.intercept1 = np.polyfit(x1, y1, 1)
        print('Linear fit 1 (Chip ', self.chip_name, '): slope = ', self.slope1, ', intercept = ', self.intercept1)

        x2 = self.df_data.SetVoltage[self.df_data.SetVoltage > 235]
        y2 = self.df_data.MeasuredCurrent[self.df_data.SetVoltage > 235]

        self.slope2, self.intercept2 = np.polyfit(x2, y2, 1)
        print('Linear fit 2(Chip ', self.chip_name, '): slope = ', self.slope2, ', intercept = ', self.intercept2)



###### IVs of all chips ######

# ivs.data_dir = '2024-01_FIBedMP1_InitialTests/FirstIVs/'
# ivs.add_chip(data_name = '2024-01-30-17_05_35', chip_name = '09-D2_IIB', curr_range= '100uA', volt_range = '200V',)
# ivs.add_chip(data_name = '2024-01-31-18_49_57', chip_name = 'FBI2_#4', curr_range= '100nA', volt_range = '200V')
# ivs.add_chip(data_name = '2024-01-31-18_44_31', chip_name = '1', curr_range= '100uA', volt_range = '200V')
# ivs.add_chip(data_name = '2024-01-31-18_31_33', chip_name = '2', curr_range= '100uA', volt_range = '200V')
# ivs.add_chip(data_name = '2024-01-31-18_21_39', chip_name = '3', curr_range= '100uA', volt_range = '200V')
# ivs.add_chip(data_name = '2024-01-31-17_46_03', chip_name = '4', curr_range= '100nA', volt_range = '200V')
# ivs.add_chip(data_name = '2024-01-31-17_56_04', chip_name = '5', curr_range= '100nA', volt_range = '200V')
# ivs.add_chip(data_name = '2024-01-31-18_07_36', chip_name = '6', curr_range= '100nA', volt_range = '200V')


###### IV biased over subchip versus subpix ######

# ivs.data_dir = '2024-02_SUBCHIPvsSUBPIX/'
# ivs.add_chip(data_name = '2024-02-21-15_46_50', chip_name = '2', curr_range= '100uA', volt_range = '200V', extra_label='SUBPIX')
# ivs.add_chip(data_name = '2024-02-21-16_04_51', chip_name = '2', curr_range= '100uA', volt_range = '200V', extra_label='SUBCHIP')


###### IVs with PCB without chip for background ######

# ivs.add_chip(data_name = '2024-02-26-17_31_41', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 3 s', marker=1, col = 'tab:green')
# ivs.add_chip(data_name = '2024-02-26-17_36_13', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 3 s', marker=1, col = 'tab:green')

# ivs.add_chip(data_name = '2024-02-26-16_07_54', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 5 s', marker=0, col = 'tab:blue')
# ivs.add_chip(data_name = '2024-02-26-15_02_23', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 5 s', marker=0, col = 'tab:blue')
# ivs.add_chip(data_name = '2024-02-26-14_54_11', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 5 s', marker=0, col = 'tab:blue')
# ivs.add_chip(data_name = '2024-02-26-17_25_30', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 5 s', marker=0, col = 'tab:blue')

# ivs.add_chip(data_name = '2024-02-26-17_40_27', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 7 s', marker=2, col = 'tab:pink')
# ivs.add_chip(data_name = '2024-02-26-17_47_40', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 7 s', marker=2, col = 'tab:pink')

# ivs.add_chip(data_name = '2024-02-26-17_55_28', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 10 s', marker=3, col = 'tab:orange')
# ivs.add_chip(data_name = '2024-02-26-18_02_44', chip_name = 'none', curr_range= '100nA', volt_range = '200V', extra_label='Settling time = 10 s', marker=3, col = 'tab:red')


###### IVs with PCB without chip for background - 5 measurements per value ######
        
ivs = PlotIVs(plot_multiple = False, plot_name = 'IVCurve', readings = 1)

ivs.plot_errorbars = False
#ivs.error_source = 'keithley'
#ivs.repeat_readings = 5
# ivs.y_max = 5
# ivs.y_min = -0.99
# ivs.x_max = 10
# ivs.x_min = -0.99
ivs.current_range   = 'nA'
# ivs.plot_logy = True

##### Thesis plots ######

### Chip 09-D2 IIB:
# ivs.add_chip(data_name = '2024-01-30-17_05_35', chip_name = '09-D2_IIB')

### Chip FBI2 #4:
# ivs.add_chip(data_name = '2024-01-31-18_49_57', chip_name = 'FBI2_#4')

### Chip 1:
#ivs.add_chip(data_name = '2024-01-31-18_44_31', chip_name = '1')

### Chip 2:
# ivs.add_chip(data_name = '2024-01-31-18_31_33', chip_name = '2')

### Chip 3:
# ivs.add_chip(data_name = '2024-01-31-18_21_39', chip_name = '3')

### Chip 4:
# ivs.add_chip(data_name = '2024-01-31-17_46_03', chip_name = '4')

### Chip 5:
ivs.add_chip(data_name = '2024-01-31-17_56_04', chip_name = '5')

### Chip 6:
#ivs.add_chip(data_name = '2024-01-31-18_07_36', chip_name = '6')

### Chip 2 - subchip vs subpix
#ivs.data_dir = '2024-02_SUBCHIPvsSUBPIX/'
# ivs.add_chip(data_name = '2024-02-21-15_46_50', chip_name = '2', extra_label='SUBPIX')
# ivs.add_chip(data_name = '2024-02-21-16_04_51', chip_name = '2', extra_label='SUBCHIP')

### Background (empty PCB)
#ivs.add_chip(data_name = '2024-02-27-14_30_19', chip_name = 'none')

ivs.plot_ivs()



