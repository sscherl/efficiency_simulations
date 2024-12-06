#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class for helper functions to print and log
#
#	Author: Sigrid Scherl
#
#	Created: March 2022
#

import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from random import randint
import sys
import yaml
from tqdm import tqdm
from statistics import mean, stdev
from scipy.stats import beta


###### Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)

###### Classes ##################################################

class Chip:

    def __init__(self):
        self.name       = config['chip']['name']
        self.cols       = config['chip']['columns']
        self.rows       = config['chip']['rows']
        self.col_width  = config['chip']['column_width'] #mm
        self.row_width  = config['chip']['row_width']    #mm
        self.daisy_chain = config['chip']['daisy_chain']
        


        self.set_parameters()

    def set_parameters(self):
        self.width      = self.cols     * self.col_width
        self.length     = self.rows     * self.row_width
        self.size       = self.width    * self.length
        self.pixels     = self.cols     * self.rows
        self.pixel_size = self.col_width* self.row_width    #mm^2

        if self.daisy_chain:
            self.first      = config['chip']['first']
            self.last       = config['chip']['last']


class Constants:

    def __init__(self):
        self.c = 2.99792458e8 #m/s


class ClopperPearson:

    def __init__(self):
        pass

    def get(self, total, detected):
        k = detected
        n = total
        alpha = 0.05
        low, high = beta.ppf([alpha/2, 1 - alpha/2], [k, k + 1], [n - k + 1, n - k])
        
        if np.isnan(low):   low = 0
        if np.isnan(high):  high = 1
        
        return low*100, high*100

    def result(self, total, detected):
        eff = (detected/total)*100

        low, high = self.get(total, detected)

        print ('='*100)
        print('Result for ', total, ' hits, out of which ', detected, ' are missing.')
        print ('-'*100)
        print('95%-CPI    = ', low, ', ', high)
        print ('-'*100)
        print('Efficiency = ', eff, ' + ', high - eff, ' - ', eff - low)
        print ('='*100)


class File:

    def __init__(self, path = '', prefix = '', base_name = '', suffix = ''):
        
        # Basic inputs
        self.path           = path
        self.prefix         = prefix                                # scifi, rnd, output, data, input, etc.
        self.file_name      = base_name
        self.suffix         = suffix                                # mostly .csv

        # Possible inputs to generate file name
        self.chip           = Chip()                                # Use chip class
        self.first_event    = 1                                     # 1
        self.last_event     = 1                                     # 500 for scifi data
        self.rate           = ''
        self.orientation    = config['orientation']
        self.position       = config['position']
        self.tot            = 2
        self.clusters       = config['clusters']['generate'][0]
        self.fsm_clk        = config['fsm_clk']
        self.repetition     = ''
        self.masked_columns = ''
        self.masked_rows    = ''
        self.extra          = ''


        if config['chip']['mask_columns']:
            if config['chip']['mask_columns_above'] < self.chip.cols:
                self.masked_columns += 'Over'+str(config['chip']['mask_columns_above'])
            if config['chip']['mask_columns_below'] > 0:
                self.masked_columns += 'Below'+str(config['chip']['mask_columns_below'])

        if config['chip']['mask_rows']:
            if config['chip']['mask_rows_above'] < self.chip.rows:
                self.masked_rows += 'Over'+str(config['chip']['mask_rows_above'])
            if config['chip']['mask_rows_below'] > 0:
                self.masked_rows += 'Below'+str(config['chip']['mask_rows_below'])


        if config['data']['type'] == 'scifi':
            self.sensor         = config['sensors']['list']         # Sensor ID (hottest one is 739)
            self.first_layer    = config['layers']['first']
            self.last_layer     = config['layers']['last']
            self.layer_mode     = config['layers']['mode']
            self.first_quadrant = config['quadrants']['first']
            self.last_quadrant  = config['quadrants']['last']
            self.quadrant_mode  = config['quadrants']['mode']
            self.secondaries    = config['secondaries']
            self.rate           = ''

        elif config['data']['type'] == 'random':
            # self.rate           = config['rates']['list']          # use for random data
            self.sensor         = 739
            self.first_layer    = config['layers']['first']
            self.last_layer     = config['layers']['last']
            self.layer_mode     = config['layers']['mode']
            self.first_quadrant = config['quadrants']['first']
            self.last_quadrant  = config['quadrants']['last']
            self.quadrant_mode  = config['quadrants']['mode']
            self.secondaries    = ''

        if self.file_name == '':
            self.generate_name()
        else: self.try_set_parameters()


    def generate_name(self):
        
        #chip_name       = config['chip']['name']
        chip_title      = 'C' + str(int(self.chip.width * 1000))     + 'x' + str(int(self.chip.length * 1000))    + 'um2'
        pixel_title     = 'P' + str(int(self.chip.col_width * 1000)) + 'x' + str(int(self.chip.row_width * 1000)) + 'um2'
        
        file_name = str(chip_title) + '_' + str(pixel_title)

        if self.sensor          != ''           : file_name += '_S'         + str(self.sensor)
        if self.first_event     != ''           : file_name += '_'          + str(self.sub_name('E', self.first_event, self.last_event))
        if self.first_layer     != ''           : file_name += '_'          + str(self.sub_name('L', self.first_layer, self.last_layer, self.layer_mode))
        if self.first_quadrant  != ''           : file_name += '_'          + str(self.sub_name('Q', self.first_quadrant, self.last_quadrant, self.quadrant_mode))
        if self.rate            != ''           : file_name += '_R'         + (str(self.rate).replace('.', '-'))
        if self.orientation     != 'original'   : file_name += '_'          + str(self.orientation)
        if self.position        != 'bottom_left': file_name += '_'          + str(self.position)
        if self.tot             != 2            : file_name += '_ToT'       + str(self.tot) + 'us'
        if not self.secondaries                 : file_name += '_wosecs'
        if self.clusters                        : file_name += '_clusters'
        if self.fsm_clk         != '40MHz'      : file_name += '_'          + str(self.fsm_clk)
        if self.repetition      != ''           : file_name += '_rep'       + str(self.repetition)
        if self.masked_columns  != ''           : file_name += '_maskCols'  + str(self.masked_columns)
        if self.masked_rows     != ''           : file_name += '_maskRows'  + str(self.masked_rows)
        if self.extra           != ''           : file_name += '_'          + str(self.extra)
        self.file_name = file_name

        self.set_parameters()


    def sub_name(self, letter, first, last, mode = ''):
        title = str(letter) + str(first)
        if first != last: title += 'to' + str(last)
        if mode == 'appended' and first != last: title += '_app'
        return title


    def set_parameters(self):
        if self.prefix != '':
            self.name_wos = self.prefix+'_'+self.file_name          # name without suffix, temporarily needed
        else:
            self.name_wos = self.file_name
        self.name = self.name_wos+self.suffix
        self.whole = self.path+self.name


    def try_set_parameters(self):
        if self.path != '' and self.file_name != '' and self.suffix != '': self.set_parameters()


class Helper:

    def __init__(self):
        self.log_file = ''
        self.summary = []

    def raise_fail(self, message):
        print ('ERROR:', message)
        sys.exit()

    def printing(self, text = '', to_log = False):
        print(text)
        if to_log:
            if type(text) == pd.core.frame.DataFrame:
                self.log_file.write(text.to_string()+'\n')
            else:
                self.log_file.write(text+'\n')

    def hline(self, symbol = '-', rep = 1, length = 66, offset = 0, to_log = False):
        for i in range(rep):
            txt = str(offset*'\t'+length*symbol)
            self.printing(txt, to_log)

    def print_debug(self, thing, symbol = '*', rep = 1, offset = 3):
        pd.set_option("display.max_rows", None, "display.max_columns", None)
        self.hline(symbol,rep)
        self.printing(offset*'\t'+'...DEBUGGING...')
        self.hline(symbol,rep, 15, offset)
        self.printing(thing)
        self.hline(symbol,rep)
        pd.reset_option("display.max_rows", "display.max_columns")

    def msg(self, text, type = 'body', to_log = False):
        if not text: self.hline(to_log = to_log)
        else:
            if type == 'head': self.hline(to_log = to_log)
            else: self.printing(to_log = to_log)
            self.printing(text,to_log)
            if type == 'foot': self.hline(to_log = to_log)

    def give_update(self, text, to_log = False):
        self.summary.append(text)
        self.msg(self.summary[-1], to_log = to_log)

    # def generate_list_from_config(name):
    #     if config[name]['list'] != []:
    #         list = config[name]['list']
    #     else:
    #         list = [*range(config[name]['first'], config[name]['last']+1)]
    #     return list
    
    def generate_list(self, first, last):
        return [*range(first, last + 1)]
    
    def add_to_name(self, old, new):
        old += '_' + str(new)

    

class AstroPixParameters:

    def __init__(self):
        self.event_length       = 50e-9    #s
        self.event_length_ns    = 50       #ns

        self.event_frequency    = 1 / (self.event_length)   #Hz

        self.readout_time_limit_s = 6.5536e-3 # s

        self.min_tof_theory = 0
        self.max_tof = 0


class LHCParameters:

    def __init__(self):

        self.event_length           = 25e-9 #s
        self.event_length_us        = self.event_length * 1e6   #us
        self.event_length_ns        = self.event_length * 1e9   #ns

        self.event_frequency        = 1. / (self.event_length)   #Hz
        self.event_frequency_MHz    = self.event_frequency / 1e6    #MHz

        self.max_event              = 3564
        self.max_time               = self.max_event * self.event_length
        self.max_time_us            = self.max_event * self.event_length_us
        self.max_time_ns            = self.max_event * self.event_length_ns


        self.bx_period              = 25e-9     #s
        self.bx_period_ns           = self.bx_period * 1e9        #ns

        self.bx_frequency = 40e6 #Hz
        self.bx_frequency_MHz = 40 #MHz

        # The minimum time of flight of a particle generated at the IP, reaching the Mighty Tracker
        self.min_tof_data = 25
        self.min_tof_theory = 26
        self.max_tof = 31

        # Limits to efficiency and readout time set by experiment specs
        self.readout_limit_mp1      = 23.75 #MHz/cm^2
        self.readout_limit_mp2      = 31.66 #MHz/cm^2
        self.readout_time_limit     = 89.1 #us
        self.readout_time_limit_s   = 89.1e-6 #s


class ProgressBar:
# Source: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters?noredirect=1&lq=1
# Alternative: https://github.com/tqdm/tqdm

    def __init__(self, total):
        self.total = total

    def PrintBar (self, iteration1, prefix = 'Progress', suffix = 'Complete', decimals = 1, length = 50, fill = '█', printEnd = "\r"):
        iteration = iteration1 + 1        
        """
        Call in a loop to create terminal progress bar
        @params:
            iteration   - Required  : current iteration (Int)
            total       - Required  : total iterations (Int)
            prefix      - Optional  : prefix string (Str)
            suffix      - Optional  : suffix string (Str)
            decimals    - Optional  : positive number of decimals in percent complete (Int)
            length      - Optional  : character length of bar (Int)
            fill        - Optional  : bar fill character (Str)
            printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(self.total)))
        filledLength = int(length * iteration // self.total)
        bar = fill * filledLength + '-' * (length - filledLength)
        print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == self.total: print()
