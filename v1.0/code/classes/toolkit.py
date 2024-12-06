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

###### Classes ##################################################

class Chip:

    def __init__(self, name = 'MightyPix1', cols = 29, rows = 320, col_width = 0.165, row_width = 0.055):
        self.name = name
        self.cols = cols
        self.rows = rows
        self.col_width = col_width #mm
        self.row_width = row_width #mm

        self.set_parameters()

    def set_parameters(self):
        self.width = self.cols * self.col_width
        self.length = self.rows * self.row_width
        self.size = self.width * self.length
        self.pixels = self.cols * self.rows
        self.pixel_size = self.col_width * self.row_width

class Constants:
    def __init__(self):
        self.c = 2.99792458e8 #m/s


class File:

    def __init__(self, path = '', type_name = '', data_name = '', suffix = ''):
        self.path = path
        self.type_name = type_name
        self.data_name = data_name
        self.suffix = suffix

        self.try_set_parameters()

    def generate_name(self, prefix, chip, first_event = '', last_event = '', first_layer = '', last_layer = '', layer_mode = '', first_quadrant = '', last_quadrant = '', quadrant_mode = '', rate = '', orientation = 'original', keep_secondaries = True, clusters = False):
        chip_title = 'C'+str(int(chip.width*1000))+'x'+str(int(chip.length*1000))+'um2'
        pixel_title = 'P'+str(int(chip.col_width*1000))+'x'+str(int(chip.row_width*1000))+'um2'

        file_name = str(chip_title)+'_'+str(pixel_title)

        if first_event != '': file_name += '_'+str(self.sub_name('E', first_event, last_event))
        if first_layer != '': file_name += '_'+str(self.sub_name('L', first_layer, last_layer, layer_mode))
        if first_quadrant != '': file_name += '_'+str(self.sub_name('Q', first_quadrant, last_quadrant, quadrant_mode))
        if rate != '': file_name += '_R'+(str(rate).replace('.', '-'))
        if orientation != 'original': file_name += '_'+str(orientation)
        if not keep_secondaries: file_name += '_wosecs'
        if clusters: file_name += '_clusters'

        self.type_name = prefix
        self.data_name = file_name

        self.set_parameters()

    def sub_name(self, letter, first, last, mode = ''):
        title = str(letter)+str(first)
        if first != last: title += 'to'+str(last)
        if mode == 'appended' and first != last: title += '_app'
        return title

    def set_parameters(self):
        if self.type_name != '':
            self.name_wos = self.type_name+'_'+self.data_name # name without suffix, temporarily needed
        else:
            self.name_wos = self.data_name
        self.name = self.name_wos+self.suffix
        self.whole = self.path+self.name

    def try_set_parameters(self):
        if self.path != '' and self.data_name != '' and self.suffix != '': self.set_parameters()


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


class LHCParameters:

    def __init__(self):
        self.bx_period = 25e-9 #s
        self.bx_period_ns = 25 #ns

        self.bx_frequency = 40e6 #Hz
        self.bx_frequency_MHz = 40 #MHz

        # The minimum time of flight of a particle generated at the IP, reaching the Mighty Tracker
        self.min_tof_data = 25
        self.min_tof_theory = 26
        self.max_tof = 31


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