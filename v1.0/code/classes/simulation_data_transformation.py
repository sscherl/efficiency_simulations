#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class to output simulation data for chosen inputs in pixel coordinates
#
#	Author: Sigrid Scherl
#
#	Created: November 2021
#

import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from random import randint
import sys

from classes.simulation_data_generator import DataGenerator
from classes.toolkit import Chip, File, LHCParameters, Helper
from classes.simulation_data_plot import Plot

###### Classes ##################################################

class DataTransformation:

	###### Initialisation ##########################################

	def __init__(self, data, chip, events = '', layers = '', quadrants = '', rate = ''):
		self.lhc = LHCParameters()
		self.helper = Helper()

		self.data_type = data
		self.chip = chip
		self.full_chip_rate = rate
		self.output_path = []
		self.saved_data_files = []

		# Physical simulation data from Zurich group, based on old SciFi tracker geometry
		if self.data_type == 'scifi':
			self.data_file = '../data/raw_data/MagUp_bs2phiphi_violaine_energytime_oscar_files_1p5e34.root_Sensor_146_12_20_20_xySym.txt'

		# Random data from DataGenerator class
		elif self.data_type == 'random':
			self.random_data = DataGenerator(chip, events, rate)

			# For random data layer and quadrant should be set to 1
			if layers == '': layers = 1
			elif layers == 1: pass
			else: self.helper.raise_fail('For random data "layers" has to be set to 1!')

			if quadrants == '': quadrants = 1
			elif quadrants == 1: pass
			else: self.helper.raise_fail('For random data "quadrants" has to be set to 1!')

		else: self.helper.raise_fail('Invalid data_type "{}"! Cannot compute. Choose "scifi" or "random".'.format(self.data_type))

		if events != '': self.set_events(1, events)
		if layers != '': self.set_layers(1, layers)
		if quadrants != '': self.set_quadrants(1, quadrants)

		self.layer_mode = 'appended'                                 # Choose if different layers will be appended or superimposed: ['appended', 'overlapped']
		self.quadrant_mode = 'overlapped'                              # Choose if different quadrants will be appended or superimposed: ['appended', 'overlapped']

		self.data_width = 20 #mm
		self.data_length = 20 #mm
		self.data_size = self.data_width * self.data_length #mm^2

		self.position = 'bottom_left'                                  # Choose subregion where chip is placed ['bottom_left', 'bottom_right', 'top_left', 'top_right']
		self.orientation = 'original'                                  # Choose orientation of data ['original', 'x_mirrored', 'y_mirrored', 'xy_mirrored'] correspond to theoretical orientation of quadrants 1, 4, 2, 3

		self.mev_to_us = 1                                             # Conversion factor of energy in MeV to ToT in us
		self.fixed_tot = 2 #us

		# Secondaries correction
		self.keep_secondaries = True

		# Clustering
		self.clusters = False
		self.cluster_percentage = 0.05

		# Time of flight distribution
		self.tof_plot_percentage = 100
		self.tof_plot_bins = 200

		# Plotting hitmap
		self.min_hits = 0
		self.max_hits = 0
		self.cutoff = False

	###### Set Values ##############################################

	def add_output_path(self, new_path):
		self.output_path.append(new_path)

	def set_events(self, first_event, last_event):                  # Choose event(s) [1 to 500]
		self.events = list(range(first_event,last_event+1))
		self.first_event = first_event
		self.last_event = last_event

	def set_layers(self, first_layer, last_layer):                  # Choose layer(s) [1 to 6]: input only valid for original data_length
		self.layers = list(range(first_layer, last_layer+1))
		self.first_layer = first_layer
		self.last_layer = last_layer

	def set_quadrants(self, first_quadrant, last_quadrant):	        # Choose quadrant(s) [1 to 4]: input only valid for single data_rate
		self.quadrants = list(range(first_quadrant, last_quadrant+1))
		self.first_quadrant = first_quadrant
		self.last_quadrant = last_quadrant

	###### Generate base name ######################################

	def generate_name(self):
		self.file = File()
		if self.data_type == 'scifi':
			self.file.generate_name('', self.chip, self.first_event, self.last_event, self.first_layer, self.last_layer, self.layer_mode, self.first_quadrant, self.last_quadrant, self.quadrant_mode, '', self.orientation, self.keep_secondaries, self.clusters)
		elif self.data_type == 'random':
			self.file.generate_name('rnd', self.chip, self.first_event, self.last_event, self.first_layer, self.last_layer, self.layer_mode, self.first_quadrant, self.last_quadrant, self.quadrant_mode, self.random_data.full_chip_rate, self.orientation, self.keep_secondaries, self.clusters)

	###### Read in data #############################################

	def get_dataframe(self):

		if self.data_type == 'scifi':
			data = pd.read_csv(self.data_file, sep='\t', header= None, comment='#')
			self.df = pd.DataFrame(data)
			self.df.columns = ['Event', 'Layer', 'Quadrant', 'X', 'Y', 'ToF', 'Energy']

		elif self.data_type == 'random':
			self.random_data.full_chip_rate = self.full_chip_rate
			self.random_data.generate_data()
			self.df = self.random_data.get_data()

	###### Transform coordinates ####################################

	def transform_coordinates(self):

		self.df['Col'] = (self.df['X'].add(10)/self.chip.col_width).astype(int)
		self.df['Row'] = (self.df['Y'].add(10)/self.chip.row_width).astype(int)
		if self.mev_to_us == 0:
			self.df['ToT'] = self.fixed_tot
		else:
			self.df['ToT'] = self.df['Energy']*self.mev_to_us

	###### Chip orientation #########################################

	def set_orientation(self):

		if self.orientation == 'x_mirrored' or self.orientation == 'xy_mirrored':
			self.df['Row'] = self.df['Row'].max() - 1 - self.df['Row']

		if self.orientation == 'y_mirrored' or self.orientation == 'xy_mirrored':
			self.df['Col'] = self.df['Col'].max() - 1 - self.df['Col']

	###### Chip size and position ###################################

	def set_size_and_position(self):

		min_col, min_row = 0, 0
		max_col = int(self.data_width/self.chip.col_width)
		max_row = int(self.data_length/self.chip.row_width)
		max_col_new = int(self.chip.width/self.chip.col_width)
		max_row_new = int(self.chip.length/self.chip.row_width)

		if not ((self.data_width == self.chip.width and self.data_length == self.chip.length) or self.chip.width == 0 or self.chip.length == 0):
			if self.position == 'bottom_left':
				max_col = max_col_new
				max_row = max_row_new
			elif self.position == 'bottom_right':
				min_col = max_col - max_col_new
				max_row = max_row_new
			elif self.position == 'top_left':
				max_col = max_col_new
				min_row = max_row - max_row_new
			elif self.position == 'top_right':
				min_col = max_col - max_col_new
				min_row = max_row - max_row_new
			else:
				self.helper.raise_fail('Position for reduced chip size not valid. Options are: "bottom_left", "bottom_right", "top_left", "top_right".')

		self.cols = list(range(min_col, max_col))
		self.rows = list(range(min_row, max_row))

	###### Output data ##############################################

	def append_or_overlap(self):

		self.df = self.df[self.df['Event'].isin(self.events) & self.df['Layer'].isin(self.layers) & self.df['Quadrant'].isin(self.quadrants) & self.df['Col'].isin(self.cols) & self.df['Row'].isin(self.rows)]

		# If more than one quandrant is chosen, they will be appended
		nr_events = self.last_event - self.first_event + 1
		nr_quadrants = self.last_quadrant - self.first_quadrant + 1

		if self.quadrant_mode == 'overlapped' and self.layer_mode == 'overlapped':
			pass
		elif self.quadrant_mode == 'appended' and self.layer_mode == 'overlapped':
			self.df['Event'] = self.df['Event'] + nr_events * (self.df['Quadrant'] - self.first_quadrant)
		elif self.quadrant_mode == 'overlapped' and self.layer_mode == 'appended':
			self.df['Event'] = self.df['Event'] + nr_events * (self.df['Layer'] - self.first_layer)
		elif self.quadrant_mode == 'appended' and self.layer_mode == 'appended':
			self.df['Event'] = self.df['Event'] + nr_events * (self.df['Quadrant'] - self.first_quadrant) + nr_events * nr_quadrants * (self.df['Layer'] - self.first_layer)
		else:
			self.helper.raise_fail('Input for quadrant_mode and/or layer_mode invalid.')

	###### Omit secondaries #########################################

	def omit_secondaries(self):

		# Hits with ToF > (bx_period + min_tof) ns will be deleted, as they are secondaries that MightyPix cannot correctly timestamp anyway
		self.df = self.df.drop(self.df[self.df.ToF > (self.lhc.bx_period_ns+self.lhc.min_tof_theory)].index)
		self.df = self.df.reset_index().drop(columns=['index'])

	###### Calculate time between hits ##############################

	def set_time_between_hits(self):

		# Time of Hit (ToH) = time when hit occured relative to time of first event (T = 0)
		self.df['ToH'] = ((self.df['Event'] - 1) * self.lhc.bx_period_ns) + self.df['ToF']
		self.df = self.df.sort_values(by=['ToH'])
		self.df = self.df.reset_index()

		# Twait = time between hits
		self.df.loc[0, 'Twait'] = self.df.loc[0, 'ToH']
		for i in range(1, len(self.df)):
			self.df.loc[i, 'Twait'] = round(self.df.loc[i, 'ToH'] - self.df.loc[i-1, 'ToH'], 4)

		# Check that wait times are non-negative
		if (self.df['Twait'] < 0).any():
			self.helper.raise_fail('Negative wait times not possible!')

		# Get rid of old index column
		self.df = self.df.drop(columns=['index'])

	###### Generate clusters #######################################

	def generate_clusters(self):

		# Randomly write 5% of the hits from self.df into new dataframe
		nr_of_hits = len(self.df)
		pct_of_hits = int(nr_of_hits * self.cluster_percentage)

		rnd_indices = [randint(0, nr_of_hits) for i in range(0, pct_of_hits)]

		df_with_index = self.df.reset_index()
		df_cluster_hits = df_with_index[df_with_index['index'].isin(rnd_indices)]

		# Change the row and/or column of these by 1
		offset = [-1,1]

		for i in range(len(df_cluster_hits)):
			rnd1 = randint(1,2) # set to randint(1,3) to also include diagonal, not just up/down or left/right

			if rnd1 != 1:
				if df_cluster_hits.Col.iloc[i] == 0:
					df_cluster_hits.Col.iloc[i] += 1
				elif df_cluster_hits.Col.iloc[i] == self.chip.cols-1:
					df_cluster_hits.Col.iloc[i] -= 1
				else:
					df_cluster_hits.Col.iloc[i] += offset[randint(0, 1)]

			if rnd1 != 2:
				if df_cluster_hits.Row.iloc[i] == 0:
					df_cluster_hits.Row.iloc[i] += 1
				elif df_cluster_hits.Row.iloc[i] == self.chip.rows-1:
					df_cluster_hits.Row.iloc[i] -= 1
				else:
					df_cluster_hits.Row.iloc[i] += offset[randint(0, 1)]

		df_cluster_hits.ToT = self.fixed_tot
		df_cluster_hits.Twait = 0

		df_with_index['Sorting'] = 0
		df_cluster_hits['Sorting'] = 1

		# Add new dataframe to old one
		df_cluster = pd.concat([df_with_index, df_cluster_hits]).sort_values(by=['index', 'Sorting']).drop(columns=['Sorting', 'index'])
		self.df = df_cluster.reset_index(drop=True)

	###### Get data rate ###########################################

	def get_data_rate(self):

		# Calculate hits per event
		df_count_events = self.df['Event'].value_counts().to_frame().reset_index()
		df_count_events.columns = ['Event', 'Count']

		# Omit contributions
		# df_count_events.drop(df_count_events.index[0], axis=0, inplace=True)

		self.data_rate = df_count_events.Count.sum()/math.ceil(self.df['Event'].iloc[-1])

	###### Terminal info ############################################

	def print_info(self):

		self.helper.hline()
		print('Data transformation successful!')
		self.helper.hline()

		if self.data_type == 'scifi': print('Data used:\tSciFi')
		elif self.data_type == 'random': print('Data used:\tRandom')
		else: self.helper.raise_fail('Invalid data_type "{}"! Cannot finish. Choose "scifi" or "random".'.format(self.data_type))

		print('Chip size:\t{} mm x {} mm'.format(self.chip.width, self.chip.length))
		print('Pixel size:\t{:2.0f} um x {:2.0f} um'.format(self.chip.col_width*1000, self.chip.row_width*1000))
		print('Columns:\t{}\nRows:\t\t{}'.format(self.chip.cols, self.chip.rows))

		if self.first_layer == self.last_layer:
			print('Layer:\t\t{}'.format(self.first_layer))
		else:
			print('Layers:\t\t{} to {} ({})'.format(self.first_layer, self.last_layer, self.layer_mode))

		if self.first_quadrant == self.last_quadrant:
			print('Quadrant:\t{}'.format(self.first_quadrant))
		else:
			print('Quadrants:\t{} to {} ({})'.format(self.first_quadrant, self.last_quadrant, self.quadrant_mode))

		print('Data rate:\t{:2.2f} Hits per Event and Chip'.format(self.data_rate))

		if self.clusters == True: print('Fake clusters:\tYes')
		else: print('Fake clusters:\tNo')

		self.helper.hline()

	###### Plot hits ###############################################

	def hitmap_hits(self, save = True, show = True):

		hits_file = File('../plots/hits/', 'hitmap_hits', self.file.name_wos, '.pdf')
		hits_plot = Plot(hits_file, save, show)

		hits_plot.set_pixels(self.chip.cols, self.chip.rows)
		hits_plot.set_hits(self.min_hits, self.max_hits, self.cutoff)

		hits_plot.hitmap(self.df)

	###### Plot Time of Flight distribution ########################

	def plot_tof_distribution(self, save = True, show = True):

		all_tofs = self.df.ToF.to_numpy()
		all_tofs.sort()

		percentage = self.tof_plot_percentage
		removals = int(len(all_tofs)*(100-percentage)/100)

		if removals > 0:
			tofs = all_tofs[0:-removals]
		else:
			tofs = all_tofs
			percentage = 100

		tof_name = 'plot_tof'

		if percentage != 100:
			tof_name += '_lower' + str(percentage) + 'pct'

		tof_file = File('../plots/tof/', tof_name, self.file.name_wos, '.pdf')
		tof_plot = Plot(tof_file, save, show)

		tof_plot.set_x(tofs, 'Time of Flight (ns)')
		tof_plot.y_label = 'Count'
		tof_plot.bins = self.tof_plot_bins

		tof_plot.histo()

	###### Save output data #########################################

	def save_data(self):

		for paths in self.output_path:
			saved_data_file = File(paths, 'data', self.file.name_wos, '.csv')
			self.df.to_csv(saved_data_file.whole, index = False)

		print('Data saved as:\n\t{}'.format(saved_data_file.name))
		print('Data saved to:')
		for paths in self.output_path:
			print('\t'+paths)
		self.helper.hline()
		self.saved_data_files.append(str(saved_data_file.name))

	###### Full transformation encompassing all functions ###########

	# Run transformation
	def transform(self, save_me):
		self.generate_name()
		self.get_dataframe()
		self.transform_coordinates()
		self.set_orientation()
		self.set_size_and_position()
		self.append_or_overlap()
		if not self.keep_secondaries: self.omit_secondaries()
		self.set_time_between_hits()
		if self.clusters: self.generate_clusters()
		self.get_data_rate()
		self.print_info()
		if save_me: self.save_data()

	# Run all plots and save them
	def plots(self, save_me, show_me):
		self.hitmap_hits(save = save_me, show = show_me)
		self.plot_tof_distribution(save = save_me, show = show_me)
