#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class to compare the data going into and coming out of MightyPix pixel matrix model
#
#	Author: Sigrid Scherl
#
#	Created: January 2022
#

import csv
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from random import randint
import sys
from os import path
import yaml
from statistics import mean

from classes.simulation_data_plot import Plot
from classes.toolkit import Chip, File, LHCParameters, Helper, ClopperPearson

# pd.set_option("display.max_rows", None, "display.max_columns", None)


########## Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


######### Classes ###################################################

class DataComparison:


	########## Initialisation #######################################

	def __init__(self, data, sensor, events, rate, layers, quadrants, secondaries, clusters, tot, extra, repetition, save_plots):

		self.chip				= Chip()

		self.data_type			= data
		self.sensor				= sensor

		self.first_event 		= 1
		self.last_event			= events

		self.first_layer		= 1
		self.last_layer			= layers
		self.layer_mode			= config['layers']['mode']

		self.first_quadrant		= 1
		self.last_quadrant		= quadrants
		self.quadrant_mode		= config['quadrants']['mode']

		self.rate 				= rate
		self.tot 				= tot
		self.secondaries 		= secondaries
		self.clusters 			= clusters
		self.extra				= extra
		self.repetition			= repetition

		self.orientation 		= config['orientation']
		self.position 			= config['position']

		self.fsm_clk 			= config['fsm_clk']

		self.print_log			= not config['mute']
		self.save_log			= not config['mute']

		self.save_plots			= save_plots
		self.show_plots			= config['plots']['show']
		
		# For logging functions
		self.missing = []

		# For plotting hitmaps
		self.min_hits = 0
		self.max_hits = 0
		self.cutoff = False

		# For secondaries correction
		self.lhc = LHCParameters()
		self.helper = Helper()
		self.ci = ClopperPearson()

		self.max_rotime = self.lhc.max_time_ns

		self.generate_name()


	########## Generate base name ###################################

	def generate_name(self):

		self.file 				= File()
		self.file.chip 			= self.chip
		self.file.first_event 	= self.first_event
		self.file.last_event 	= self.last_event
		self.file.first_layer 	= self.first_layer
		self.file.last_layer	= self.last_layer
		self.file.layer_mode	= self.layer_mode
		self.file.first_quadrant= self.first_quadrant
		self.file.last_quadrant = self.last_quadrant
		self.file.quadrant_mode = self.quadrant_mode
		self.file.orientation 	= self.orientation
		self.file.position 		= self.position
		self.file.tot 			= self.tot
		self.file.secondaries 	= self.secondaries
		self.file.clusters 		= self.clusters
		self.file.fsm_clk 		= self.fsm_clk
		self.file.rate			= ''
		self.file.repetition	= self.repetition
		self.file.extra			= self.extra

		if self.data_type == 'random':
			self.file.prefix = 'rnd'

		self.file.generate_name()
		self.rate_file				= File('../data/rates/',		'rates',		self.file.name,	'.csv')

		self.file.rate				= self.rate
		self.file.sensor			= ''
		self.file.generate_name()

		self.empty_efficiency_file	= File('../data/efficiencies/',	'efficiencies', 'empty',		'.csv')
		self.efficiency_file		= File('../data/efficiencies/',	'efficiencies', self.file.name, '.csv')
			
		self.file.sensor 			= self.sensor
		self.file.generate_name()

		self.input_file				= File('../data/input_data/',	'data',			self.file.name, '.csv')
		self.output_file			= File('../data/output_data/',	'output',		self.file.name, '.csv')
		self.output_extra_file		= File('../data/output_data/',	'output_extra',	self.file.name, '.csv')
		self.log_file				= File('../data/log_files/',	'log',			self.file.name, '.txt')
		self.file_name				= self.file.name_wos


	########## Functions ############################################

	def update_missing(self, new_input, first = False):

		if first:
			self.missing.append(new_input)
		else:
			self.missing.append(-new_input)

		if self.print_log:
			if sum(self.missing) == 0:
				self.helper.msg('>>> All hits accounted for!', to_log = self.save_log)
			elif sum(self.missing) == 1 or self.missing == -1:
				self.helper.msg('>>> {:d} hit not yet accounted for'.format(sum(self.missing)), to_log = self.save_log)
			else:
				self.helper.msg('>>> {:d} hits not yet accounted for'.format(sum(self.missing)), to_log = self.save_log)


	########## Log file #############################################

	def start_log(self):

		self.log = open(self.log_file.whole, 'w')
		self.helper.log_file = self.log


	########## Read data ############################################

	def read_data(self):

		# Read output data (from pixel matrix simulation)
		self.data_out_all = pd.read_csv(self.output_file.whole)
		self.df_out_all = pd.DataFrame(self.data_out_all)
		self.df_out_all.drop(self.df_out_all.index[0], axis=0, inplace=True)
		self.df_out_all = self.df_out_all.sort_values(by=['TS']).reindex()
		self.df_out_all = self.df_out_all.rename(columns={'RO_Time': 'T_out', 'ToT': 'ToT_out'})

		# Read output data containing readin times
		self.data_readin_all = pd.read_csv(self.output_extra_file.whole)
		self.df_readin_all = pd.DataFrame(self.data_readin_all)
		self.df_readin_all.drop(columns=['#','ToT', 'TS', 'Amp'])

		# Read full input data (which is sent into pixel matrix simulation) and add readin times
		self.data_in_all = pd.read_csv(self.input_file.whole)
		self.df_in_all = pd.DataFrame(self.data_in_all)
		self.df_in_all['T_in'] = self.df_readin_all['RO_Time']
		self.df_in_all = self.df_in_all.rename(columns={'Event': 'TS', 'ToT': 'ToT_in'})


	########## TS correction for spill over #########################

	def ts_correction(self):

		# Original timestamps
		self.df_in_all['TS_orig'] = self.df_in_all['TS']

		# TS offsets (for spillover)
		self.df_in_all['TS_offset'] = ((self.df_in_all['ToF'] - self.lhc.min_tof_data)/self.lhc.event_length_ns).astype(int)

		# TS MightyPix is expected to show
		self.df_in_all['TS'] +=  self.df_in_all['TS_offset']


	########## Print header #########################################

	def print_header(self):

		self.helper.msg('ANALYSIS: Hits in the MightyPix Pixelmatrix model', type = 'head', to_log = self.save_log)
		self.helper.msg('FSM Clock speed: 40 MHz\nSensor: '+str(self.sensor), type = 'head', to_log = self.save_log)
		self.helper.msg('Input  file:\t'+self.input_file.name+'\nOutput file:\t'+self.output_file.name, type = 'head', to_log = self.save_log)


	########## Compare in- and outgoing hits ########################

	def check_in_out_hits(self):

		# Select input values to compare
		self.df_in = self.df_in_all[['Col','Row','TS', 'ToT_in', 'T_in']]
		self.df_in = self.df_in.round({'ToT_in': 2})

		# Select output values to compare
		self.df_out = self.df_out_all[['Col', 'Row', 'TS', 'ToT_out', 'T_out']]
		self.df_out['TS'] -= (self.df_out['TS'].min() - self.df_in['TS'].min())
		self.df_out['ToT_out'] = self.df_out['ToT_out'] / 1000
		self.df_out = self.df_out.round({'ToT_out': 2})

		# Merge in- and outgoing hits
		self.df = pd.merge(self.df_in, self.df_out, on = ['Col', 'Row', 'TS'], how = 'outer', indicator = 'Exist')
		self.df['Exist'].replace({'left_only': 'in', 'right_only': 'out'}, inplace = True)
		# self.df = self.df.sort_values(by=['TS']).reindex()

		# T_readout is time to readout the hit in ns
		self.df['T_readout'] = self.df['T_out'] - self.df['T_in']

		# TS_read is the timestamp when the hit was read out (point in time), TS_readout is how many TS it took (duration)
		self.df['TS_readout'] = self.df['T_readout'] / self.lhc.event_length_ns
		self.df['TS_read'] = self.df['TS_readout'] + self.df['TS']

		# Difference between ingoing and outgoing hits
		missing_hits = len(self.df_in) - len(self.df_out)

		# Report results
		if self.print_log:
			self.helper.msg('NUMBER OF HITS...', type = 'head', to_log = self.save_log)
			self.helper.give_update('Incoming hits: {:d}\nOutgoing hits: {:d}\nDifference in hits: {:d}'.format(len(self.df_in), len(self.df_out), missing_hits), to_log = self.save_log)
		self.update_missing(missing_hits, first = True)

		print(self.df)


	########## Check hits with too long readout times ###############

	# Readout time correction:
	#	For LHCb count hits as missing that have a readout time over 89.1 us (max_rotime)
	
	def check_rotimes_over_max(self):
		self.df['T_readout_over_max'] = self.df['T_readout'] > self.max_rotime
		self.df_rotimes_over_max = self.df[self.df['T_readout_over_max'] == True]

		# Report results
		if self.print_log:
			self.helper.msg('READOUT TIME EXCEEDING MAX...', type = 'head', to_log = self.save_log)
			self.helper.give_update('Hits taking too long to read out: {:d}'.format(len(self.df_rotimes_over_max), to_log = self.save_log))


	########## Check hits that came at end of simulation ############

	# End of simulation correction:
	# 	Hits in the last 3564 (for LHCb) events have to be removed from analysis
	# 	as they have a higher chance to be detected and be read out quickly

	def check_end_of_sim_hits(self):

		# Last valid event
		self.lastTS = self.df_in['TS'].max() - self.lhc.max_event

		# Flag hits after last valid event
		self.df['EndOfSim'] = self.df['TS'] > self.lastTS
		self.df_endofsim = self.df[self.df['EndOfSim'] == True]
		self.df_endofsim_detected = self.df_endofsim[self.df_endofsim['Exist'] == 'both']

		# Remove end of sim hits from data set
		#self.df = self.df[self.df['EndOfSim'] == False]

		# Report results
		if self.print_log:
			self.helper.msg('END OF SIM HITS...', type = 'head', to_log = self.save_log)
			self.helper.msg('Hits in the last 3564 events:)', to_log = self.save_log)
			self.helper.printing(self.df_endofsim, to_log = self.save_log)
			self.helper.give_update('Hits read out in the last {:d} events: {:d}'.format(self.lhc.max_event, len(self.df_endofsim), to_log = self.save_log))

	########## Check hits read out after sim finished ###############

	def check_post_sim_hits(self):

		# Last event of simulation
		self.maxTS = self.df_in['TS'].max()

		#  Flag hits read out after sim finished
		self.df['PostSim'] = self.df['TS_read'] > self.maxTS
		self.df_postsim = self.df[self.df['PostSim'] == True]

		# Report results
		if self.print_log:
			self.helper.msg('POST SIM HITS...', type = 'head', to_log = self.save_log)
			self.helper.give_update('Hits read out after simulation finished: {:d}'.format(len(self.df_postsim), to_log = self.save_log))


	########## Check repeat hits ####################################

	def check_repeat_hits(self):

		# Count how often a hit occurs in same pixel at same TS
		df_count = self.df[['Col', 'Row', 'TS']].value_counts().to_frame().reset_index()
		df_count.columns = ['Col', 'Row', 'TS', 'Count']
		df_count = df_count.sort_values(by = ['Col', 'Row', 'TS'])

		# Pixels hit more than once in same event (same TS)
		df_repeats = df_count[df_count['Count'] > 1].reset_index(drop = True)
		repeats = len(df_repeats)

		# Report results
		if self.print_log:
			self.helper.msg('REPEAT HITS...', type = 'head', to_log = self.save_log)
			self.helper.msg('Repeat hits in same pixel at same time: \n(Will not be counted as separate hits by MightyPix)', to_log = self.save_log)
			self.helper.printing(df_repeats, to_log = self.save_log)
			self.helper.give_update('Repeat hits that were not measured: {:d}'.format(repeats), to_log = self.save_log)
		self.update_missing(repeats)


	########## Check pixels hit in multiple events ##################

	def check_multi_hits(self):

		# Check pixels hit more than once over all events (different TS)
		df_multi_pixels = self.df[['Col', 'Row']].value_counts().to_frame().reset_index()
		df_multi_pixels.columns = ['Col', 'Row', 'Count']
		df_multi_pixels = df_multi_pixels[df_multi_pixels.Count > 1]
		df_multi_pixels = df_multi_pixels.sort_values(by = ['Col', 'Row']).reset_index(drop = True)

		# Add info on multi hits to original dataset
		df_multi_hits = pd.merge(self.df[['Col','Row','TS', 'TS_read', 'T_in', 'T_out', 'T_readout', 'TS_readout', 'ToT_in', 'ToT_out', 'Exist']], df_multi_pixels, on = ['Col', 'Row'], how = 'outer', indicator = 'Multi')
		df_multi_hits['Multi'].replace({'left_only': 'no', 'right_only': 'yes'}, inplace=True)
		df_multi_hits = df_multi_hits[df_multi_hits['Multi'] == 'both']
		df_multi_hits = df_multi_hits.sort_values(by = ['Col', 'Row', 'TS']).reset_index(drop = True)
		df_multi_hits = df_multi_hits.drop(columns = ['Count', 'Multi', 'T_in', 'T_out'])

		# Calculate TS difference between consecutive hits
		df_multi_hits['TS_diff_to_prev'] = df_multi_hits.TS.diff()

		# Flag hits that were the first to arrive in a pixel, TS_diff_to_prev not valid
		df_multi_hits['FirstHit'] = np.select([df_multi_hits.Col.diff() != 0, df_multi_hits.Row.diff() != 0], [True, True], default = False)
		df_multi_hits.loc[df_multi_hits['FirstHit'] == True, 'TS_diff_to_prev'] = np.nan

		# Report results
		if self.print_log:
			self.helper.msg('MULTIPLE HITS...', type='head', to_log = self.save_log)
			self.helper.msg('Pixels hit in more than one event:', to_log = self.save_log)
			self.helper.printing(df_multi_hits, to_log = self.save_log)		# Quite a long list of all hits that occur in pixels that are hit multiple times
			self.helper.give_update('Number of pixels hit multiple times: {:d}\nNumber of hits not detected: {:d}'.format(len(df_multi_pixels), len(df_multi_hits[df_multi_hits['Exist'] != 'both'])), to_log = self.save_log)

		self.df_multi_hits = df_multi_hits


	########## Chek hits that came within previous hit's ToT ########

	def check_within_prev_tot(self):

		# Flag hits that arrived within ToT of prev hit
		self.df_withinToT = self.df_multi_hits[abs(self.df_multi_hits.TS_diff_to_prev) < (self.tot / (self.lhc.event_length_ns / 1000))]

		if self.print_log:
			self.helper.msg('HITS ARRIVING WITHIN TOT OF PREVIOUS HIT...', type='head', to_log = self.save_log)
			if len(self.df_withinToT) > 0:
				self.helper.msg('Hits that arrived within the previous hit ToT:', to_log = self.save_log)
				self.helper.printing(self.df_withinToT, to_log = self.save_log)
			self.helper.give_update('Hits within previous ToT: {:d}'.format(len(self.df_withinToT)), to_log = self.save_log)


	########## Check hits with readout time < ToT ###################

	def check_rotime_below_tot(self):

		df_rotime_below_tot = self.df[self.df['T_readout'] < 2000]
		
		if self.print_log:
			self.helper.msg('HITS WITH READOUT TIME BELOW TOT...', type='head', to_log = self.save_log)
			if len(df_rotime_below_tot) > 0:
				self.helper.msg('Hits that were readout faster than their tot:', to_log = self.save_log)
				self.helper.printing(df_rotime_below_tot, to_log = self.save_log)
			self.helper.give_update('Hits with readout time below ToT: {:d}'.format(len(df_rotime_below_tot)), to_log = self.save_log)


	########## Check missing hits ###################################

	def check_missing_hits(self):

		# Hits that only show in in- or outgoing dataset
		df_mb_missing = self.df.loc[self.df['Exist'] != 'both', ('Col', 'Row', 'TS', 'Exist')].reset_index(drop=True)
		df_mb_missing_in = df_mb_missing[df_mb_missing['Exist'] == 'in'].drop(columns=['Exist'])
		df_mb_missing_out = df_mb_missing[df_mb_missing['Exist'] == 'out'].drop(columns=['Exist'])
		df_mb_missing = df_mb_missing.drop(columns=['Exist'])

		# Pixels that show up in both sets of missing hits and probably have wrong TS
		df_wrong = pd.merge(df_mb_missing_in, df_mb_missing_out, on = ['Col', 'Row'], how = 'inner', indicator = 'Exist', suffixes=('_in', '_out'))
		df_wrong = df_wrong[['Col', 'Row', 'TS_in', 'TS_out']]
		df_wrong['TS_diff_in_out'] = df_wrong['TS_out'] - df_wrong['TS_in']
		df_wrong = df_wrong[df_wrong['TS_diff_in_out'] > 0]
		df_wrong = df_wrong[df_wrong['TS_diff_in_out'] < 10]

		df_wrong_in = df_wrong[['Col', 'Row', 'TS_in']].rename(columns={'TS_in': 'TS'})
		df_wrong_out = df_wrong[['Col', 'Row', 'TS_out']].rename(columns={'TS_out': 'TS'})
		df_wrong_all = pd.concat([df_wrong_in, df_wrong_out])

		# Drop hits with wrong TS from missing list
		df_missing = pd.merge(df_mb_missing, df_wrong_all, on = ['Col', 'Row', 'TS'], how = 'outer', indicator = 'Exist')
		df_missing = df_missing[df_missing['Exist'] == 'left_only'].sort_values(by=['Col', 'Row', 'TS']).drop (columns=['Exist']).reset_index(drop=True)
		df_missing.drop_duplicates(inplace=True) # Remove repeat hits from list

		# Report results
		if self.print_log:
			self.helper.msg('MISSING HITS...', type='head', to_log = self.save_log)
			if len(df_wrong) > 0:
				self.helper.msg('Hits with probably just wrong timestamp:', to_log = self.save_log)
				self.helper.printing(df_wrong, to_log = self.save_log)
			self.helper.give_update('Hits with wrong timestamp: {:d}'.format(len(df_wrong)), to_log = self.save_log)
			self.helper.msg('Hits that are actually missing:', to_log = self.save_log)
			self.helper.printing(df_missing, to_log = self.save_log)
			self.helper.give_update('Missing hits: {:d}'.format(len(df_missing)), to_log = self.save_log)
		self.update_missing(len(df_missing))

		self.df_missing = df_missing
		self.df_wrongTS = df_wrong

	########## Report different efficiencies ########################

	def report_efficiencies(self):

		eff_total = len(self.df_in)
		eff_detected = len(self.df_in) - len(self.df_missing)
		eff = 100*(eff_detected/eff_total)
		eff_low, eff_high = self.ci.get(total = eff_total, detected = eff_detected)

		rt_corr_eff_total = len(self.df_in)
		rt_corr_eff_detected = len(self.df_in) - len(self.df_missing) - len(self.df_rotimes_over_max)
		rt_corr_eff = 100*(rt_corr_eff_detected/rt_corr_eff_total)
		rt_corr_eff_low, rt_corr_eff_high = self.ci.get(total = rt_corr_eff_total, detected = rt_corr_eff_detected)

		eos_corr_eff_total = len(self.df_in) - len(self.df_endofsim)
		eos_corr_eff_detected = len(self.df_in) - len(self.df_missing) - len(self.df_endofsim_detected)

		# print(eos_corr_eff_total, eos_corr_eff_detected)
		# eos_corr_eff = 100*(eos_corr_eff_detected/eos_corr_eff_total)
		# eos_corr_eff_low, eos_corr_eff_high = self.ci.get(total = eos_corr_eff_total, detected = eos_corr_eff_detected)

		ro_eff_total = len(self.df_in) - len(self.df_withinToT)
		ro_eff_detected = len(self.df_in) - len(self.df_missing)
		ro_eff = 100*(ro_eff_detected/ro_eff_total)
		ro_eff_low, ro_eff_high = self.ci.get(total=ro_eff_total, detected=ro_eff_detected)

		rt_corr_ro_eff_total = len(self.df_in) - len(self.df_withinToT)
		rt_corr_ro_eff_detected = len(self.df_in) - len(self.df_missing) - len(self.df_rotimes_over_max)
		rt_corr_ro_eff = 100*(rt_corr_ro_eff_detected/rt_corr_ro_eff_total)
		rt_corr_ro_eff_low, rt_corr_ro_eff_high = self.ci.get(total = rt_corr_ro_eff_total, detected = rt_corr_ro_eff_detected)

		if self.print_log:
			self.helper.give_update('Efficiency: {:f} [{:f}, {:f}]'.format(eff, eff_low, eff_high), to_log = self.save_log)
			self.helper.give_update('RT corrected efficiency: {:f} [{:f}, {:f}]'.format(rt_corr_eff, rt_corr_eff_low, rt_corr_eff_high), to_log = self.save_log)
			# self.helper.give_update('EoS corrected efficiency: {:f} [{:f}, {:f}]'.format(eos_corr_eff, eos_corr_eff_low, eos_corr_eff_high), to_log = self.save_log)			
			self.helper.give_update('Readout efficiency: {:f} [{:f}, {:f}]'.format(ro_eff, ro_eff_low, ro_eff_high), to_log = self.save_log)
			self.helper.give_update('Corrected readout efficiency: {:f} [{:f}, {:f}]'.format(rt_corr_ro_eff, rt_corr_ro_eff_low, rt_corr_ro_eff_high), to_log = self.save_log)


	########## Save missing hits per sensor in extra file ###########

	def save_sensor_efficiencies(self):

		if not path.exists(self.efficiency_file.whole):
			data = pd.read_csv(self.empty_efficiency_file.whole)
		else:
			data = pd.read_csv(self.efficiency_file.whole)
			
		df_sensors = pd.DataFrame(data)

		df_sensors = df_sensors.set_index('Sensor')
		df_sensors.loc[self.sensor, 'MissingHits'] = len(self.df_missing)
		df_sensors.loc[self.sensor, 'TotalHits'] = len(self.df_in)

		df_sensors.to_csv(self.efficiency_file.whole)


	########## Save hit rates for random data #######################

	def save_rates(self):
		df = pd.DataFrame(columns = ['Rate', 'TotalHits', 'MissingHits', 'WithinPrevToT', 'OverMaxROTime','PostSim'])

		df.Rate = [float(self.rate)]
		df.TotalHits = [int(len(self.df_in))]
		df.MissingHits = [int(len(self.df_missing))]
		df.WithinPrevToT = [int(len(self.df_withinToT))]
		df.OverMaxROTime = [int(len(self.df_rotimes_over_max))]
		df.PostSim = [int(len(self.df_postsim))]
		df = df.set_index('Rate')

		if path.exists(self.rate_file.whole):
			data = pd.read_csv(self.rate_file.whole)
			df_old = pd.DataFrame(data)
			df_old = df_old.set_index('Rate')
			df_new = pd.concat([df_old, df])
			df_new.to_csv(self.rate_file.whole)
		else:
			df.to_csv(self.rate_file.whole)


	########## Compare ToTs #########################################

	def check_tots(self):

		# ToTs rounded down to multiples of 25
		self.df_in['ToT'] = (self.df_in['ToT'] * 1000 - self.df_in['ToT'] * 1000 % self.lhc.event_length_ns)

		# ToTs rounded up/down to nearest multiple of 25
		# df_in['ToT'] = round((df_in['ToT'] * 1000)/self.lhc.event_length_ns)*self.lhc.event_length_ns

		# Merge data sets with ToT
		df_tot = pd.merge(self.df_in, self.df_out, on = ['Col', 'Row', 'TS'], how = 'inner', indicator = 'Exist', suffixes=('_in', '_out'))
		df_tot['Exist'].replace({'left_only': 'in', 'right_only': 'out'}, inplace=True)
		df_tot['ToT_diff_in_out'] = df_tot['ToT_in'] - df_tot['ToT_out']

		correct_tot = (df_tot.ToT_diff == 0).sum()
		wrong_tot = (df_tot.ToT_diff != 0).sum()

		# Count how many wrong ToTs - check diff
		df_wrong_tot = df_tot[df_tot['ToT_diff_in_out'] != 0]
		df_count_wrong_tot = df_wrong_tot['ToT_diff_in_out'].value_counts().to_frame().reset_index()
		df_count_wrong_tot.columns = ['ToT_diff_in_out', 'Count']

		# Count how many negative ToTs
		neg_tot_in = len(df_tot[df_tot['ToT_in'] < 0])
		neg_tot_out = len(df_tot[df_tot['ToT_out'] < 0])

		# Report results
		if self.print_log:
			self.helper.msg('HITS WITH WRONG TOT...', type='head', to_log = self.save_log)
			self.helper.give_update('Hits with correct ToT: {:d}\nHits with wrong ToT: {:d}'.format(correct_tot, wrong_tot), to_log = self.save_log)
			self.helper.give_update('Ingoing hits with negative ToT: {:d}\nOutgoing hits with negative ToT: {:d}'.format(neg_tot_in, neg_tot_out), to_log = self.save_log)

			# Print all wrong ToTs if there are less than 10
			if (wrong_tot < 10):
				self.helper.msg('Hits with wrong ToT:', to_log = self.save_log)
				self.helper.printing(df_wrong_tot, to_log = self.save_log)


	########## Check data rate ######################################

	def check_data_rate(self):

		events = self.last_event
		event_length = self.lhc.event_length	# in s
		chip_area = self.chip.size / 100 # in cm^2
		in_rate  = len(self.df_in)  / (chip_area * events * event_length) / 1000000 # MHz/cm^2
		out_rate = len(self.df_out) / (chip_area * events * event_length) / 1000000 # MHz/cm^

		# Report results
		if self.print_log:
			self.helper.msg('DATA RATE...', type='head', to_log = self.save_log)
			self.helper.give_update('Ingoing hit rate: {:2.4f} MHz/cm^2\nOutgoing hit rate: {:2.4f} MHz/cm^2'.format(in_rate, out_rate))


###################################################################################################

	########## Plot hits ############################################

	def hitmap_hits(self, save = True, show = True):

		hits_file = File('../plots/hits/', 'hitmap_hits', self.file_name, '.pdf')
		hits_plot = Plot(hits_file, save, show)

		hits_plot.set_pixels(self.chip.cols, self.chip.rows)
		hits_plot.set_hits(self.min_hits, self.max_hits, self.cutoff)

		hits_plot.hitmap(self.df_in)


	########## Plot missing hits ####################################

	def hitmap_missing_hits(self, save = True, show = True):

		missing_hits_file = File('../plots/missing_hits/', 'hitmap_missing_hits', self.file_name, '.pdf')
		missing_hits_plot = Plot(missing_hits_file, save, show)

		missing_hits_plot.set_pixels(self.chip.cols, self.chip.rows)
		missing_hits_plot.set_hits(self.min_hits, self.max_hits, self.cutoff)

		missing_hits_plot.hitmap(self.df_missing)


	########## Plot hits vs column and row ##########################

	def plot_colrow_vs_hits(self, col = True, row = True, save = True, show = True):

		# Plot hits vs column
		if col:
			col_vs_hits_file = File('../plots/colrow_vs_hits/', 'plot_col_vs_hits', self.file_name, '.pdf')
			col_vs_hits_plot = Plot(col_vs_hits_file, save, show)

			col_vs_hits_plot.set_x(self.df_in['Col'].to_numpy(), 'Column')
			col_vs_hits_plot.y_label = 'Hits'

			col_vs_hits_plot.histo()

		# Plot hits vs row
		if row:
			row_vs_hits_file = File('../plots/colrow_vs_hits/', 'plot_row_vs_hits', self.file_name, '.pdf')
			row_vs_hits_plot = Plot(row_vs_hits_file, save, show)

			row_vs_hits_plot.set_x(self.df_in['Row'].to_numpy(), 'Row')
			row_vs_hits_plot.y_label = 'Hits'

			row_vs_hits_plot.histo()

		if col and row:

			colrow_vs_hits_file = File('../plots/hits/', 'hitmap_colrow_vs_hits', self.file_name, '.pdf')
			colrow_vs_hits_plot = Plot(colrow_vs_hits_file, save, show)

			if self.rate == 40:
				vmax = 77#80
				col_max = 16000
				row_max = 1500
			elif self.rate == 17:
				vmax = ''#40
				col_max = 7000
				row_max = 700
			elif self.data_type == 'scifi' and self.last_layer == 1:
				vmax = 2
				col_max = 15
				row_max = 5
				ticks = ticks = np.arange(0,vmax+1)
			elif self.data_type == 'scifi' and self.last_layer == 6:
				vmax = 3
				col_max = 60
				row_max = 12
				ticks = ticks = np.arange(0,vmax+1)


			try:

				colrow_vs_hits_plot.hitmap_side_histos(x = self.df_in['Col'].to_numpy(),
										   			   y = self.df_in['Row'].to_numpy(),
													   xbins = np.arange(0,30),
													   ybins = np.arange(0,321),
													   vmax = vmax,
													   col_max = col_max,
													   row_max = row_max,
													   ticks=ticks)
			
			except:

				colrow_vs_hits_plot.hitmap_side_histos(x = self.df_in['Col'].to_numpy(),
										   			   y = self.df_in['Row'].to_numpy(),
													   xbins = np.arange(0,30),
													   ybins = np.arange(0,321))




	########## Plot missing hits vs column and row ##################

	def plot_colrow_vs_missing_hits(self, col = True, row = True, save = True, show = True):

		# Plot missing hits vs column
		if col:
			col_vs_missinghits_file = File('../plots/colrow_vs_missing_hits/', 'plot_col_vs_missing_hits', self.file_name, '.pdf')
			col_vs_missinghits_plot = Plot(col_vs_missinghits_file, save, show)

			col_vs_missinghits_plot.set_x(self.df_in['Col'].to_numpy(), 'Column')
			col_vs_missinghits_plot.set_x2(self.df_missing['Col'].to_numpy(), 'Column')
			col_vs_missinghits_plot.y_label = 'Total Hits'
			col_vs_missinghits_plot.y2_label = 'Missing Hits'
			col_vs_missinghits_plot.bins = self.chip.cols
			# col_vs_missinghits_plot.y_max = 200

			col_vs_missinghits_plot.histo_two_in_one()

		# Plot missing hits vs row
		if row:
			row_vs_missinghits_file = File('../plots/colrow_vs_missing_hits/', 'plot_row_vs_missing_hits', self.file_name, '.pdf')
			row_vs_missinghits_plot = Plot(row_vs_missinghits_file, save, show)

			row_vs_missinghits_plot.set_x(self.df_in['Row'].to_numpy(), 'Row')
			row_vs_missinghits_plot.set_x2(self.df_missing['Row'].to_numpy(), 'Row')
			row_vs_missinghits_plot.y_label = 'Total Hits'
			row_vs_missinghits_plot.y2_label = 'Missing Hits'
			row_vs_missinghits_plot.bins = self.chip.rows
			# row_vs_missinghits_plot.y_max = 30

			row_vs_missinghits_plot.histo_two_in_one()

		if col and row:

			colrow_vs_missinghits_file = File('../plots/missing_hits/', 'hitmap_colrow_vs_missing_hits', self.file_name, '.pdf')
			colrow_vs_missinghits_plot = Plot(colrow_vs_missinghits_file, save, show)

			if self.rate == 40:
				vmax = 77#80
				col_max = 16000
				row_max = 1500
				ticks = []
			elif self.rate == 17:
				vmax = 3
				col_max = 45
				row_max = 9
				ticks = np.arange(0,vmax+1)
			elif self.data_type == 'scifi' and self.last_layer == 1:
				vmax = 1
				col_max = 1.5
				row_max = 1.5
				ticks = np.arange(0,vmax+1)
			elif self.data_type == 'scifi' and self.last_layer == 6:
				vmax = 3
				col_max = 3
				row_max = 2
				ticks = np.arange(0,vmax+1)

			try:
			
				colrow_vs_missinghits_plot.hitmap_side_histos(x = self.df_missing['Col'].to_numpy(),
															  y = self.df_missing['Row'].to_numpy(),
															  xbins = np.arange(0,30),
															  ybins = np.arange(0,321),
															  vmax = vmax,
															  col_max = col_max,
															  row_max = row_max,
															  ticks = ticks)
				
			except:

				colrow_vs_missinghits_plot.hitmap_side_histos(x = self.df_missing['Col'].to_numpy(),
															  y = self.df_missing['Row'].to_numpy(),
															  xbins = np.arange(0,30),
															  ybins = np.arange(0,321))


	########## Plot hits over all events incl missing hits ##########

	def plot_events_vs_hits(self, save = True, show = True):

		events_vs_hits_file = File('../plots/events_vs_hits/', 'plot_events_vs_hits', self.file_name, '.pdf')
		events_vs_hits_plot = Plot(events_vs_hits_file, save, show)

		events_vs_hits_plot.set_x( self.df_in['TS'].to_numpy(), 'Event')
		events_vs_hits_plot.set_x2(self.df_missing['TS'].to_numpy(), 'Event')

		events_vs_hits_plot.y_label = 'Simulated Hits'
		events_vs_hits_plot.y2_label = 'Missing Hits'

		events_vs_hits_plot.bins = 100

		events_vs_hits_plot.histo_two_in_one()
		#events_vs_hits_plot.histo_two_in_one_wide()



	# OLD VERSION WITH PLOT INSTEAD OF HISTO

	def old_plot_events_vs_hits(self, save = True, show = True):

		# Count total hits per event
		measured_events = self.df_in['TS'].to_numpy()
		all_events = np.arange(1, measured_events.max()+1, 1)
		all_counts = np.zeros(measured_events.max())

		for i in range(len(measured_events)):
			all_counts[measured_events[i]-1] += 1

		# Count missing hits per event
		missing_events = self.df_missing['TS'].to_numpy()
		missing_counts = np.zeros(measured_events.max())

		for i in range(len(missing_events)):
			missing_counts[missing_events[i]-1] += 1

		# Plot both in same plot
		events_vs_hits_file = File('../plots/events_vs_hits/', 'plot_events_vs_hits', self.file_name, '.pdf')
		events_vs_hits_plot = Plot(events_vs_hits_file, save, show)

		events_vs_hits_plot.set_x(all_events, 'Event')
		events_vs_hits_plot.set_y(all_counts, 'Incoming Hits')
		events_vs_hits_plot.set_y2(missing_counts, 'Missing Hits')

		events_vs_hits_plot.bins = np.arange(0,500000,5000)

		events_vs_hits_plot.two_in_one()


	########## Plot hits per event ##################################

	def plot_hits_per_event(self, save = True, show = True):

		# Count total hits per event
		measured_events = self.df_in['TS'].to_numpy()
		all_events = np.arange(1, measured_events.max()+1, 1)
		all_counts = np.zeros(measured_events.max())

		for i in range(len(measured_events)):
			all_counts[measured_events[i]-1] += 1

		# Plot both in same plot
		hits_per_event_file = File('../plots/events_vs_hits/', 'plot_hits_per_event', self.file_name, '.pdf')
		hits_per_event_plot = Plot(hits_per_event_file, save, show)

		hits_per_event_plot.set_x(all_counts, 'Hits per event')
		hits_per_event_plot.y_label = 'Counts'

		hits_per_event_plot.bins = np.amax(all_counts.astype(int))

		hits_per_event_plot.histo()


	########## Plot hits over all events incl wrong TS ##############

	def plot_events_vs_wrongTS(self, save = True, show = True):

		# Count total hits per event
		measured_events = self.df_in['TS'].to_numpy()
		all_events = np.arange(1, measured_events.max()+1, 1)
		all_counts = np.zeros(measured_events.max())

		for i in range(len(measured_events)):
			all_counts[measured_events[i]-1] += 1

		# Count missing hits per event
		missing_events = self.df_wrongTS['TS_in'].to_numpy()
		missing_counts = np.zeros(measured_events.max())

		for i in range(len(missing_events)):
			missing_counts[missing_events[i]-1] += 1

		# Plot both in same plot
		events_vs_wrongTS_file = File('../plots/events_vs_hits/', 'plot_events_vs_wrongTS', self.file_name, '.pdf')
		events_vs_wrongTS_plot = Plot(events_vs_wrongTS_file, save, show)

		events_vs_wrongTS_plot.set_x(all_events, 'Event')
		events_vs_wrongTS_plot.set_y(all_counts, 'Incoming Hits')
		events_vs_wrongTS_plot.set_y2(missing_counts, 'Hits with Wrong Time Stamp')

		events_vs_wrongTS_plot.two_in_one()


########## Plot readout times over all events #######################

	def plot_events_vs_rotime(self, save = True, show = True):

		measured_df = self.df[self.df['T_readout_over_max'] == False]
		measured_events = measured_df['TS'].to_numpy()
		measured_rotimes = measured_df['T_readout'].to_numpy()
		measured_rotimes = measured_rotimes / 1000.
		all_events = np.arange(1, measured_events.max()+1, 1)
		rotimes = np.zeros(measured_events.max())
		nr_of_rotimes = np.zeros(measured_events.max())

		for i in range(len(measured_events)):
			if measured_rotimes[i] > 0. and not np.isnan(measured_rotimes[i]):
				rotimes[measured_events[i]-1] += measured_rotimes[i]
				nr_of_rotimes[measured_events[i]-1] += 1

		mean_rotimes = np.divide(rotimes, nr_of_rotimes, out=np.zeros_like(rotimes), where=nr_of_rotimes!=0)

		events_vs_rotimes_file = File('../plots/events_vs_hits/', 'plot_events_vs_rotime', self.file_name, '.pdf')
		events_vs_rotimes_plot = Plot(events_vs_rotimes_file, save, show)
		events_vs_rotimes_plot.set_x(all_events, 'Event')
		events_vs_rotimes_plot.set_y(mean_rotimes, 'Mean Readout Time (us)')
		events_vs_rotimes_plot.y_min = 2
		events_vs_rotimes_plot.scatter()


	########## Plot readout times ###################################

	def plot_rotime_vs_colrow(self, col = True, row = True, save = True, show = True):

		# Readout times in us
		df_ro = self.df[['Col','Row','T_readout']].dropna()
		cols = df_ro['Col'].to_numpy()
		rows = df_ro['Row'].to_numpy()
		ro_times = df_ro['T_readout'].to_numpy()
		ro_times = ro_times/1000

		# Plot readout time vs column
		if col:

			rotime_vs_col_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_vs_col', self.file_name, '.pdf')
			rotime_vs_col_plot = Plot(rotime_vs_col_file, save, show)

			rotime_vs_col_plot.set_x(ro_times, 'Readout time (us)')
			# rotime_vs_col_plot.set_x2(mean_cols, 'Readout time (us)')
			rotime_vs_col_plot.set_y(cols, 'Column')
			# rotime_vs_col_plot.x_max = 10

			rotime_vs_col_plot.sub_histo()
			# rotime_vs_col_plot.sub_histo2()

		# Plot readout time vs row
		if row:

			rotime_vs_row_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_vs_row', self.file_name, '.pdf')
			rotime_vs_row_plot = Plot(rotime_vs_row_file, save, show)

			rotime_vs_row_plot.set_x(ro_times, 'Readout time (us)')
			# rotime_vs_row_plot.set_x2(mean_rows, 'Readout time (us)')

			rotime_vs_row_plot.set_y(rows, 'Row')
			# rotime_vs_row_plot.x_max = 200

			rotime_vs_row_plot.sub_histo()
			# rotime_vs_row_plot.sub_histo2()


########## Plot readout times above max #############################

	def plot_rotime_over_max_vs_colrow(self, col = True, row = True, save = True, show = True):

		# Readout times in us
		df_ro = self.df[['Col','Row','T_readout']].dropna()
		cols = df_ro['Col'].to_numpy()
		rows = df_ro['Row'].to_numpy()
		ro_times = df_ro['T_readout'].to_numpy()
		ro_times = ro_times/1000

		df_ro_over_max = self.df[self.df['T_readout_over_max'] == True]
		df_ro_over_max = df_ro_over_max[['Col','Row','T_readout']].dropna()
		cols_ro_over_max = df_ro_over_max['Col'].to_numpy()
		rows_ro_over_max = df_ro_over_max['Row'].to_numpy()
		ro_times_over_max = df_ro_over_max['T_readout'].to_numpy()
		ro_times_over_max = ro_times_over_max/1000

		bins = np.logspace(start = np.log10(df_ro.T_readout.min()), stop = np.log10(df_ro.T_readout.max()), num = 1300)


		# Plot readout time vs column
		if col:

			rotime_vs_col_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_over_max_vs_col', self.file_name, '.pdf')
			rotime_vs_col_plot = Plot(rotime_vs_col_file, save, show)

			rotime_vs_col_plot.set_x(ro_times, 'Readout time (us)')
			rotime_vs_col_plot.set_x2(ro_times_over_max, 'Readout time over max (us)')

			rotime_vs_col_plot.set_y(cols, 'Column')
			rotime_vs_col_plot.set_y2(cols_ro_over_max, 'Column')
			# rotime_vs_col_plot.x_max = 10
			rotime_vs_col_plot.bins = bins

			rotime_vs_col_plot.sub_histo4()

		# Plot readout time vs row
		if row:

			rotime_vs_row_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_over_max_vs_row', self.file_name, '.pdf')
			rotime_vs_row_plot = Plot(rotime_vs_row_file, save, show)

			rotime_vs_row_plot.set_x(ro_times, 'Readout time (us)')
			rotime_vs_row_plot.set_x2(ro_times_over_max, 'Readout time over max (us)')

			rotime_vs_row_plot.set_y(rows, 'Row')
			rotime_vs_row_plot.set_y2(rows_ro_over_max, 'Row')
			# rotime_vs_row_plot.x_max = 200

			rotime_vs_col_plot.bins = bins

			rotime_vs_row_plot.sub_histo4()


########## Plot readout times only below max #############################

	def plot_rotime_below_max_vs_colrow(self, col = True, row = True, save = True, show = True):

		# Readout times below max (89.1 us for lhc) in us
		df_ro_below_max = self.df[self.df['T_readout_over_max'] == False]
		df_ro_below_max = df_ro_below_max[['Col','Row','T_readout']].dropna()
		cols_ro_below_max = df_ro_below_max['Col'].to_numpy()
		rows_ro_below_max = df_ro_below_max['Row'].to_numpy()
		ro_times_below_max = df_ro_below_max['T_readout'].to_numpy()
		ro_times_below_max = ro_times_below_max/1000


		# Plot readout time vs column
		if col:
			rotime_vs_col_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_below_max_vs_col', self.file_name, '.pdf')
			rotime_vs_col_plot = Plot(rotime_vs_col_file, save, show)

			rotime_vs_col_plot.set_x(ro_times_below_max, 'Readout time (us)')
			rotime_vs_col_plot.set_y(cols_ro_below_max, 'Column')

			rotime_vs_col_plot.y_min = 0
			rotime_vs_col_plot.y_max = config['chip']['columns']

			if self.data_type == 'scifi':
				rotime_vs_col_plot.bins = np.arange(2, 10.1, 0.05)

				if self.last_layer == 1:
					rotime_vs_col_plot.x_min = 2
					rotime_vs_col_plot.x_max = 6
					rotime_vs_col_plot.y2_min = 0.8
					rotime_vs_col_plot.y2_max = 3e1
					vmax = 3
					title = ''
					ticks = np.arange(1,vmax+1)

				elif self.last_layer == 6:
					rotime_vs_col_plot.x_min = 2
					rotime_vs_col_plot.x_max = 8
					rotime_vs_col_plot.y2_min = 0.8
					rotime_vs_col_plot.y2_max = 5e2
					vmax = 8
					title = ''
					ticks = []

				else:
					print('WARNING: No contraints defined. Plotting with automatic settings.')
			
			elif self.data_type == 'random':
				rotime_vs_col_plot.y2_min = 0.8
				rotime_vs_col_plot.y2_max = 1e6

				if self.rate == 17 and self.extra != 'newfsm':
					rotime_vs_col_plot.x_min = 2 - 0.02 * max(ro_times_below_max)
					rotime_vs_col_plot.x_max = max(ro_times_below_max) + 0.05 * max(ro_times_below_max)
					rotime_vs_col_plot.bins = np.arange(0, 95, 0.1)

				elif self.rate == 17 and self.extra == 'newfsm':
					rotime_vs_col_plot.x_min = 1.8
					rotime_vs_col_plot.x_max = 4.5
					rotime_vs_col_plot.bins = np.arange(1.5, 5.5, 0.01)
					
				else:
					rotime_vs_col_plot.plot_readout_time_limit = True
					rotime_vs_col_plot.x_min = 0#2 - 0.02 * max(ro_times_below_max) #0
					rotime_vs_col_plot.x_max = 95#max(ro_times_below_max) + 0.05 * max(ro_times_below_max) #95
					rotime_vs_col_plot.bins = np.arange(0.1, 95.1, 0.5)

					if self.extra != 'newfsm' and self.rate in [20,21,22]:
						vmax = 3000
						title = str(self.rate)+' MHz/cm²'
						ticks = []

					elif self.extra == 'newfsm' and self.rate in [31,32,33]:
						vmax = 6500
						title = str(self.rate)+' MHz/cm²'
						ticks = []

			else:
				print('ERROR: Data type not defined. Choose "scifi" or "random".')
				sys.exit()

			#rotime_vs_col_plot.sub_histo()

			rotime_vs_col_file2 = File('../plots/readout_time_vs_colrow/', 'hitmap_rotime_below_max_vs_col', self.file_name, '.pdf')
			rotime_vs_col_plot.file = rotime_vs_col_file2

			try:
				rotime_vs_col_plot.hitmap_sub_histo(ybins=np.arange(0,30), vmax = vmax, title = title, ticks = ticks)
			except:
				rotime_vs_col_plot.hitmap_sub_histo(ybins=np.arange(0,30))



		# Plot readout time vs row
		if row:

			rotime_vs_row_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_below_max_vs_row', self.file_name, '.pdf')
			rotime_vs_row_plot = Plot(rotime_vs_row_file, save, show)

			rotime_vs_row_plot.set_x(ro_times_below_max, 'Readout time (us)')
			rotime_vs_row_plot.set_y(rows_ro_below_max, 'Row')

			rotime_vs_row_plot.y_min = 0
			rotime_vs_row_plot.y_max = config['chip']['rows']

			if self.data_type == 'scifi':
				rotime_vs_row_plot.bins = np.arange(2, 10.1, 0.05)

				if self.last_layer == 1:
					rotime_vs_row_plot.x_min = 2
					rotime_vs_row_plot.x_max = 6
					rotime_vs_row_plot.y2_min = 0.8
					rotime_vs_row_plot.y2_max = 3e1
					vmax = 2
					title = ''
					ticks = np.arange(0,vmax+1)

				elif self.last_layer == 6:
					rotime_vs_row_plot.x_min = 2
					rotime_vs_row_plot.x_max = 8
					rotime_vs_row_plot.y2_min = 0.8
					rotime_vs_row_plot.y2_max = 5e2
					vmax = 3
					title = ''
					ticks = np.arange(0,vmax+1)

				else:
					print('WARNING: No contraints defined. Plotting with automatic settings.')

			elif self.data_type == 'random':
				rotime_vs_row_plot.y2_min = 0.8
				rotime_vs_row_plot.y2_max = 1e6

				if self.rate == 17 and self.extra != 'newfsm':
					rotime_vs_row_plot.x_min = 2 - 0.02 * max(ro_times_below_max)
					rotime_vs_row_plot.x_max = max(ro_times_below_max) + 0.05 * max(ro_times_below_max)
					rotime_vs_row_plot.bins = np.arange(0, 95, 0.1)

				elif self.rate == 17 and self.extra == 'newfsm':
					rotime_vs_row_plot.x_min = 1.8
					rotime_vs_row_plot.x_max = 4.5
					rotime_vs_row_plot.bins = np.arange(1.5, 5.5, 0.01)

				else:
					rotime_vs_row_plot.plot_readout_time_limit = True
					rotime_vs_row_plot.x_min = 0#2 - 0.02 * max(ro_times_below_max) #0
					rotime_vs_row_plot.x_max = 95#max(ro_times_below_max) + 0.05 * max(ro_times_below_max) #95
					rotime_vs_row_plot.bins = np.arange(0.1, 95.1, 0.5)

					if self.extra != 'newfsm' and self.rate in [20,21,22]:
						vmax = 300
						title = str(self.rate)+' MHz/cm²'
						ticks = []

					elif self.extra == 'newfsm' and self.rate in [31,32,33]:
						vmax = 550
						title = str(self.rate)+' MHz/cm²'
						ticks = []

			else:
				print('ERROR: Data type not defined. Choose "scifi" or "random".')
				sys.exit()


			#rotime_vs_row_plot.sub_histo()

			rotime_vs_row_file2 = File('../plots/readout_time_vs_colrow/', 'hitmap_rotime_below_max_vs_row', self.file_name, '.pdf')
			rotime_vs_row_plot.file = rotime_vs_row_file2

			try:
				rotime_vs_row_plot.hitmap_sub_histo(ybins=np.arange(0,321), vmax = vmax, title = title, ticks = ticks)
			except:
				rotime_vs_row_plot.hitmap_sub_histo(ybins=np.arange(0,321))


		#if col and row:
		if False:

			rotime_vs_colrow_file = File('../plots/readout_time_vs_colrow/', 'hitmap_rotime_below_max_vs_colrow', self.file_name, '.pdf')
			rotime_vs_colrow_plot = Plot(rotime_vs_colrow_file, save, show)

			if self.rate == 40:
				vmax = 80
				col_max = 16000
				row_max = 1500
				ticks = []
			elif self.rate == 17:
				vmax = 3
				col_max = 45
				row_max = 9
				ticks = np.arange(0,vmax+1)
			
			rotime_vs_colrow_plot.hitmap_side_histos(x = cols_ro_below_max,
												     y = rows_ro_below_max,
													 xbins = np.arange(0,30),
													 ybins = np.arange(0,321)#,
													 #vmax = vmax,
													 #col_max = col_max,
													 #row_max = row_max,
													 #ticks = ticks
													 )



	########## Plot after sim readout times #########################

	def plot_postsim_rotime_vs_colrow(self, col = True, row = True, save = True, show = True):

		# Readout times in us
		df_ro = self.df[['Col','Row','T_readout']].dropna()
		cols = df_ro['Col'].to_numpy()
		rows = df_ro['Row'].to_numpy()
		ro_times = df_ro['T_readout'].to_numpy()
		ro_times = ro_times/1000

		df_ro_postsim = self.df[self.df['PostSim'] == True]
		df_ro_postsim = df_ro_postsim[['Col','Row','T_readout']].dropna()
		cols_postsim = df_ro_postsim['Col'].to_numpy()
		rows_postsim = df_ro_postsim['Row'].to_numpy()
		ro_times_postsim = df_ro_postsim['T_readout'].to_numpy()
		ro_times_postsim = ro_times_postsim/1000


		# Plot readout time vs column
		if col:

			rotime_vs_col_file = File('../plots/readout_time_vs_colrow/', 'plot_postsim_rotime_vs_col', self.file_name, '.pdf')
			rotime_vs_col_plot = Plot(rotime_vs_col_file, save, show)

			rotime_vs_col_plot.set_x(ro_times, 'Readout time (us)')
			rotime_vs_col_plot.set_x2(ro_times_postsim, 'Readout time for post sim hits (us)')

			rotime_vs_col_plot.set_y(cols, 'Column')
			rotime_vs_col_plot.set_y2(cols_postsim, 'Column')
			# rotime_vs_col_plot.x_max = 10

			rotime_vs_col_plot.sub_histo3()

		# Plot readout time vs row
		if row:

			rotime_vs_row_file = File('../plots/readout_time_vs_colrow/', 'plot_postsim_rotime_vs_row', self.file_name, '.pdf')
			rotime_vs_row_plot = Plot(rotime_vs_row_file, save, show)

			rotime_vs_row_plot.set_x(ro_times, 'Readout time (us)')
			rotime_vs_row_plot.set_x2(ro_times_postsim, 'Readout time for post sim hits (us)')

			rotime_vs_row_plot.set_y(rows, 'Row')
			rotime_vs_row_plot.set_y2(rows_postsim, 'Row')
			# rotime_vs_row_plot.x_max = 200

			rotime_vs_row_plot.sub_histo3()


###################################################################################################

	########## Conclusion ###########################################

	def print_conclusion(self):

		self.helper.msg('ANALYSIS CONCLUDED...', type = 'head', to_log = self.save_log)

		if sum(self.missing) == 0:
			self.helper.msg('Success! All hits were accounted for.', to_log = self.save_log)
		elif abs(sum(self.missing)) == 1:
			self.helper.msg(str(sum(self.missing))+' hit still not accounted for.', to_log = self.save_log)
		else:
			self.helper.msg(str(sum(self.missing))+' hits still not accounted for.', to_log = self.save_log)


	########## Summary ##############################################

	def print_summary(self):

		self.helper.msg('SUMMARY...\n', type='head', to_log = self.save_log)

		for item in self.helper.summary:
			self.helper.printing(item, to_log = self.save_log)

		self.helper.msg('', to_log = self.save_log)


	########## Close log file #######################################

	def end_log(self):

		self.log.close()

		print('Log saved as:\n\t{}'.format(self.log_file.name))
		print('Log saved to:\n\t{}'.format(self.log_file.path))

		self.helper.hline()

	########## Full comparison encompassing all functions ###########

	# Run analysis
	def compare(self):
		# Initialising functions
		if self.save_log: self.start_log()
		self.read_data()
		#self.ts_correction()
		if self.print_log: self.print_header()

		# Analysis
		self.check_in_out_hits()
		self.check_rotimes_over_max()
		self.check_end_of_sim_hits()
		self.check_post_sim_hits()
		self.check_repeat_hits()
		self.check_multi_hits()
		self.check_within_prev_tot()
		self.check_rotime_below_tot()
		self.check_missing_hits()
		self.report_efficiencies()
		if self.data_type == 'scifi':	self.save_sensor_efficiencies()
		if self.data_type == 'random':	self.save_rates()
		#self.check_tots()
		self.check_data_rate()

		# Concluding functions
		if self.save_log:
			self.print_conclusion()
			self.print_summary()

	# Run all plots and save them
	def plots(self, save_me, show_me):
		self.hitmap_hits(save = save_me, show = show_me)
		self.hitmap_missing_hits(save = save_me, show = show_me)
		#self.plot_colrow_vs_hits(save = save_me, show = show_me)
		self.plot_colrow_vs_missing_hits(save = save_me, show = show_me)
		self.plot_events_vs_hits(save = save_me, show = show_me)
		self.plot_hits_per_event(save = save_me, show = show_me)
		#self.plot_events_vs_wrongTS(save = save_me, show = show_me) # Don't have wrong TS so not of interest
		#self.plot_events_vs_rotime(save = save_me, show = show_me) # No interested in
		#self.plot_rotime_vs_colrow(save = save_me, show = show_me) # Includes readout times above 89.1 us, with we are not interested in
		#self.plot_rotime_over_max_vs_colrow(save = save_me, show = show_me)
		self.plot_rotime_below_max_vs_colrow(save = save_me, show = show_me)
		#self.plot_postsim_rotime_vs_colrow(save = save_me, show = show_me)

	def plot_selected(self, save_me, show_me):
		#pass
		#self.plot_colrow_vs_hits(save = save_me, show = show_me)
		#self.plot_colrow_vs_missing_hits(save = save_me, show = show_me)
		#self.plot_events_vs_hits(save = save_me, show = show_me)
		self.plot_rotime_below_max_vs_colrow(save = save_me, show = show_me)

	def end(self):
		if self.save_log: self.end_log()

	# Run everything
	def run(self):
		self.compare()
		if self.show_plots or self.save_plots:
			if config['plots']['selected']:
				self.plot_selected(save_me=self.save_plots, show_me=self.show_plots)
			else:
				self.plots(save_me=self.save_plots, show_me=self.show_plots)
		self.end()