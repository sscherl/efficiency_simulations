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

from classes.simulation_data_plot import Plot
from classes.toolkit import Chip, File, LHCParameters, Helper

# pd.set_option("display.max_rows", None, "display.max_columns", None)

###### Classes ##################################################

class DataComparison:

	###### Initialisation ##########################################

	def __init__(self, chip, file_name, fsm_clock):
		self.chip = chip
		self.file_name = file_name
		self.fsm_clock = fsm_clock

		# Files
		self.input_file = File('../data/input_data/', 'data', self.file_name, '.csv')
		self.output_file = File('../data/output_data/', 'output', self.file_name+'_'+self.fsm_clock, '.csv')
		self.output_extra_file = File('../data/output_data/', 'output_extra', self.file_name+'_'+self.fsm_clock, '.csv')
		self.log_file = File('../data/log_files/', 'log', self.file_name+'_'+self.fsm_clock, '.txt')

		# For logging functions
		self.missing = []

		# For plotting hitmaps
		self.min_hits = 0
		self.max_hits = 0
		self.cutoff = False

		# For secondaries correction
		self.lhc = LHCParameters()
		self.helper = Helper()

	###### Functions ###############################################

	def update_missing(self, new_input, first = False):

		if first: self.missing.append(new_input)
		else: self.missing.append(-new_input)

		if sum(self.missing) == 0:
			self.helper.msg('>>> All hits accounted for!', to_log = True)
		elif sum(self.missing) == 1 or self.missing == -1:
			self.helper.msg('>>> {:d} hit not yet accounted for'.format(sum(self.missing)), to_log = True)
		else:
			self.helper.msg('>>> {:d} hits not yet accounted for'.format(sum(self.missing)), to_log = True)

	###### Log file ################################################

	def start_log(self):

		self.log = open(self.log_file.whole, 'w')
		self.helper.log_file = self.log

	###### Read data ###############################################

	def read_data(self):

		# Read output data (from pixelmatrix)
		self.data_out_all = pd.read_csv(self.output_file.whole)
		self.df_out_all = pd.DataFrame(self.data_out_all)
		self.df_out_all.drop(self.df_out_all.index[0], axis=0, inplace=True)
		self.df_out_all = self.df_out_all.sort_values(by=['TS']).reindex()
		self.df_out_all = self.df_out_all.rename(columns={'RO_Time': 'Treadout'})

		# Read output data containing readin times
		self.data_readin_all = pd.read_csv(self.output_extra_file.whole)
		self.df_readin_all = pd.DataFrame(self.data_readin_all)
		self.df_readin_all.drop(columns=['#','ToT', 'TS', 'Amp'])

		# Read full input data (to pixelmatrx) and add readin times
		self.data_in_all = pd.read_csv(self.input_file.whole)
		self.df_in_all = pd.DataFrame(self.data_in_all)
		self.df_in_all['Treadin'] = self.df_readin_all['RO_Time']
		self.df_in_all = self.df_in_all.rename(columns={'Event': 'TS'})

	###### TS correction for secondaries ###########################

	def ts_correction(self):

		# Original timestamps
		self.df_in_all['TS_orig'] = self.df_in_all['TS']

		# TS offsets (for secondaries)
		self.df_in_all['TS_offset'] = ((self.df_in_all['ToF'] - self.lhc.min_tof_data)/self.lhc.bx_period_ns).astype(int)

		# TS MightyPix is expected to show
		self.df_in_all['TS'] +=  self.df_in_all['TS_offset']

	###### Print header ############################################

	def print_header(self):

		self.helper.msg('ANALYSIS: Hits in the MightyPix Pixelmatrix model', type = 'head', to_log = True)
		self.helper.msg('FSM Clock speed: 40 MHz', type = 'head', to_log = True)
		self.helper.msg('Input  file:\t'+self.input_file.name+'\nOutput file:\t'+self.output_file.name, type = 'head', to_log = True)

	###### Compare in- and outgoing hits ###########################

	def check_in_out_hits(self):

		# Select input values to compare
		self.df_in = self.df_in_all[['Col','Row','TS', 'ToT', 'Treadin']]

		# Select output values to compare
		self.df_out = self.df_out_all[['Col', 'Row', 'ToT', 'TS', 'Treadout']]
		self.df_out['TS'] -= (self.df_out['TS'].min() - self.df_in['TS'].min())

		# Merge in- and outgoing hits
		self.df = pd.merge(self.df_in, self.df_out, on = ['Col', 'Row', 'TS'], how = 'outer', indicator = 'Exist')
		self.df['Exist'].replace({'left_only': 'in', 'right_only': 'out'}, inplace=True)
		# self.df = self.df.sort_values(by=['TS']).reindex()

		# Tread is time to readout the hit in the pixel in ns
		self.df['Tread'] = self.df['Treadout'] - self.df['Treadin']

		# TSread is the timestamp when the hit was read out
		self.df['TSread'] = self.df['Tread']/self.lhc.bx_period_ns + self.df['TS']

		# Difference between ingoing and outgoing hits
		missing_hits = len(self.df_in) - len(self.df_out)

		# Report results
		self.helper.msg('NUMBER OF HITS...', type = 'head', to_log = True)
		self.helper.give_update('Incoming hits: {:d}\nOutgoing hits: {:d}\nDifference in hits: {:d}'.format(len(self.df_in), len(self.df_out), missing_hits), to_log = True)
		self.update_missing(missing_hits, first=True)

	###### Check repeat hits #######################################

	def check_repeat_hits(self):

		# Count how often a hit occurs in same pixel at same TS
		df_count = self.df[['Col', 'Row', 'TS']].value_counts().to_frame().reset_index()
		df_count.columns = ['Col', 'Row', 'TS', 'Count']
		df_count = df_count.sort_values(by=['Col', 'Row', 'TS'])

		# Count hits that occure more than once in same pixel at same TS
		df_repeats = df_count[df_count['Count'] > 1].reset_index(drop=True)
		repeats = len(df_repeats)

		# Report results
		self.helper.msg('REPEAT HITS...', type = 'head', to_log = True)
		self.helper.msg('Repeat hits in same pixel at same time: \n(Will not be counted as separate hits by MightyPix)', to_log = True)
		self.helper.printing(df_repeats, to_log = True)
		self.helper.give_update('Repeat hits that were not measured: {:d}'.format(repeats), to_log = True)
		self.update_missing(repeats)

	###### Check pixels hit in multiple events #####################

	def check_multi_hits(self):

		df_multi_pixels = self.df[['Col', 'Row']].value_counts().to_frame().reset_index()
		df_multi_pixels.columns = ['Col', 'Row', 'Count']
		df_multi_pixels = df_multi_pixels[df_multi_pixels.Count > 1]
		df_multi_pixels = df_multi_pixels.sort_values(by=['Col', 'Row']).reset_index(drop=True)

		df_multi_hits = pd.merge(self.df[['Col','Row','TS', 'Treadin', 'Treadout', 'Tread', 'TSread', 'Exist']], df_multi_pixels, on = ['Col', 'Row'], how = 'outer', indicator = 'Multi')
		df_multi_hits['Multi'].replace({'left_only': 'no', 'right_only': 'yes'}, inplace=True)
		df_multi_hits = df_multi_hits[df_multi_hits['Multi'] == 'both']
		df_multi_hits = df_multi_hits.sort_values(by=['Col', 'Row', 'TS']).reset_index(drop=True)
		df_multi_hits = df_multi_hits.drop(columns=['Count', 'Multi', 'Treadin', 'Treadout'])

		# Report results
		self.helper.msg('MULTIPLE HITS...', type='head', to_log = True)
		self.helper.msg('Pixels hit in more than one event:', to_log = True)
		self.helper.printing(df_multi_hits, to_log = True)		# Quite a long list of all hits that occur in pixels that are hit multiple times
		self.helper.give_update('Number of pixels hit multiple times: {:d}\nNumber of hits not detected: {:d}'.format(len(df_multi_pixels), len(df_multi_hits[df_multi_hits['Exist'] != 'both'])), to_log = True)

	###### Check missing hits ######################################

	def check_missing_hits(self):

		# Hits that only show in in- or outgoing dataset
		df_mb_missing = self.df.loc[self.df['Exist'] != 'both', ('Col', 'Row', 'TS', 'Exist')].reset_index(drop=True)
		df_mb_missing_in = df_mb_missing[df_mb_missing['Exist'] == 'in'].drop(columns=['Exist'])
		df_mb_missing_out = df_mb_missing[df_mb_missing['Exist'] == 'out'].drop(columns=['Exist'])
		df_mb_missing = df_mb_missing.drop(columns=['Exist'])

		# Pixels that show up in both sets of missing hits and probably have wrong TS
		df_wrong = pd.merge(df_mb_missing_in, df_mb_missing_out, on = ['Col', 'Row'], how = 'inner', indicator = 'Exist', suffixes=('_in', '_out'))
		df_wrong = df_wrong[['Col', 'Row', 'TS_in', 'TS_out']]
		df_wrong['TS_diff'] = df_wrong['TS_out'] - df_wrong['TS_in']
		df_wrong = df_wrong[df_wrong['TS_diff'] > 0]
		df_wrong = df_wrong[df_wrong['TS_diff'] < 10]

		df_wrong_in = df_wrong[['Col', 'Row', 'TS_in']].rename(columns={'TS_in': 'TS'})
		df_wrong_out = df_wrong[['Col', 'Row', 'TS_out']].rename(columns={'TS_out': 'TS'})
		df_wrong_all = pd.concat([df_wrong_in, df_wrong_out])

		# Drop hits with wrong TS from missing list
		df_missing = pd.merge(df_mb_missing, df_wrong_all, on = ['Col', 'Row', 'TS'], how = 'outer', indicator = 'Exist')
		df_missing = df_missing[df_missing['Exist'] == 'left_only'].sort_values(by=['Col', 'Row', 'TS']).drop(columns=['Exist']).reset_index(drop=True)

		# Report results
		self.helper.msg('MISSING HITS...', type='head', to_log = True)
		if len(df_wrong) > 0:
			self.helper.msg('Hits with probably just wrong timestamp:', to_log = True)
			self.helper.printing(df_wrong, to_log = True)
		self.helper.give_update('Hits with wrong timestamp: {:d}'.format(len(df_wrong)), to_log = True)
		self.helper.msg('Hits that are actually missing:', to_log = True)
		self.helper.printing(df_missing, to_log = True)
		self.helper.give_update('Missing hits: {:d}'.format(len(df_missing)), to_log = True)
		self.update_missing(len(df_missing))

		self.df_missing = df_missing

	###### Compare ToTs ############################################

	def check_tots(self):

		# ToTs rounded down to multiples of 25
		self.df_in['ToT'] = (self.df_in['ToT'] * 1000 - self.df_in['ToT'] * 1000 % self.lhc.bx_period_ns)

		# ToTs rounded up/down to nearest multiple of 25
		# df_in['ToT'] = round((df_in['ToT'] * 1000)/self.lhc.bx_period_ns)*self.lhc.bx_period_ns

		# Merge data sets with ToT
		df_tot = pd.merge(self.df_in, self.df_out, on = ['Col', 'Row', 'TS'], how = 'inner', indicator = 'Exist', suffixes=('_in', '_out'))
		df_tot['Exist'].replace({'left_only': 'in', 'right_only': 'out'}, inplace=True)
		df_tot['ToT_diff'] = df_tot['ToT_in'] - df_tot['ToT_out']

		correct_tot = (df_tot.ToT_diff == 0).sum()
		wrong_tot = (df_tot.ToT_diff != 0).sum()

		# Count how many wrong ToTs - check diff
		df_wrong_tot = df_tot[df_tot['ToT_diff'] != 0]
		df_count_wrong_tot = df_wrong_tot['ToT_diff'].value_counts().to_frame().reset_index()
		df_count_wrong_tot.columns = ['ToT_diff', 'Count']

		# Count how many negative ToTs
		neg_tot_in = len(df_tot[df_tot['ToT_in'] < 0])
		neg_tot_out = len(df_tot[df_tot['ToT_out'] < 0])

		# Report results
		self.helper.msg('HITS WITH WRONG TOT...', type='head', to_log = True)
		self.helper.give_update('Hits with correct ToT: {:d}\nHits with wrong ToT: {:d}'.format(correct_tot, wrong_tot), to_log = True)
		self.helper.give_update('Ingoing hits with negative ToT: {:d}\nOutgoing hits with negative ToT: {:d}'.format(neg_tot_in, neg_tot_out), to_log = True)

		# Print all wrong ToTs if there are less than 10
		if (wrong_tot < 10):
			self.helper.msg('Hits with wrong ToT:', to_log = True)
			self.helper.printing(df_wrong_tot, to_log = True)

		# self.helper.msg('Wrongness of ToTs: \n(ToT_diff = ToT_in - ToT_out)', to_log = True)
		# self.helper.printing(df_count_wrong_tot, to_log = True)

	###### Check data rate #########################################

	def check_data_rate(self):

		# Calculate hits per event
		df_count_events_in = self.df_in['TS'].value_counts().to_frame().reset_index()
		df_count_events_in.columns = ['Event', 'Count']

		df_count_events_out = self.df_out['TS'].value_counts().to_frame().reset_index()
		df_count_events_out.columns = ['Event', 'Count']

		# Omit contributions
		# df_count_events.drop(df_count_events.index[0], axis=0, inplace=True)

		data_rate_in = df_count_events_in.Count.sum()/(500 * math.ceil(self.df_in['TS'].iloc[-1]/500))
		data_rate_out = df_count_events_out.Count.sum()/(500 * math.ceil(self.df_out['TS'].iloc[-1]/500))

		# Report results
		self.helper.msg('DATA RATE...', type='head', to_log = True)
		self.helper.give_update('Ingoing data rate: {:2.2f}\nOutgoing data rate: {:2.2f}'.format(data_rate_in, data_rate_out), to_log = True)

	###### Plot hits ###############################################

	def hitmap_hits(self, save = True, show = True):

		hits_file = File('../plots/hits/', 'hitmap_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
		hits_plot = Plot(hits_file, save, show)

		hits_plot.set_pixels(self.chip.cols, self.chip.rows)
		hits_plot.set_hits(self.min_hits, self.max_hits, self.cutoff)

		hits_plot.hitmap(self.df_in)

	###### Plot missing hits #######################################

	def hitmap_missing_hits(self, save = True, show = True):

		missing_hits_file = File('../plots/missing_hits/', 'hitmap_missing_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
		missing_hits_plot = Plot(missing_hits_file, save, show)

		missing_hits_plot.set_pixels(self.chip.cols, self.chip.rows)
		missing_hits_plot.set_hits(self.min_hits, self.max_hits, self.cutoff)

		missing_hits_plot.hitmap(self.df_missing)

	###### Plot hits vs column and row #############################

	def plot_colrow_vs_hits(self, col = True, row = True, save = True, show = True):

		# Plot hits vs column
		if col:
			col_vs_hits_file = File('../plots/colrow_vs_hits/', 'plot_col_vs_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
			col_vs_hits_plot = Plot(col_vs_hits_file, save, show)

			col_vs_hits_plot.set_x(self.df_in['Col'].to_numpy(), 'Column')
			col_vs_hits_plot.y_label = 'Hits'

			col_vs_hits_plot.histo()

		# Plot hits vs row
		if row:
			row_vs_hits_file = File('../plots/colrow_vs_hits/', 'plot_row_vs_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
			row_vs_hits_plot = Plot(row_vs_hits_file, save, show)

			row_vs_hits_plot.set_x(self.df_in['Row'].to_numpy(), 'Row')
			row_vs_hits_plot.y_label = 'Hits'

			row_vs_hits_plot.histo()

	###### Plot missing hits vs column and row #####################

	def plot_colrow_vs_missing_hits(self, col = True, row = True, save = True, show = True):

		# Plot missing hits vs column
		if col:
			col_vs_missinghits_file = File('../plots/colrow_vs_missing_hits/', 'plot_col_vs_missing_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
			col_vs_missinghits_plot = Plot(col_vs_missinghits_file, save, show)

			col_vs_missinghits_plot.set_x(self.df_missing['Col'].to_numpy(), 'Column')
			col_vs_missinghits_plot.y_label = 'Missing Hits'

			col_vs_missinghits_plot.histo()

		# Plot missing hits vs row
		if row:
			row_vs_missinghits_file = File('../plots/colrow_vs_missing_hits/', 'plot_row_vs_missing_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
			row_vs_missinghits_plot = Plot(row_vs_missinghits_file, save, show)

			row_vs_missinghits_plot.set_x(self.df_missing['Row'].to_numpy(), 'Row')
			row_vs_missinghits_plot.y_label = 'Missing Hits'

			row_vs_missinghits_plot.histo()

	###### Plot data rate (hits per event) with missing hits #######

	def plot_events_vs_hits(self, save = True, show = True):

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
		events_vs_hits_file = File('../plots/events_vs_hits/', 'plot_events_vs_hits', self.file_name+'_'+self.fsm_clock, '.pdf')
		events_vs_hits_plot = Plot(events_vs_hits_file, save, show)

		events_vs_hits_plot.set_x(all_events, 'Event')
		events_vs_hits_plot.set_y(all_counts, 'Incoming Hits')
		events_vs_hits_plot.set_y2(missing_counts, 'Missing Hits')

		events_vs_hits_plot.two_in_one()

	###### Plot readout times ######################################

	def plot_rotime_vs_colrow(self, col = True, row = True, save = True, show = True):

		# Readout times in us
		df_ro = self.df[['Col','Row','Tread']].dropna()
		cols = df_ro['Col'].to_numpy()
		rows = df_ro['Row'].to_numpy()
		ro_times = df_ro['Tread'].to_numpy()
		ro_times = ro_times/1000

		# Plot readout time vs column
		if col:
			rotime_vs_col_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_vs_col', self.file_name+'_'+self.fsm_clock, '.pdf')
			rotime_vs_col_plot = Plot(rotime_vs_col_file, save, show)

			rotime_vs_col_plot.set_x(ro_times, 'Readout time (us)')
			rotime_vs_col_plot.set_y(cols, 'Column')

			rotime_vs_col_plot.sub_histo()

		# Plot readout time vs row
		if row:
			rotime_vs_row_file = File('../plots/readout_time_vs_colrow/', 'plot_rotime_vs_row', self.file_name+'_'+self.fsm_clock, '.pdf')
			rotime_vs_row_plot = Plot(rotime_vs_row_file, save, show)

			rotime_vs_row_plot.set_x(ro_times, 'Readout time (us)')
			rotime_vs_row_plot.set_y(rows, 'Row')

			rotime_vs_row_plot.sub_histo()

	###### Conclusion ##############################################

	def print_conclusion(self):

		self.helper.msg('ANALYSIS CONCLUDED...', type = 'head', to_log = True)

		if sum(self.missing) == 0:
			self.helper.msg('Success! All hits were accounted for.', to_log = True)
		elif abs(sum(self.missing)) == 1:
			self.helper.msg(str(sum(self.missing))+' hit still not accounted for.', to_log = True)
		else:
			self.helper.msg(str(sum(self.missing))+' hits still not accounted for.', to_log = True)

	###### Summary #################################################

	def print_summary(self):

		self.helper.msg('SUMMARY...\n', type='head', to_log = True)

		for item in self.helper.summary:
			self.helper.printing(item, to_log = True)

		self.helper.msg('', to_log = True)

	###### Close log file ##########################################

	def end_log(self):

		self.log.close()

		print('Log saved as:\n\t{}'.format(self.log_file.name))
		print('Log saved to:\n\t{}'.format(self.log_file.path))

		self.helper.hline()

	###### Full comparison encompassing all functions ##############

	# Run analysis
	def compare(self):
		# Initialising functions
		self.start_log()
		self.read_data()
		self.ts_correction()
		self.print_header()

		# Analysis
		self.check_in_out_hits()
		self.check_repeat_hits()
		self.check_multi_hits()
		self.check_missing_hits()
		self.check_tots()
		self.check_data_rate()

		# Concluding functions
		self.print_conclusion()
		self.print_summary()

	# Run all plots and save them
	def plots(self, save_me, show_me):
		#self.hitmap_hits(save = save_me, show = show_me)
		self.hitmap_missing_hits(save = save_me, show = show_me)
		self.plot_colrow_vs_hits(save = save_me, show = show_me)
		self.plot_colrow_vs_missing_hits(save = save_me, show = show_me)
		self.plot_events_vs_hits(save = save_me, show = show_me)
		self.plot_rotime_vs_colrow(save = save_me, show = show_me)

	def end(self):
		self.end_log()
