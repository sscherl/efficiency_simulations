#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Class to plot hitmap (heatmap) for a specific dataframe
#
#	Author: Sigrid Scherl
#
#	Created: February 2021
#

from ROOT import TFile, TTree, TH1D, TCanvas

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import sys
import csv
from matplotlib.lines import Line2D
import yaml
from obspy.imaging.cm import viridis_white_r

from classes.toolkit import Helper, LHCParameters

########## Load configs #############################################

with open("config.yaml", "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


###### Classes ##################################################

class Plot:

	def __init__(self, file, save, show):
		self.file = file
		self.save = save
		self.show = show

		# For histograms
		self.bins = 200
		self.y_min = ''
		self.y_max = ''

		self.y2_min = ''
		self.y2_max = ''

		self.x_min = ''
		self.x_max = ''

		self.x2_min = ''
		self.x2_max = ''

		# For hitmaps
		self.min_hits = 0
		self.max_hits = 0
		self.cutoff = False
		self.colour_scheme = 'viridis'

		self.plot_readout_time_limit = False

		self.helper = Helper()
		self.col = MyColours()
		self.lhc = LHCParameters()

		matplotlib.rcParams.update({'font.size': 19}) # Before: 14

	def set_x(self, value, label):
		self.x = value
		self.x_label = label

	def set_x2(self, value, label):
		self.x2 = value
		self.x2_label = label

	def set_y(self, value, label):
		self.y = value
		self.y_label = label

	def set_y2(self, value, label):
		self.y2 = value
		self.y2_label = label

	def set_pixels(self, cols, rows):
		self.cols = cols
		self.rows = rows

	def set_hits(self, min, max, cutoff = False):
		self.min_hits = min
		self.max_hits = max
		self.cutoff = cutoff

	def show_plot(self):
		plt.show()

	def save_plot(self):
		plt.savefig(self.file.whole)
		print('Plot saved as:\n\t{}'.format(self.file.name))
		print('Plot saved to:\n\t{}'.format(self.file.path))
		self.helper.hline()

	def end(self, tight=True):
		if tight: plt.tight_layout()
		if self.save: self.save_plot()
		if self.show: self.show_plot()
		plt.clf()
		plt.close('all')
		try: fig.clf()
		except: pass

	def histo(self):
		plt.hist(self.x, bins = self.bins)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label)
		if self.x_min != '' and self.x_max != '':
			plt.xlim(self.x_min,self.x_max)
		# plt.yscale('log')
		self.end()

	def histo_two_in_one(self):
		fig, sub = plt.subplots()

		hy, hx, _ = sub.hist(self.x, bins = self.bins, color = 'tab:blue')
		sub.set_xlabel(self.x_label)
		sub.set_ylabel(self.y_label, color = 'tab:blue')
		# sub.set_ylim(0,30)
		if self.y_max == '':
			sub.set_ylim(0, int(hy.max()*1.2))
		else:
			sub.set_ylim(0, self.y_max)
		

		sub2 = sub.twinx()
		sub2.hist(self.x2, bins = self.bins, color = 'tab:orange')
		sub2.set_ylabel(self.y2_label, color = 'tab:orange')
		# sub2.set_ylim(0,30)
		if self.y_max == '':
			sub2.set_ylim(0, int(hy.max()*1.2))
		else:
			sub2.set_ylim(0, self.y_max)

		self.end()

	def histo_two_in_one_wide(self):

		fig, sub = plt.subplots(figsize=(20,4.8))

		hy, hx, _ = sub.hist(self.x, bins = self.bins, color = 'tab:blue')
		sub.set_xlabel(self.x_label)
		sub.set_ylabel(self.y_label, color = 'tab:blue')
		# sub.set_ylim(0,30)
		if self.y_max == '':
			sub.set_ylim(0, int(hy.max()*1.2))
		else:
			sub.set_ylim(0, self.y_max)
		

		sub2 = sub.twinx()
		sub2.hist(self.x2, bins = self.bins, color = 'tab:orange')
		sub2.set_ylabel(self.y2_label, color = 'tab:orange')
		# sub2.set_ylim(0,30)
		if self.y_max == '':
			sub2.set_ylim(0, int(hy.max()*1.2))
		else:
			sub2.set_ylim(0, self.y_max)

		self.end()

	def two_in_one(self):
		fig, sub = plt.subplots()

		sub.plot(self.x, self.y, color = 'tab:blue')
		sub.set_xlabel(self.x_label)
		sub.set_ylabel(self.y_label, color = 'tab:blue')
		if self.y_max == '':
			sub.set_ylim(0, self.y.max()+1)
		else:
			sub.set_ylim(0, self.y_max)

		sub2 = sub.twinx()
		sub2.plot(self.x, self.y2, color = 'tab:orange')
		sub2.set_ylabel(self.y2_label, color = 'tab:orange')
		if self.y_max == '':
			sub2.set_ylim(0, self.y.max()+1)
		else:
			sub2.set_ylim(0, self.y_max)

		self.end()



	def sub_histo(self):
		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [2, 1]})
		fig.subplots_adjust(hspace=0)

		sub[0].plot(self.x, self.y, '.', color = 'tab:blue')
		sub[0].set_ylabel(self.y_label)

		if self.plot_readout_time_limit:
			sub[0].plot([self.lhc.readout_time_limit, self.lhc.readout_time_limit],[self.y_min,self.y_max], linewidth = 1, linestyle = '--', color = 'black')
			sub[0].set_ylim(self.y_min, self.y_max)
			sub[0].set_xlim(self.x_min, self.x_max)


		sub[1].hist(self.x, bins = self.bins, color = 'tab:blue')
		sub[1].set_xlabel(self.x_label)
		if self.x_max != '': sub[1].set_xlim(0,self.x_max)
		sub[1].set_ylabel('Count')
		sub[1].set_yscale('log')

		if self.plot_readout_time_limit:
			sub[1].plot([self.lhc.readout_time_limit, self.lhc.readout_time_limit],[self.y2_min, self.y2_max], linewidth = 1, linestyle = '--', color = 'black')
			sub[1].set_ylim(self.y2_min, self.y2_max)
			sub[1].set_xlim(self.x_min, self.x_max)

		self.end()



	def sub_histo2(self):
		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [2, 1]})
		fig.subplots_adjust(hspace=0)

		sub[0].plot(self.x, self.y, '.', color = 'tab:blue')
		sub[0].plot(self.x2, self.y, 'v', color = 'tab:orange')
		sub[0].set_ylabel(self.y_label)

		sub[0].legend(handles = [Line2D([0],[0], lw = 0, color = 'tab:blue', marker = '.', label = 'Data'), Line2D([0],[0], lw = 0, color = 'tab:orange', marker = 'v', label = 'Mean')])

		sub[1].hist(self.x, bins = self.bins, color = 'tab:blue')
		sub[1].set_xlabel(self.x_label)
		sub[1].set_ylabel('Count')

		self.end()



	def sub_histo3(self):
		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [2, 1]})
		fig.subplots_adjust(hspace=0)

		sub[0].plot(self.x, self.y, '.', color = 'tab:blue')
		sub[0].plot(self.x2, self.y2, '.', color = 'tab:orange')
		sub[0].set_ylabel(self.y_label)

		sub[0].legend(handles = [Line2D([0],[0], lw = 0, color = 'tab:blue', marker = '.', label = 'Read during simulation'), Line2D([0],[0], lw = 0, color = 'tab:orange', marker = '.', label = 'Read after simulation')], loc = 'upper right')

		# sub[1].hist(self.x, bins = self.bins, color = 'tab:blue')
		sub[1].hist(self.x, bins = self.bins, color = 'tab:blue')
		sub[1].hist(self.x2, bins = self.bins, color = 'tab:orange')
		sub[1].set_xlabel(self.x_label)
		sub[1].set_ylabel('Count')
		sub[1].set_yscale('log')

		self.end()



	def sub_histo4(self):
		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [2, 1]})
		fig.subplots_adjust(hspace=0)

		sub[0].plot(self.x, self.y, '.', color = 'tab:blue')
		sub[0].plot(self.x2, self.y2, '.', color = 'tab:orange')
		sub[0].set_ylabel(self.y_label)
		# sub[0].set_xscale('log')

		sub[0].legend(handles = [Line2D([0],[0], lw = 0, color = 'tab:blue', marker = '.', label = 'Below max'), Line2D([0],[0], lw = 0, color = 'tab:orange', marker = '.', label = 'Above max')], loc = 'upper right')

		bins = np.histogram(np.hstack((self.x,self.x2)), bins = self.bins)[1]

		# sub[1].hist(self.x, bins = self.bins, color = 'tab:blue')
		sub[1].hist(self.x, bins = bins, color = 'tab:blue')
		sub[1].hist(self.x2, bins = bins, color = 'tab:orange')
		sub[1].set_xlabel(self.x_label)
		sub[1].set_ylabel('Count')
		sub[1].set_yscale('log')
		# sub[1].set_xscale('log')

		self.end()



	def scatter(self):
		plt.plot(self.x, self.y, color = 'tab:blue', marker = '.', linewidth = 0)
		plt.xlabel(self.x_label)
		plt.ylabel(self.y_label, color = 'tab:blue')
		if self.y_max == '':
			if self.y_min == '':
				plt.ylim(0, self.y.max()+1)
			else:
				plt.ylim(self.y_min, self.y.max()+1)
		else:
			plt.ylim(0, self.y_max)

		if self.x_min != '':
			plt.xlim(self.x_min, self.x.max())

		self.end()



	def lines(self):
		pass



	def hitmap(self, df):
		hit_grid  = np.zeros((self.rows, self.cols))

		# Fill hit grid
		for i in range(len(df)):
			curr_col = int(df.loc[i, 'Col'])
			curr_row = int(df.loc[i, 'Row'])
			hit_grid[-curr_row][curr_col] += 1

		# Plot grid in certain range:
		y_dim, x_dim = hit_grid.shape

		fig = plt.figure(figsize=(6,5),dpi=100)

		if self.min_hits == self.max_hits:
			plt.imshow(hit_grid, cmap=self.colour_scheme, interpolation="nearest", aspect="auto", extent = [0, x_dim, 0, y_dim])
		elif self.min_hits < self.max_hits:
			if self.max_hits >= np.amax(hit_grid) or self.cutoff:
				plt.imshow(hit_grid, cmap=self.colour_scheme, vmin = self.min_hits, vmax = self.max_hits, interpolation="nearest", aspect="auto", extent = [0, x_dim, 0, y_dim])
			else:
				self.helper.raise_fail("Error! Max number of plotted hits smaller than detected hits. Choose higher max, enable cutoff, or set range to automatic.")
		else:
			self.helper.raise_fail("Error! Max number of plotted hits smaller than min. Choose higher max or set range to automatic.")

		plt.xlabel('Pixel columns')
		plt.ylabel('Pixel rows')
		
		cbar = plt.colorbar()
		cbar.set_label('Hits')

		self.end()


	# Heat map with one histo below
	def hitmap_sub_histo(self, ybins, vmax='', title='', ticks = []):

		matplotlib.rcParams.update({'font.size': 14})

		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(5, 5))
		#fig.subplots_adjust(hspace=0)

		xbins = self.bins

		if vmax == '':
			hitmap = sub[0].hist2d(self.x, self.y, bins=[xbins,ybins], cmin = 1, vmin = 1, cmap= matplotlib.colormaps.get_cmap('viridis_r')) # viridis_white_r
		else:
			hitmap = sub[0].hist2d(self.x, self.y, bins=[xbins,ybins], cmin = 1, vmin = 1, vmax = vmax, cmap= matplotlib.colormaps.get_cmap('viridis_r')) #viridis_white_r
		#sub[0].plot(self.x, self.y, '.', color = 'tab:blue')
		sub[0].set_ylabel(self.y_label)
		if title != '': sub[0].set_title(title)


		if self.plot_readout_time_limit:
			sub[0].plot([self.lhc.readout_time_limit, self.lhc.readout_time_limit],[self.y_min,self.y_max], linewidth = 1, linestyle = '--', color = 'black')
		
		sub[0].set_ylim(self.y_min, self.y_max)
		sub[0].set_xlim(self.x_min, self.x_max)


		sub[1].hist(self.x, bins = xbins, color = 'tab:blue')
		sub[1].set_xlabel(self.x_label)
		if self.x_max != '': sub[1].set_xlim(0,self.x_max)
		sub[1].set_ylabel('Count')
		sub[1].set_yscale('log')

		if self.plot_readout_time_limit:
			sub[1].plot([self.lhc.readout_time_limit, self.lhc.readout_time_limit],[self.y2_min, self.y2_max], linewidth = 1, linestyle = '--', color = 'black')
		
		sub[1].set_ylim(self.y2_min, self.y2_max)
		sub[1].set_xlim(self.x_min, self.x_max)

		bottom = 0.12
		left = 0.15
		right = 0.75
		top = 0.92
		spacing = 0.02

		plt.subplots_adjust(bottom=bottom, left=left, right= right, top=top, hspace = 3*spacing)

		cax = plt.axes([right+spacing, bottom, 0.05, top-bottom]) # Left, bottom, width, height

		if ticks == []:
			fig.colorbar(mappable=hitmap[3], cax = cax, label = 'Counts')
		else:
			fig.colorbar(mappable=hitmap[3], cax = cax, ticks = ticks, label = 'Counts')

		self.end(tight=False)


	# Heat map with one histo below
	def hitmap_sub_histo_wide(self, ybins, vmax='', title='', ticks = []):

		matplotlib.rcParams.update({'font.size': 14})

		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(10, 10))
		#fig.subplots_adjust(hspace=0)

		xbins = self.bins

		if vmax == '':
			hitmap = sub[0].hist2d(self.x, self.y, bins=[xbins,ybins], cmin = 0, vmin = 0, cmap= viridis_white_r)
		else:
			hitmap = sub[0].hist2d(self.x, self.y, bins=[xbins,ybins], cmin = 0, vmin = 0, vmax = vmax, cmap= viridis_white_r)
		#sub[0].plot(self.x, self.y, '.', color = 'tab:blue')
		sub[0].set_ylabel(self.y_label)
		if title != '': sub[0].set_title(title)


		if self.plot_readout_time_limit:
			sub[0].plot([self.lhc.readout_time_limit, self.lhc.readout_time_limit],[self.y_min,self.y_max], linewidth = 1, linestyle = '--', color = 'black')
		
		sub[0].set_ylim(self.y_min, self.y_max)
		sub[0].set_xlim(self.x_min, self.x_max)


		sub[1].hist(self.x, bins = xbins, color = 'tab:blue')
		sub[1].set_xlabel(self.x_label)
		if self.x_max != '': sub[1].set_xlim(0,self.x_max)
		sub[1].set_ylabel('Count')
		sub[1].set_yscale('log')

		if self.plot_readout_time_limit:
			sub[1].plot([self.lhc.readout_time_limit, self.lhc.readout_time_limit],[self.y2_min, self.y2_max], linewidth = 1, linestyle = '--', color = 'black')
		
		sub[1].set_ylim(self.y2_min, self.y2_max)
		sub[1].set_xlim(self.x_min, self.x_max)

		bottom = 0.12
		left = 0.15
		right = 0.75
		top = 0.92
		spacing = 0.02

		plt.subplots_adjust(bottom=bottom, left=left, right= right, top=top, hspace = 3*spacing)

		cax = plt.axes([right+spacing, bottom, 0.05, top-bottom]) # Left, bottom, width, height

		if ticks == []:
			fig.colorbar(mappable=hitmap[3], cax = cax, label = 'Counts')
		else:
			fig.colorbar(mappable=hitmap[3], cax = cax, ticks = ticks, label = 'Counts')

		self.end(tight=False)

	# NOT USED: Long heat map with one histo below
	def hitmap_sub_histo_long(self, ybins, vmax = '', col_max = '', row_max = '', ticks = []):
		matplotlib.rcParams.update({'font.size': 14})

		x = self.x
		y = self.y
		xbins = self.bins

		# definitions for the axes
		left, bottom, width, height, small_height = 0.2, 0.07, 0.6, 0.81, 0.1
		spacing = 0.01

		# start with a rectangular figure
		fig = plt.figure(figsize=(5, 8))

		ax_map = plt.axes([left, bottom+small_height+spacing, width, height])
		ax_map.tick_params(direction='inout', right=True)

		ax_rotime = plt.axes([left, bottom, width, small_height])
		ax_rotime.tick_params(direction='inout', labelleft=True, labelbottom = False)

		# the histograms
		if vmax == '':
			mappable = ax_map.hist2d(x, y, bins=[xbins,ybins], cmin = 1, vmin = 1, cmap= matplotlib.colormaps.get_cmap('viridis_r')) # For 40 MHz: vmin = 0, vmax = 80
		else: 
			mappable = ax_map.hist2d(x, y, bins=[xbins,ybins], cmin = 1, vmin = 1, vmax = vmax, cmap= matplotlib.colormaps.get_cmap('viridis_r'))

		ax_rotime.hist(x, bins=xbins, orientation='vertical')
		if col_max != '':
			ax_rotime.set_ylim(ymin=0, ymax=col_max) # For 40 MHz: ymin = 0, ymax = 16000

		# align axes
		ax_rotime.set_xlim(ax_map.get_xlim())

		cax = plt.axes([left + width + spacing, bottom, 0.05, height + spacing + small_height])

		if ticks == []:
			fig.colorbar(mappable[3], ax = ax_map, cax = cax, label= 'Counts')
		else:
			fig.colorbar(mappable[3], ax = ax_map, cax = cax, label= 'Counts', ticks = ticks)

		# axis labels
		ax_map.set_ylabel('Rows')
		ax_map.set_xlabel('Columns')
		ax_rotime.set_ylabel('Counts')
		#ax_rows.set_xscale('log')
		#ax_rows.set_title('Missing Hits per Row')

		self.end(tight=False)
	

	# Heat map with one histo on top and one to thr left
	def hitmap_side_histos(self,x,y,xbins, ybins, vmax = '', col_max = '', row_max = '', ticks = [], cmap = viridis_white_r, cmin = 0, vmin = 0):

		matplotlib.rcParams.update({'font.size': 14})

		# definitions for the axes
		left, bottom, width, height, small_height, small_width = 0.2, 0.07, 0.4, 0.81, 0.1, 0.19
		spacing = 0.01

		# start with a rectangular figure
		fig = plt.figure(figsize=(5, 8))

		ax_map = plt.axes([left, bottom, width, height])
		ax_map.tick_params(direction='inout', right=True)

		ax_rows = plt.axes([left + width + spacing * 1.9, bottom, small_width, height])
		ax_rows.tick_params(direction='inout', labelleft=False)

		ax_cols = plt.axes([left, bottom + height + spacing, width, small_height])
		ax_cols.tick_params(direction='inout', labelleft=True, labelbottom= False)

		# the histograms
		if vmax == '':
			mappable = ax_map.hist2d(x, y, bins=[xbins,ybins], cmin = cmin, vmin = vmin, cmap= cmap) # For 40 MHz: vmin = 0, vmax = 80
		else:
			mappable = ax_map.hist2d(x, y, bins=[xbins,ybins], cmin = cmin, vmin = vmin, vmax = vmax, cmap= cmap)

		ax_rows.hist(y, bins=ybins, orientation='horizontal')
		if row_max != '':
			ax_rows.set_xlim(xmin=0, xmax=row_max) # For 40 MHz: xmin = 0, xmax = 1500
		ax_cols.hist(x, bins=xbins, orientation='vertical')
		if col_max != '':
			ax_cols.set_ylim(ymin=0, ymax=col_max) # For 40 MHz: ymin = 0, ymax = 16000

		# align axes
		ax_rows.set_ylim(ax_map.get_ylim())
		ax_cols.set_xlim(ax_map.get_xlim())

		cax = plt.axes([left + width + 2*1.9 * spacing + small_width, bottom, 0.05, height])

		if ticks == []:
			fig.colorbar(mappable[3], ax = ax_map, cax = cax, label= 'Counts')
		else:
			fig.colorbar(mappable[3], ax = ax_map, cax = cax, label= 'Counts', ticks = ticks)

		# axis labels
		ax_map.set_ylabel('Rows')
		ax_map.set_xlabel('Columns')
		ax_rows.set_xlabel('Counts')
		#ax_rows.set_xscale('log')
		#ax_rows.set_title('Missing Hits per Row')
		ax_cols.set_ylabel('Counts')

		self.end(tight=False)



class RootHisto:

	def __init__(self, dir, plot_name, plot_title, x_axis_title, x_values, bins, min, max):
		self.directory = dir
		self.plot_name = plot_name
		self.plot_title = plot_title
		self.x_axis_title = x_axis_title
		self.x_values = x_values
		self.bins = bins
		self.min = min
		self.max = max

		self.create_histo()

	def create_histo(self):

		self.canv = TCanvas('Canvas', 'Canvas')
		self.histo = TH1D('Stats', self.plot_title, self.bins, self.min, self.max)

		for value in self.x_values:
			self.histo.Fill(value)

		self.histo.SetStats(0)
		self.histo.GetXaxis().SetTitle(self.x_axis_title)
		self.histo.GetYaxis().SetTitle('Counts')
		self.histo.Draw()

		self.canv.Update()
		self.canv.Print('../plots/'+self.directory+'/'+self.plot_name+'.pdf')
		RootFile = TFile.Open('../plots/'+self.directory+'/'+self.plot_name+'.root', 'RECREATE')
		RootFile.WriteObject(self.histo, self.plot_title)

		del self.canv
		del self.histo


class MyColours:

	def __init__(self):
		self.dark_blue		= '#6C8EBF'
		self.mid_blue		= '#A3BBDE'
		self.light_blue		= '#DAE8FC'

		self.dark_green		= '#82B366'
		self.mid_green		= '#ABCD9D'
		self.light_green	= '#D5E8D4'
		
		self.dark_yellow 	= '#E9CA2C'
		self.mid_yellow		= '#F4DE7C'
		self.light_yellow	= '#FFF2CC'

		self.dark_orange	= '#D79B00'
		self.mid_orange		= '#EBC066'
		self.light_orange	= '#FFE6CC'

		self.dark_red		= '#B85450'
		self.mid_red		= '#D8918E'
		self.light_red		= '#F8CECC'

		self.dark_purple	= '#9673A6'
		self.mid_purple		= '#BBA4C6'
		self.light_purple	= '#E1D5E7'

		self.dark_grey		= '#999999'
		self.mid_grey		= '#C7C7C7'
		self.light_grey		= '#F5F5F5'

		self.black			= '#000000'
		self.grey			= '#7F7F7F'
		self.white			= '#FFFFFF'