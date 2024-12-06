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
import matplotlib.pyplot as plt
import sys
import csv

from classes.toolkit import Helper

###### Classes ##################################################

class Plot:

	def __init__(self, file, save, show):
		self.file = file
		self.save = save
		self.show = show

		# For histograms
		self.bins = 200

		# For hitmaps
		self.min_hits = 0
		self.max_hits = 0
		self.cutoff = False
		self.colour_scheme = 'viridis'

		self.helper = Helper()

	def set_x(self, value, label):
		self.x = value
		self.x_label = label

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

	def end(self):
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
		self.end()

	def two_in_one(self):
		fig, sub = plt.subplots()

		sub.plot(self.x, self.y, color = 'C0')
		sub.set_xlabel(self.x_label)
		sub.set_ylabel(self.y_label, color = 'C0')
		sub.set_ylim(0, self.y.max()+1)

		sub2 = sub.twinx()
		sub2.plot(self.x, self.y2, color = 'C1')
		sub2.set_ylabel(self.y2_label, color = 'C1')
		sub2.set_ylim(0, self.y.max()+1)

		self.end()

	def sub_histo(self):
		fig, sub = plt.subplots(2,1, sharex = True, gridspec_kw={'height_ratios': [2, 1]})
		fig.subplots_adjust(hspace=0)

		sub[0].plot(self.x, self.y, '.')
		sub[0].set_ylabel(self.y_label)

		sub[1].hist(self.x, bins = self.bins)
		sub[1].set_xlabel(self.x_label)
		sub[1].set_ylabel('Count')

		self.end()

	def scatter(self):
		pass

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

		plt.colorbar()
		plt.xlabel('Pixel columns')
		plt.ylabel('Pixel rows')

		self.end()


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

		self.histo.GetXaxis().SetTitle(self.x_axis_title)
		self.histo.GetYaxis().SetTitle('Counts')
		self.histo.Draw()

		self.canv.Update()
		self.canv.Print('../plots_v2/'+self.directory+'/'+self.plot_name+'.pdf')

		del self.canv
		del self.histo