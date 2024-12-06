#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Extract data from .root file and plot histograms
#
#	Author: Sigrid Scherl
#
#	Created: June 2022
#

import ROOT
from ROOT import TFile, TTree, TH1D, TCanvas

import csv
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
import sys
import time

from toolkit import ProgressBar
from simulation_data_plot import RootHisto

# For all histos
folder = 'root_data_histograms'
data_type = 'LHCb Simulation Data'

# Read file and get tree
file = TFile('../data_v2/raw_data/MagUp_bs2phiphi_violaine_energytime_oscar_files_1p5e34.root_Sorted_Sensors.root')
tree = file.Get('SensorHits')
entries = tree.GetEntries()

# Initialise variables for leaves
run_num, evt_num, layer, x, y, z, x_symloc, y_symloc, x0, y0  = [], [], [], [], [], [], [], [], [], []

# Get tree content
print('Reading tree...')
progress_bar = ProgressBar(entries)

for i in range(0, entries):
    tree.GetEntry(i)
    run_num.append(tree.run_num)      # collection of simulated events
    evt_num.append(tree.evt_num)      # identifier of event within a given run
    layer.append(tree.layer)          # layer number of hit, only x-aligned layers included, 1-6 in ascending order with z position
    x.append(tree.x)                  # x position of hit in SciFi in global detector coordinate system
    y.append(tree.y)                  # y position of hit in SciFi in global detector coordinate system
    z.append(tree.z)                  # z position of hit in SciFi in global detector coordinate system
    x_symloc.append(tree.x_symloc)    # symmetrised local sensor coordinates, described above in preamble
    y_symloc.append(tree.y_symloc)    # symmetrical local sensor coordinates, described above in preamble
    x0.append(tree.x0)                # sensor x centre position
    y0.append(tree.y0)                # sensor y centre position

    progress_bar.PrintBar(i)

print('Creating histograms...')
histo_run_num  = RootHisto(directory = folder, plot_name = 'histo_run_num', plot_title = data_type, x_axis_title = 'Run Number',   x_values = run_num, bins = 100, min = 4990,  max = 5110)
histo_evt_num  = RootHisto(directory = folder, plot_name = 'histo_evt_num', plot_title = data_type, x_axis_title = 'Event Number', x_values = evt_num, bins = 60,  min = 0,     max = 6   )
histo_layer    = RootHisto(directory = folder, plot_name = 'histo_layer',   plot_title = data_type, x_axis_title = 'Layer',        x_values = layer,   bins = 70,  min = 0,     max = 7   )
histo_x        = RootHisto(directory = folder, plot_name = 'histo_x',       plot_title = data_type, x_axis_title = 'X (mm)',       x_values = x,       bins = 100, min = -3000, max = 3000)
histo_y        = RootHisto(directory = folder, plot_name = 'histo_y',       plot_title = data_type, x_axis_title = 'Y (mm)',       x_values = y,       bins = 100, min = -700,  max = 700 )
histo_z        = RootHisto(directory = folder, plot_name = 'histo_z',       plot_title = data_type, x_axis_title = 'Z (mm)',       x_values = z,       bins = 100, min = 7700,  max = 9500)
histo_x_symloc = RootHisto(directory = folder, plot_name = 'histo_x_symloc',plot_title = data_type, x_axis_title = 'X_symloc (mm)',x_values = x_symloc,bins = 100, min = -15,   max = 15  )
histo_y_symloc = RootHisto(directory = folder, plot_name = 'histo_y_symloc',plot_title = data_type, x_axis_title = 'Y_symloc (mm)',x_values = y_symloc,bins = 100, min = -15,   max = 15  )
histo_x0       = RootHisto(directory = folder, plot_name = 'histo_x0',      plot_title = data_type, x_axis_title = 'X0 (mm)',      x_values = x0,      bins = 100, min = 0,     max = 2500)
histo_y0       = RootHisto(directory = folder, plot_name = 'histo_y0',      plot_title = data_type, x_axis_title = 'Y0 (mm)',      x_values = y0,      bins = 100, min = 0,     max = 700 )
