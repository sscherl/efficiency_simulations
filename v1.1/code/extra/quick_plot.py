#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Quickly plot various things
#
#	Author: Sigrid Scherl
#
#	Created: February 2021
#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import csv

###### EFFICIENCY PLOTS - Hit rate vs missing hit percentage #######

###### Inputs #####################################################

tot = 5
clusters = False

###### Data #####################################################
theo_full_hit_rates = [1.7, 3.4, 5.1, 6.8]

### Missing hits at Varying Data Rates ###

# Readout FSM = 40 MHz, ToT = 2 us
if tot == 2 and not clusters:
    theo_hit_rates = [0.37, 0.75, 1.12, 1.49]
    lhcb_hit_rates = [0.39, 0.77, 1.18, 1.58]
    rnd_hit_rates = [0.36, 0.73, 1.10, 1.50]

    lhcb_missing_hits_prc = [0.77, 4.52, 9.85, 14.87]
    lhcb_missing_hits_prc_err = [0.04, 0.14, 0.24, 0.31]
    rnd_missing_hits_prc = [0.37, 4.26, 9.34, 14.86]
    rnd_missing_hits_prc_err = [0.02, 0.13, 0.23, 0.32]

# Readout FSM = 40 MHz, ToT = 5 us
elif tot == 5 and not clusters:
    theo_hit_rates = [0.37, 0.75, 1.12, 1.49]
    lhcb_hit_rates = [0.39, 0.77, 1.18, 1.58]
    rnd_hit_rates = [0.37, 0.73, 1.12, 1.49]

    lhcb_missing_hits_prc = [1.37, 5.25, 10.59, 15.42]
    lhcb_missing_hits_prc_err = [0.06, 0.16, 0.26, 0.32]
    rnd_missing_hits_prc = [0.71, 4.45, 10.61, 14.44]
    rnd_missing_hits_prc_err = [0.03, 0.14, 0.26, 0.31]

# Readout FSM = 40 MHz, ToT = 2 us, with clusters
elif tot == 2 and clusters:
    theo_hit_rates = [0.39, 0.78, 1.17, 1.57]
    lhcb_hit_rates = [0.41, 0.81, 1.23, 1.66]
    rnd_hit_rates = [0.40, 0.78, 1.18, 1.57]

    lhcb_missing_hits_prc = [0.74, 5.01, 10.56, 15.86]
    lhcb_missing_hits_prc_err = [0.03, 0.15, 0.25, 0.32]
    rnd_missing_hits_prc = [0.25, 5.53, 10.40, 15.55]
    rnd_missing_hits_prc_err = [0.02, 0.17, 0.25, 0.32]


else: print('ERROR: Invalid inputs!')

###### Main #####################################################

hit_rates_label = 'Hit rate'
hit_rates_mhz_label = 'Hit rate (MHz/cm^2)'
missing_hits_label_prc = 'Missing hits (%)'

x0 = lhcb_hit_rates
y0 = lhcb_missing_hits_prc
y0err = lhcb_missing_hits_prc_err
label0 = 'LHCb simulation'

x1 = rnd_hit_rates
y1 = rnd_missing_hits_prc
y1err = rnd_missing_hits_prc_err
label1 = 'Random'

x_label = hit_rates_label
y_label = missing_hits_label_prc

if clusters: file_name = 'plot_missing_per_hit_rate_tot'+str(tot)+'us_40mhz_var_rate_clusters.pdf'
else: file_name = 'plot_missing_per_hit_rate_tot'+str(tot)+'us_40mhz_var_rate.pdf'

file_path = '../plots/hit_rate_vs_missing_hits/'

# Fits
m0, b0 = np.polyfit(x0, y0, 1)
m1, b1 = np.polyfit(x1, y1, 1)

l0 = m0 * np.array(x0) + b0
l1 = m1 * np.array(x1) + b1

plt.errorbar(x0, y0, yerr = y0err, color = 'C0', fmt = '.', label = label0)
plt.errorbar(x1, y1, yerr = y1err, color = 'C1', fmt = '.', label = label1)
plt.plot(x0, l0, color = 'C0')
plt.plot(x1, l1, color = 'C1')
# plt.plot(x2, y2, color = 'C2')

plt.xlabel(x_label)
plt.ylabel(y_label)
# plt.set_ylim(0, y.max())

plt.legend()

plt.savefig(file_path+file_name)
print('Plot saved as:\n\t{}'.format(file_name))
print('Plot saved to:\n\t{}'.format(file_path))
plt.show()
