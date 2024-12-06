#
#	LHCb Verification Framework for the MightyPix
#
#	Task: Quickly calculate Clopper-Pearson interval for manually entered values
#
#	Author: Sigrid Scherl
#
#	Created: March 2023
#

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import csv
from statistics import mean, stdev
from scipy.stats import beta


def cloppear(total, detected):
    k = detected
    n = total
    alpha = 0.05
    low, high = beta.ppf([alpha/2, 1 - alpha/2], [k, k + 1], [n - k + 1, n - k])
    
    if np.isnan(low):   low = 0
    if np.isnan(high):  high = 1
    
    return low*100, high*100


def result(total, detected):
    eff = (detected/total)*100

    low, high = cloppear(total, detected)

    print ('='*100)
    print('Result for ', total, ' hits, out of which ', total - detected, ' are missing.')
    print ('-'*100)
    print('95%-CPI    = ', low, ', ', high)
    print ('-'*100)
    print('Efficiency = ', eff, ' + ', high - eff, ' - ', eff - low)
    print ('='*100)

#####################################################################

# 40 MHZ/cm^2 
# tot = 420715
# miss = 190831 + 12768

# 17 MHZ/cm^2

tot = 1183
miss = 14
det = tot - miss

result(total = tot, detected = det)