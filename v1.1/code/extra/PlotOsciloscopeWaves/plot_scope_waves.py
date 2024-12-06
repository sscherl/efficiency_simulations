import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scipy.signal as sig
import sys
import csv

matplotlib.rcParams.update({'font.size': 12})

dir = '/Users/sigridscherl/Documents/PhD/Projects/MightyPix1/Measurements/2023_09_KIT/AmpOutHitBus/2023_09_13/'
ch1 = ['230523_235656004_Ch1.csv', 'Injection']
ch2 = ['230523_235656004_Ch2.csv', 'AmpOut']
ch3 = ['230523_235656004_Ch4.csv', 'HitBus']


def plot_all():

    fig = plt.subplots(figsize = (6,4), layout='constrained')

    for ch in [ch1, ch2, ch3]:

        df = pd.DataFrame(pd.read_csv(dir+ch[0], header = None))

        x, y = df[3] * 1e6, df[4]
        min_x = min(x)

        plt.plot(x-min_x,y, label = ch[1])

    plt.xlabel('Time (µs)')
    plt.ylabel('Voltage (V)')
    plt.legend()
    plt.savefig('inj_ampout_hitbus_-50HV.pdf')
    plt.show()



def plot_ampout():

    fig = plt.subplots(figsize = (6,4), layout='constrained')

    for ch in [ch2]:

        df = pd.DataFrame(pd.read_csv(dir+ch[0], header = None))

        x, y = df[3] * 1e6, df[4]
        min_x = min(x)

        plt.plot(x-min_x,y, label = ch[1], color = 'tab:orange')

    plt.arrow(9.2, 0.06, -1.2, 0.12, head_width = 0.014, head_length = 0.08, color = 'black') #length_includes_head = True,
    plt.arrow(9.5, 0.06, 0.9, 0.03, head_width = 0.014, head_length = 0.08, color = 'black')
    plt.Circle((2, 2), 3, color='b', fill=False)

    plt.text(8.5, 0.03, 'Feedback')

    plt.xlim(5,15)
    plt.ylim(-0.04,0.4)

    plt.xlabel('Time (µs)')
    plt.ylabel('Voltage (V)')
    plt.legend()
    plt.savefig('ampout_-50HV.pdf')
    plt.show()

plot_ampout()