import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.signal as sig
import sys
import csv

dir = '/Users/sigridscherl/Documents/PhD/Projects/MightyPix1/Measurements/2023_09_KIT/TimeWalk/'
file = 'timewalk_parameter_test1.txt'

df = pd.DataFrame(pd.read_csv(dir+file, sep='\t'))

class PlotIt:

    def __init__(self,df,x,x_name,x_label, text = ''):
        self.df = df
        self.x = self.df[x]
        self.x_name = x_name
        self.x_label = x_label
        self.y2 = ''
        self.yerr = ''
        self.x_ticks = []
        self.text = text

    def fig_name(self):
        return 'TW_param_scan_'+self.x_name+'_vs_'+self.y_name+'.pdf'

    def plot(self):

        fig = plt.subplots(figsize = (6,4), layout='constrained')

        plt.plot(self.x,self.y, '.')
        if len(self.y2) > 0: plt.plot(self.x,self.y2,'.')
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        plt.title(self.text)
        if len(self.x_ticks) > 0: plt.xticks(self.x_ticks)

        plt.savefig(self.fig_name())
        plt.show()

    def error_plot(self):

        fig = plt.subplots(figsize = (6,4), layout='constrained')

        plt.errorbar(self.x,self.y,yerr = self.yerr, linestyle='', marker='.')
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        plt.title(self.text)
        if len(self.x_ticks) > 0: plt.xticks(self.x_ticks)

        plt.savefig(self.fig_name())
        plt.show()

    def plot_vs_mean_toa(self):
        self.y = self.df['mean toa>25']
        self.y_name = 'mean_toa'
        self.y_label = 'Mean ToA > 25 (ns)'
        self.plot()

    def plot_vs_std_toa(self):
        self.y = self.df['std toa>25']
        self.y_name = 'std_toa'
        self.y_label = 'Std ToA > 25 (ns)'
        self.plot()

    def plot_vs_max_toa(self):
        self.y = self.df['max toa']
        self.y_name = 'max_toa'
        self.y_label = 'Max ToA (ns)'
        self.plot()

    def plot_vs_mu(self):
        self.y = self.df['mu']
        self.y_name = 'mu'
        self.y_label = 'Mu (V)'
        self.plot()

    def plot_vs_sigma(self):
        self.y = self.df['sigma']
        self.y_name = 'sigma'
        self.y_label = 'Sigma (V)'
        self.plot()
    
    def plot_vs_mu_w_sigma(self):
        self.y = self.df['mu']
        self.yerr = self.df['sigma']
        self.y_name = 'mu_w_sigma'
        self.y_label = 'Mu (V)'
        self.error_plot()

    def plot_all(self):
        self.plot_vs_mean_toa()
        self.plot_vs_std_toa()
        self.plot_vs_max_toa()
        self.plot_vs_mu()
        self.plot_vs_sigma()
        self.plot_vs_mu_w_sigma()


###### Go through thpix from 572 to 592 in steps of 5, keeping blpix = 512, vn = 10 and vncomp = 10 ######
plot_thpix = PlotIt(df = df[(df['vn'] == 10) & (df['vncomp'] == 10)],
                    x = 'thpix',
                    x_name = 'thpix',
                    x_label = 'ThPix (DAC)',
                    text = 'BlPix = 512, VN = 10, VNComp = 10')
plot_thpix.x_ticks = [572,577,582,587,592]
#plot_thpix.plot_all()


###### Go through vncomp from 1 to 30 in steps of 5, keeping thpix = 582, blpix = 512 and vn = 10 ######
plot_vncomp = PlotIt(df = df[(df['vn'] == 10) & (df['thpix'] == 582)],
                    x = 'vncomp',
                    x_name = 'vncomp',
                    x_label = 'VNComp (DAC)',
                    text = 'ThPix = 582, BlPix = 512, VN = 10')
#plot_vncomp.plot_all()

###### Go through vn from 1 to 50 in steps of 5, keeping thpix = 582, blpix = 512 and vncomp = 1 ######
plot_vn = PlotIt(df = df[(df['vncomp'] == 1) & (df['thpix'] == 582)],
                 x = 'vn',
                 x_name = 'vn',
                 x_label = 'VN (DAC)',
                 text = 'ThPix = 582, BlPix = 512, VNComp = 1')
plot_vn.plot_all()