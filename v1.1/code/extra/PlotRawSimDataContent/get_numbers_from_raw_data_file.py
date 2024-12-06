import csv
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib
from mpl_toolkits import mplot3d
import numpy as np
import pandas as pd
from pathlib import Path
import random
from random import randint
import sys
import time
import math
from scipy.optimize import curve_fit



class ExtractSimData:


    def __init__(self, sensor_id):

        self.raw_data_file = '../../../data/raw_data/MagUp_bs2phiphi_violaine_energytime_oscar_files_1p5e34_Sorted_Sensors.csv'
        self.plot_dir = '../../../plots/root_data_plots/'

        self.sensor = sensor_id
        self.hottest_sensor = 739

        matplotlib.rcParams.update({'font.size': 18})


    def get_raw_data(self):
        data = pd.read_csv(self.raw_data_file)
        self.df_raw = pd.DataFrame(data)

    def data_for_sensor(self):
        if self.sensor == '':
            self.df_sensor = self.df_raw
        else:
            self.df_sensor = self.df_raw[self.df_raw.Sensor == self.sensor]

    def distance(self,x,y):
        return np.sqrt(x**2 + y**2)

    def get_closest_sensor(self):
        self.df_raw.Distance = self.distance(self.df_raw.X0, self.df_raw.Y0)
        min = self.df_raw.Distance.min()
        self.closest_sensor = self.df_raw[self.df_raw.Distance == min].Sensor.iloc[0]
        print('Hottest sensor is ', self.closest_sensor)

    def get_farthest_sensor(self):
        self.df_raw.Distance = self.distance(self.df_raw.X0, self.df_raw.Y0)
        max = self.df_raw.Distance.max()
        self.farthest_sensor = self.df_raw[self.df_raw.Distance == max].Sensor.iloc[0]
        print('Coldest sensor is ', self.farthest_sensor)

    def get_hits(self, layer, quadrant):
        return len(self.df_sensor[(self.df_sensor.Layer == layer) & (self.df_sensor.Quadrant == quadrant)])
    
    def get_rate(self,hits):
        if self.sensor == '':
            rate = round(hits/ (2964 * 500 * 25e-9 * 4 * 1e6),2)
        else:
            rate = round(hits/ (500 * 25e-9 * 4 * 1e6),4)
        return rate


    def lq_data(self):
        layers = []
        quadrants = []
        hits = []
        rates = []
        for l in range(1,7):
            for q in range (1,5):
                layers.append(l)
                quadrants.append(q)
                hits.append(self.get_hits(l,q))
                rates.append(self.get_rate(self.get_hits(l,q)))

        self.df_lq = pd.DataFrame(columns=['Layer','Quadrant','Hits','Rate'])
        self.df_lq.Layer = layers
        self.df_lq.Quadrant = quadrants
        self.df_lq.Hits = hits
        self.df_lq.Rate = rates


    def mean_per_layer(self):

        layers = []
        mean_hits = []
        std_hits = []
        mean_rates = []
        std_rates = []

        for l in range(1,7):
            this_df = data.df_lq[data.df_lq.Layer == l]
            df_mean = this_df.mean()
            df_std = this_df.std()

            layers.append(l)
            mean_hits.append(round(df_mean.Hits,2))
            std_hits.append(round(df_std.Hits,2))
            mean_rates.append(round(df_mean.Rate,3))
            std_rates.append(round(df_std.Rate,3))

        self.df_mean_per_layer = pd.DataFrame(columns=['Layer','MeanHits', 'StdHits', 'MeanRate', 'StdRate'])
        self.df_mean_per_layer.Layer = layers
        self.df_mean_per_layer.MeanHits = mean_hits
        self.df_mean_per_layer.StdHits = std_hits
        self.df_mean_per_layer.MeanRate = mean_rates
        self.df_mean_per_layer.StdRate = std_rates


    def mean_per_quadrant(self):

        quadrants = []
        mean_hits = []
        std_hits = []
        mean_rates = []
        std_rates = []

        for q in range(1,5):
            this_df = data.df_lq[data.df_lq.Quadrant == q]
            df_mean = this_df.mean()
            df_std = this_df.std()

            quadrants.append(q)
            mean_hits.append(round(df_mean.Hits,2))
            std_hits.append(round(df_std.Hits,2))
            mean_rates.append(round(df_mean.Rate,3))
            std_rates.append(round(df_std.Rate,3))

        self.df_mean_per_quadrant = pd.DataFrame(columns=['Quadrant','MeanHits', 'StdHits', 'MeanRate', 'StdRate'])
        self.df_mean_per_quadrant.Quadrant = quadrants
        self.df_mean_per_quadrant.MeanHits = mean_hits
        self.df_mean_per_quadrant.StdHits = std_hits
        self.df_mean_per_quadrant.MeanRate = mean_rates
        self.df_mean_per_quadrant.StdRate = std_rates

    def mean_overall(self):

        df_mean = self.df_lq.mean()
        df_std = self.df_lq.std()

        self.df_mean_overall = pd.DataFrame(columns=['MeanHits', 'StdHits', 'MeanRate', 'StdRate'])

        self.df_mean_overall.MeanHits = [round(df_mean.Hits,2)]
        self.df_mean_overall.StdHits = [round(df_std.Hits,2)]
        self.df_mean_overall.MeanRate = [round(df_mean.Rate,3)]
        self.df_mean_overall.StdRate = [round(df_std.Rate,3)]


    def linearFunc(self, x,intercept,slope):
        y = intercept + slope * x
        return y
    

    def perform_fit(self,whatplot):

        layer_list = [1,2,3,4,5,6]
        quadrant_list = [1,2,3,4]

        if whatplot == 'Layer':
            xlist = layer_list
            mean_rates = self.df_mean_per_layer.MeanRate
            std_rates = self.df_mean_per_layer.StdRate
        elif whatplot == 'Quadrant':
            xlist = quadrant_list
            mean_rates = self.df_mean_per_quadrant.MeanRate
            std_rates = self.df_mean_per_quadrant.StdRate

        self.xdata = np.array(xlist)
        self.mean_rates = np.array(mean_rates)
        self.std_rates = np.array(std_rates)

        a_fit,cov=curve_fit(self.linearFunc,self.xdata,self.mean_rates,sigma=std_rates,absolute_sigma=True)

        inter = a_fit[0]
        slope = a_fit[1]
        d_inter = np.sqrt(cov[0][0])
        d_slope = np.sqrt(cov[1][1])

        self.yfit = inter + slope * self.xdata

        print('d = ', inter, ', k = ', slope)


    def create_plot(self, whatplot):

        plt.errorbar(self.xdata,self.mean_rates,yerr=self.std_rates,color='tab:blue',fmt = 'o',label='Simulation Data', linewidth = 2)
        plt.plot(self.xdata,self.yfit,label='Linear Fit',color='tab:orange', linewidth = 2)

        plt.xlabel(whatplot)
        plt.ylabel('Hit Rate (MHz/cm²)')

        plt.legend()
        plt.tight_layout()

        if self.sensor == '':
            sensor_name = 'AllSensors'
        else:
            sensor_name = 'Sensor'+str(self.sensor)

        plt.savefig('SimulatedHitRatePer'+str(whatplot)+sensor_name+'.pdf')
        plt.show()


    def run(self):
        self.get_raw_data()
        self.data_for_sensor()
        self.get_closest_sensor()
        self.get_farthest_sensor()
        self.lq_data()
        self.mean_per_layer()
        self.mean_per_quadrant()
        self.mean_overall()
        self.perform_fit('Layer')
        self.create_plot('Layer')
        self.perform_fit('Quadrant')
        self.create_plot('Quadrant')


data = ExtractSimData(sensor_id = 979) # 739 1935
data.run()

print(data.df_lq)
print(data.df_mean_per_layer)
print(data.df_mean_per_quadrant)
print(data.df_mean_overall)








    





