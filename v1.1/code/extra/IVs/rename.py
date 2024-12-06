from os import listdir
from os.path import isfile, join
import os

mypath = '/Users/sigridscherl/Documents/GitRepos/mightypix_measurements/my_gecco-daq/output/iv_scan/'

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

# print(onlyfiles)

for f in onlyfiles:
    filename = str(f)
    if filename[0:4] != '2024':
        day = filename[0:2]
        month = filename[3:5]
        year = filename[6:10]
        rest = filename[10:]

        old = mypath+str(day)+'-'+str(month)+'-'+str(year) +rest
        new = mypath+str(year)+'-'+str(month)+'-'+str(day) +rest
        try:
            os.rename(old, new)
        except: pass
        # new_name = mypath + 
