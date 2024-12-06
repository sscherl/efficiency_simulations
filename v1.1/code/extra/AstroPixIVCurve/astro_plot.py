import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import csv

df = pd.DataFrame(pd.read_csv('ui_test_20230323-151008.csv'))
df_lowres = pd.DataFrame(pd.read_csv('ui_test_lowRes_20230323-155750.csv'))

df.U = df.U.apply(lambda x: x*-1)
df.I = df.I.apply(lambda x: x*-1000000)

df_lowres.U = df_lowres.U.apply(lambda x: x*-1)
df_lowres.I = df_lowres.I.apply(lambda x: x*-1000000)

df_sub = df
df_sub.I = df_sub.I - df_lowres.I

plt.plot(df_sub.U, df_sub.I)
plt.plot(df.U, df.I)

plt.title('AstroPix2 Breakdown Voltage HiRes Wafer')
plt.xlabel('High Voltage (- V)')
plt.ylabel('Current (- uA)')
plt.yscale('log')


plt.tight_layout()


plt.savefig('IVcurve_ui_test_20230323-151008.pdf')
plt.show()
