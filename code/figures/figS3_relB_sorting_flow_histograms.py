import os
import glob
import pickle
import re

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys
sys.path.insert(0, '../')
import NB_sortseq_utils as utils

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from IPython.core.pylabtools import figsize

# Seaborn, useful for graphics
import seaborn as sns

sns.set_palette("deep", color_codes=True)
utils.set_plotting_style1()

#===============================================================================
# Set output directory
#===============================================================================
output = 'output_figs/'

#===============================================================================
# Read the data
#===============================================================================

df = pd.read_csv('../flow/20160311_flow_hist.csv', comment='#')

# histogram # 1 details
date = 20160311
promoter = 'relB'
strain = 'MG1655'
region = 'mut1'
bin = 'na'
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist1 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 1 details
date = 20160311
promoter = 'relB'
strain = 'MG1655'
region = 'na'
bin = 'na'
sequence = 'wild-type'
media = 'M9glucose'
condition = 'na'

df_hist2 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

#===============================================================================
# identify gate windows (approximate)
#===============================================================================
df_hist1['cum_sum'] = df_hist1.fraction.cumsum()

#===============================================================================
# plot the data
#===============================================================================
palette = sns.color_palette()

fig1 = plt.figure(figsize = (4,2))

plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=2,
            label=r'$relB$ library')
plt.semilogx(df_hist2['fluorescence'], df_hist2['fraction'], linewidth=2,
            label=r'$relB$ wild-type', alpha = 0.4, color = 'k')

plt.yticks([])
plt.xlim(1,1E5)
plt.xlabel('Fluorescence (a.u.)')
plt.ylabel('Frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig1.savefig(output + 'fig_SI_relBE_binning_WT_lib.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')

#===============================================================================
fig2 = plt.figure(figsize = (4,2))

plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=1, color = palette[0])
x_gate = [[df_hist1[df_hist1.cum_sum<=df_hist1.cum_sum.min()].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.15].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.2833].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.4333].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.56663].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.71663].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.85].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=df_hist1.cum_sum.max()].fluorescence.values[-1]]]
x = df_hist2['fluorescence']
for i in range(0,4):
    plt.fill_between(x,df_hist1['fraction'],color = palette[0],where =(x > x_gate[i][0]) & (x < x_gate[i][1]))

plt.yticks([])
plt.xlim(1,1E5)
plt.xlabel('Fluorescence (a.u.)')
plt.ylabel('Frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig2.savefig(output + 'figS3_relBE_sorting_hist_4bin_15percent.pdf',
            bbox_extra_artists=(lgd,), bbox_inches='tight')

#===============================================================================
fig3 = plt.figure(figsize = (4,2))

plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=1, color = palette[0])
x_gate = [[df_hist1[df_hist1.cum_sum<=df_hist1.cum_sum.min()].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.22].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.26].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.48].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.52].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.74].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.78].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=df_hist1.cum_sum.max()].fluorescence.values[-1]]]
x = df_hist2['fluorescence']
for i in range(0,4):
    plt.fill_between(x,df_hist1['fraction'],color = palette[0],where =(x > x_gate[i][0]) & (x < x_gate[i][1]))

plt.yticks([])
plt.xlim(1,1E5)
plt.xlabel('Fluorescence (a.u.)')
plt.ylabel('Frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig3.savefig(output + 'figS3_relBE_sorting_hist_4bin_22percent.pdf',
            bbox_extra_artists=(lgd,), bbox_inches='tight')

#===============================================================================
fig4 = plt.figure(figsize = (4,2))

plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=1, color = palette[0])
x_gate = [[df_hist1[df_hist1.cum_sum<=df_hist1.cum_sum.min()].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.11].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.12714].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.23714].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.25428].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.36428].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.38142].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.49142].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.50856].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.61856].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.6357].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.7457].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.76284].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=0.87284].fluorescence.values[-1]],
            [df_hist1[df_hist1.cum_sum>=0.88].fluorescence.values[0],
            df_hist1[df_hist1.cum_sum<=df_hist1.cum_sum.max()].fluorescence.values[-1]]]
x = df_hist2['fluorescence']
for i in range(0,8):
    plt.fill_between(x,df_hist1['fraction'],color = palette[0],where =(x > x_gate[i][0]) & (x < x_gate[i][1]))

plt.yticks([])
plt.xlim(1,1E5)
plt.xlabel('Fluorescence (a.u.)')
plt.ylabel('Frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig4.savefig(output + 'figS3_relBE_sorting_hist_8bin_10percent.pdf',
                bbox_extra_artists=(lgd,), bbox_inches='tight')
