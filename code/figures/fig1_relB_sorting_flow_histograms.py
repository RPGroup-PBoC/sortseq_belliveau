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

df = pd.read_csv('../flow/20150209_flow_hist.csv', comment='#')

# histogram # 1 details
date = 20150209
promoter = 'relB'
strain = 'MG1655'
region = 'na'
bin = 'na'
sequence = 'wild-type'
media = 'M9glucose'
condition = 'na'

df_wt_hist1 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 1 details
date = 20150209
promoter = 'relB'
strain = 'MG1655'
region = 'mut1'
bin = 'na'
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_lib_hist1 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

df = pd.read_csv('../flow/20150210_flow_hist.csv', comment='#')

# histogram # 1 details
date = 20150210
promoter = 'relB'
strain = 'MG1655'
region = 'mut1'
bin = 1
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist1 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 1 details
date = 20150210
promoter = 'relB'
strain = 'MG1655'
region = 'mut1'
bin = 2
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist2 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 1 details
date = 20150210
promoter = 'relB'
strain = 'MG1655'
region = 'mut1'
bin = 3
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist3 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 1 details
date = 20150210
promoter = 'relB'
strain = 'MG1655'
region = 'mut1'
bin = 4
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist4 = df[(df.date == date) & (df.promoter == promoter) \
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

fig1 = plt.figure(figsize = (5,4))

plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=2,
            label=r'$relB$ library', alpha = 1, color = '#EFF6EF')
plt.fill_between(df_hist1['fluorescence'], df_hist1['fraction'], 0,  color = '#EFF6EF')
plt.semilogx(df_hist2['fluorescence'], df_hist2['fraction'], linewidth=2,
            label=r'$relB$ wild-type', alpha = 1, color = '#C6E1CD')
plt.fill_between(df_hist2['fluorescence'], df_hist2['fraction'], 0,  color = '#C6E1CD')
plt.semilogx(df_hist3['fluorescence'], df_hist3['fraction'], linewidth=2,
            label=r'$relB$ wild-type', alpha = 1, color = '#A0CCA9')
plt.fill_between(df_hist3['fluorescence'], df_hist3['fraction'], 0,  color = '#A0CCA9')
plt.semilogx(df_hist4['fluorescence'], df_hist4['fraction'], linewidth=2,
            label=r'$relB$ wild-type', alpha = 1, color = '#67AC74')
plt.fill_between(df_hist4['fluorescence'], df_hist4['fraction'], 0,  color = '#67AC74')

plt.semilogx(df_lib_hist1['fluorescence'], df_lib_hist1['fraction'], linewidth=2,
            label=r'$relB$ lib', alpha = 0.4, color = 'k', ls='--')

# plt.yticks([])
plt.xlim(1,1E5)
plt.ylim(bottom=0)
plt.xlabel('Fluorescence (a.u.)')
plt.ylabel('Frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig1.savefig(output + 'fig1_relB_histograms.pdf', bbox_extra_artists=(lgd,), bbox_inches='tight')
