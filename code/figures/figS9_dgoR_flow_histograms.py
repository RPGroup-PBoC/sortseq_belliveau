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

df = pd.read_csv('../flow/20160325_flow_hist.csv', comment='#')

# histogram # 1 details
date = 20160325
promoter = 'na'
strain = 'MG1655'
region = 'na'
bin = 'na'
sequence = 'na'
media = 'M9glucose'
condition = 'na'

df_hist1 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 2 details
date = 20160325
promoter = 'dgoR'
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

# histogram # 3 details
date = 20160325
promoter = 'dgoR'
strain = 'MG1655'
region = 'na'
bin = 'na'
sequence = 'wild-type'
media = 'M9galactose'
condition = 'na'

df_hist3 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram #4 details
date = 20160325
promoter = 'dgoR'
strain = 'MG1655'
region = 'na'
bin = 'na'
sequence = 'wild-type'
media = 'M9galactonate'
condition = 'na'

df_hist4 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]



#===============================================================================
# plot the data
#===============================================================================

fig = plt.figure(figsize = (4,2))

plt.semilogx(df_hist2['fluorescence'], df_hist2['fraction'], linewidth=2,
            label=r'$dgoR$, glucose')
plt.semilogx(df_hist3['fluorescence'], df_hist3['fraction'], linewidth=2,
            label=r'$dgoR$, galactose')
plt.semilogx(df_hist4['fluorescence'], df_hist4['fraction'], linewidth=2,
            label=r'$dgoR$, galactonate')
plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=2,
            color = '#CC79A7', ls = '-.', label='autofluorescence')
plt.yticks([])
plt.xlim(10,1E5)
plt.xlabel('Fluorescence (a.u.)')
plt.ylabel('Frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig.savefig(output + 'figS9_dgoR_flow_histograms_carbonsource.pdf',
bbox_extra_artists=(lgd,), bbox_inches='tight')
