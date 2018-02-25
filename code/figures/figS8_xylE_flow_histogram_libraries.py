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

df = pd.read_csv('../flow/20160519_flow_hist_corrected.csv', comment='#')

# histogram # 1 details
date = 20160519
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
date = 20160519
promoter = 'xylE'
strain = 'MG1655'
region = 'mut1'
bin = 'na'
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist2 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]


# histogram # 3 details
date = 20160519
promoter = 'xylE'
strain = 'MG1655'
region = 'mut2'
bin = 'na'
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist3 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 3 details
date = 20160519
promoter = 'xylE'
strain = 'MG1655'
region = 'mut3'
bin = 'na'
sequence = 'library'
media = 'M9glucose'
condition = 'na'

df_hist4 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]


#===============================================================================
# Read the data for growth in xlyose
#===============================================================================

df = pd.read_csv('../flow/20160527_flow_hist.csv', comment='#')

# histogram # 5 details
date = 20160527
promoter = 'xylE'
strain = 'MG1655'
region = 'mut1'
bin = 'na'
sequence = 'library'
media = 'M9xylose'
condition = 'na'

df_hist5 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]


# histogram # 6 details
date = 20160527
promoter = 'xylE'
strain = 'MG1655'
region = 'mut2'
bin = 'na'
sequence = 'library'
media = 'M9xylose'
condition = 'na'

df_hist6 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]

# histogram # 7 details
date = 20160527
promoter = 'xylE'
strain = 'MG1655'
region = 'mut3'
bin = 'na'
sequence = 'library'
media = 'M9xylose'
condition = 'na'

df_hist7 = df[(df.date == date) & (df.promoter == promoter) \
                & (df.strain == strain) & (df.region == region) \
                & (df.bin == bin) & (df.sequence == sequence) \
                & (df.media == media) & (df.condition == condition)]


#===============================================================================
# plot the data
#===============================================================================

fig = plt.figure(figsize = (4,2))

plt.semilogx(df_hist2['fluorescence'], df_hist2['fraction'], linewidth=2,
            color='#0072B2', label='library 1,  glucose')
plt.semilogx(df_hist3['fluorescence'], df_hist3['fraction'], linewidth=2,
            color='#0072B2', ls = ':', label='library 2,  glucose')
plt.semilogx(df_hist4['fluorescence'], df_hist4['fraction'], linewidth=2,
            color='#0072B2', ls = '-.', label='library 3,  glucose')
plt.semilogx(df_hist5['fluorescence'], df_hist5['fraction'], linewidth=2,
            color ='#D55E00' ,label='library 1,  xylose')
plt.semilogx(df_hist6['fluorescence'], df_hist6['fraction'], linewidth=2,
            color ='#D55E00', ls = ':', label='library 2,  xylose')
plt.semilogx(df_hist7['fluorescence'], df_hist7['fraction'], linewidth=2,
            color ='#D55E00', ls = '-.', label='library 3,  xylose')
plt.semilogx(df_hist1['fluorescence'], df_hist1['fraction'], linewidth=2,
            color = '#009E73', ls = '-.', label='autofluorescence')
plt.yticks([])
plt.xlim(10,1E5)
plt.xlabel('fluorescence (a.u.)')
plt.ylabel('frequency')

lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.tight_layout()

fig.savefig(output + 'figS8_xylE_flow_histograms_libraries.pdf',
bbox_extra_artists=(lgd,), bbox_inches='tight')
