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
utils.set_plotting_style4()

#===============================================================================
# Set output directory
#===============================================================================
output = 'input_data/'
#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#

data_dir_1 = '../sortseq/20160707_purT_xylE_dgoR/'
fname_1 = '20160710_dgoR_MG1655_M9glucose_na_mut1_4bins_summary.csv'
df1 = pd.read_csv(data_dir_1 + fname_1)
df1['position'] = df1['position']

data_dir_2 = '../sortseq/20160707_purT_xylE_dgoR/'
fname_2 = '20160710_dgoR_MG1655_M9glucose_na_mut2_4bins_summary.csv'
df2 = pd.read_csv(data_dir_2 + fname_2)

data_dir_3 = '../sortseq/20160707_purT_xylE_dgoR/'
fname_3 = '20160710_dgoR_MG1655_M9glucose_na_mut3_4bins_summary.csv'
df3 = pd.read_csv(data_dir_3 + fname_3)

df = df1.append(df2)
df = df.append(df3)
# df = df[df.position>=9]
print(len(df))

seqLength = 147

#------------------------------------------------------------------------------#
# make plot of expression shifts
#------------------------------------------------------------------------------#
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

# fig, (ax1, ax1) = plt.subplots(2,1,figsize=(12.5,6), gridspec_kw = {'height_ratios':[1.75, 1]})
# make array for x positions - note that this will be relative to TSS
ind = np.arange(seqLength) - 161 # assumed to be at position 6 relative to -10, so approx

fig = plt.figure(figsize=utils.cm2inch((9,3)))
ax1 = fig.add_subplot(111)

# delta bin
ax1.bar(ind, df[df['expshift']!=0].groupby(['position'])['expshift'].mean(), width=0.8, \
        linewidth = 0, color = 'k')
ax1.set_xlim(ind[0],ind[-1])
ax1.set_xlabel('position (relative to $dgoR$ gene)', fontsize=8)
ax1.grid(b=False)
ax1.set_facecolor('white')

# Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

# Annotate y axis
ylim = ax1.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax1.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax1.set_yticklabels(['-','+'], fontname='Arial')
ax1.set_ylabel('expression\nshift')

# Make grid lines
xticks = [x for x in ind if x%20 == 0]
for n, x in enumerate(xticks):
        ax1.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)

plt.tight_layout()

figname_out = 'fig7_dgoR_expshfit_glucose.pdf'
fig.savefig(output + figname_out, format='pdf')


# lets also save the values to a txt file; I want to compare between conditions

df[df['expshift']!=0].groupby(['position'])[['expshift']].mean().to_csv(output + 'fig7_dgoR_expshfit_glucose.csv')
