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
output = 'output_figs/'

#------------------------------------------------------------------------------#
# Load in the summary csv file
#------------------------------------------------------------------------------#

data_dir = '../sortseq/20170717_yebG/'
fname = '20170717_yebG_MG1655_M9glucose_na_mut2_4bins_summary.csv'
df = pd.read_csv(data_dir + fname)

seqLength = 60
#------------------------------------------------------------------------------#
# make shared axis plots of information footprint and delta bin shift
#------------------------------------------------------------------------------#
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

# fig, (ax1, ax2) = plt.subplots(2,1,figsize=(5,5), gridspec_kw = {'height_ratios':[1.75, 1]})
fig, (ax1, ax2) = plt.subplots(2,1,figsize=(2,2), gridspec_kw = {'height_ratios':[1, 1]})

# make array for x positions - note that this will be relative to TSS
ind = np.arange(seqLength) - 74 # assumed to be at position 6 relative to -10, so approx


# footprint
ax1.bar(ind, df[df['MI']!=0].sort_values(by='position').MI.values[::-1], alpha=1, width=0.8, \
        linewidth = 0, color=colours[2])
ax1.set_xlim(ind[0],ind[-1])
ax1.set_ylabel('mutual\ninformation')
ax1.grid(b=False)

# # Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)


# delta bin
ax2.bar(ind, df[df['MI']!=0].sort_values(by='position').expshift.values[::-1], width=0.8, \
        linewidth = 0, color = 'k')
ax2.set_xlim(ind[0],ind[-1])
ax2.set_ylabel('fluorescence\nbin shift')
ax2.set_xlabel('position')
ax2.grid(b=False)
ax2.set_facecolor('white')

# Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')

# make these tick labels invisible
#plt.setp(ax1.get_xticklabels(), fontsize=6)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.set_xlabel('position\n(relative to ${\it yebG}$ gene)')

fig.subplots_adjust(hspace=-1)
plt.tight_layout()

figname_out = 'figS8_yebG_footprint_expshift.pdf'
fig.savefig(output + figname_out, format='pdf')
