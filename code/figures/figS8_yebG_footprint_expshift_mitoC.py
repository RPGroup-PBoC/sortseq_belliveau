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
fname = '20170717_yebG_MG1655_M9glucose_mitoC_mut2_4bins_summary.csv'
df = pd.read_csv(data_dir + fname)

seqLength = 60
#------------------------------------------------------------------------------#
# make shared axis plots of information footprint and delta bin shift
#------------------------------------------------------------------------------#
colours = ["#348ABD", "#A60628","#87181C","#E8DCD2"]

fig, (ax1, ax2) = plt.subplots(2,1,figsize=(2,2), gridspec_kw = {'height_ratios':[1, 1]})
# make array for x positions - note that this will be relative to TSS
ind = np.arange(seqLength) - 74


# footprint
ax1.bar(ind, df[df['MI']!=0].sort_values(by='position').MI.values, alpha=1, width=0.8, \
        linewidth = 0, color=colours[2])
ax1.set_xlim(ind[0],ind[-1])
ax1.grid(b=False)
ax1.tick_params(labelsize=8)
# # Hide the right and top spines
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
#
# # Only show ticks on the left and bottom spines
ax1.yaxis.set_ticks_position('left')
ax1.xaxis.set_ticks_position('bottom')

# Annotate y axis
ylim = ax1.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax1.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax1.set_yticklabels([' ',' '], fontname='Arial')
ax1.set_ylabel('mutual\ninformation', fontsize=8)

# Make grid lines
xticks = [x for x in ind if x%20 == 0]
for n, x in enumerate(xticks):
        ax1.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)


# delta bin
ax2.bar(ind, df[df['MI']!=0].sort_values(by='position').expshift.values, width=0.8, \
        linewidth = 0, color = 'k')
ax2.set_xlim(ind[0],ind[-1])
ax2.set_xlabel('position (relative to ${\it yebG}$ gene)', fontsize=8)
ax2.grid(b=False)
ax2.set_facecolor('white')
ax2.tick_params(labelsize=8)

# Hide the right and top spines
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax2.yaxis.set_ticks_position('left')
ax2.xaxis.set_ticks_position('bottom')


# Annotate y axis
ylim = ax2.get_ylim()
#ax.set_yticks([])
yticks = np.linspace(ylim[0],ylim[1],5)[[1,3]]
ax2.set_yticks(yticks)
# fontProperties = {'family':'sans-serif', 'sans-serif':['Helvetica'], 'weight' : 'normal', 'size' : 14}
ax2.set_yticklabels(['-','+'], fontname='Arial')
ax2.set_ylabel('expression\nshift', fontsize=8)

# Make grid lines
xticks = [x for x in ind if x%20 == 0]
for n, x in enumerate(xticks):
        ax2.axvline(x, color='lightgray', linewidth=0.36, zorder=-1)

plt.tight_layout()

figname_out = 'figS8_yebG_footprint_expshift_mitoC.pdf'
fig.savefig(output + figname_out, format='pdf')
